//! GPU solver for the 2D grid (feature `gpu`, native only): matrix-free Jacobi-preconditioned
//! conjugate gradients in f64 through wgpu (Vulkan). The whole CG loop runs on the GPU; the CPU
//! only reads the residual every `CHECK_EVERY` steps. Same DOF layout, element matrix and
//! fixed-DOF rule as `mg2d::Level::apply`, so results compare one to one with the CPU solver.
#![allow(clippy::indexing_slicing)] // sizes are fixed at construction

use wgpu::util::DeviceExt;

const WG: u32 = 256;
const CHECK_EVERY: usize = 25;

#[repr(C)]
#[derive(Clone, Copy, bytemuck::Pod, bytemuck::Zeroable)]
struct Params {
    nx: u32,
    ny: u32,
    ndof: u32,
    npart: u32,
}


/// Opens the GPU: the adapter whose name contains `GDE_GPU` if set, else a discrete GPU with
/// 64-bit float shaders (the AMD card preferred here).
fn open_device() -> Result<(wgpu::Device, wgpu::Queue, String), String> {
        let instance = wgpu::Instance::new(wgpu::InstanceDescriptor { backends: wgpu::Backends::VULKAN, ..wgpu::InstanceDescriptor::new_without_display_handle() });
        let want = std::env::var("GDE_GPU").ok();
        let adapters = pollster::block_on(instance.enumerate_adapters(wgpu::Backends::VULKAN));
        let adapter = adapters
            .into_iter()
            .filter(|a| a.features().contains(wgpu::Features::SHADER_F64))
            .filter(|a| match &want {
                Some(w) => a.get_info().name.contains(w.as_str()),
                None => a.get_info().device_type == wgpu::DeviceType::DiscreteGpu,
            })
            .max_by_key(|a| (a.get_info().vendor == 0x1002) as u8) // prefer the AMD card here
            .ok_or("no GPU with 64-bit float shaders found")?;
        let adapter_name = adapter.get_info().name;
        let mut limits = wgpu::Limits::default();
        let al = adapter.limits();
        limits.max_storage_buffers_per_shader_stage = al.max_storage_buffers_per_shader_stage.min(32);
        limits.max_storage_buffer_binding_size = al.max_storage_buffer_binding_size;
        limits.max_buffer_size = al.max_buffer_size;
        let (device, queue) = pollster::block_on(adapter.request_device(&wgpu::DeviceDescriptor {
            label: Some("gde"),
            required_features: wgpu::Features::SHADER_F64 | (adapter.features() & wgpu::Features::TIMESTAMP_QUERY),
            required_limits: limits,
            ..Default::default()
        }))
        .map_err(|e| format!("GPU device: {e}"))?;
        Ok((device, queue, adapter_name))
}

pub struct GpuPcg {
    device: wgpu::Device,
    queue: wgpu::Queue,
    pub adapter_name: String,
    nx: usize,
    ny: usize,
    ndof: usize,
    npart: u32,
    pipes: Pipes,
    layout: wgpu::BindGroupLayout,
    params: wgpu::Buffer,
    ke: wgpu::Buffer,
    fixed: wgpu::Buffer,
}

struct Pipes {
    apply: wgpu::ComputePipeline,
    copy_u_p: wgpu::ComputePipeline,
    init_r: wgpu::ComputePipeline,
    reduce_rz: wgpu::ComputePipeline,
    reduce_pq: wgpu::ComputePipeline,
    final_init: wgpu::ComputePipeline,
    final_alpha: wgpu::ComputePipeline,
    update_ur: wgpu::ComputePipeline,
    final_beta: wgpu::ComputePipeline,
    update_p: wgpu::ComputePipeline,
}

/// The result of one solve.
pub struct GpuSolve {
    pub u: Vec<f64>,
    pub iterations: usize,
    pub converged: bool,
    /// Multigrid levels (1 for the Jacobi solver).
    pub levels: usize,
}

fn f64_bytes(v: &[f64]) -> &[u8] {
    bytemuck::cast_slice(v)
}

impl GpuPcg {
    /// Opens the GPU (the adapter whose name contains `GDE_GPU` if set, else the first discrete
    /// GPU with f64) and uploads what does not change between solves.
    pub fn new(nx: usize, ny: usize, ke: &[[f64; 8]; 8], fixed: &[bool]) -> Result<Self, String> {
        let (device, queue, adapter_name) = open_device()?;
        let module = device.create_shader_module(wgpu::ShaderModuleDescriptor { label: Some("pcg"), source: wgpu::ShaderSource::Wgsl(include_str!("gpu_pcg.wgsl").into()) });
        let entry = |b: u32, ty: wgpu::BufferBindingType| wgpu::BindGroupLayoutEntry {
            binding: b,
            visibility: wgpu::ShaderStages::COMPUTE,
            ty: wgpu::BindingType::Buffer { ty, has_dynamic_offset: false, min_binding_size: None },
            count: None,
        };
        let ro = wgpu::BufferBindingType::Storage { read_only: true };
        let rw = wgpu::BufferBindingType::Storage { read_only: false };
        let mut entries = vec![entry(0, wgpu::BufferBindingType::Uniform)];
        for b in 1..=5 {
            entries.push(entry(b, ro));
        }
        for b in 6..=12 {
            entries.push(entry(b, rw));
        }
        let layout = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor { label: Some("pcg"), entries: &entries });
        let pl = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor { label: Some("pcg"), bind_group_layouts: &[Some(&layout)], immediate_size: 0 });
        let pipe = |name: &str| {
            device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor { label: Some(name), layout: Some(&pl), module: &module, entry_point: Some(name), compilation_options: Default::default(), cache: None })
        };
        let pipes = Pipes {
            apply: pipe("apply_pq"),
            copy_u_p: pipe("copy_u_p"),
            init_r: pipe("init_r"),
            reduce_rz: pipe("reduce_rz"),
            reduce_pq: pipe("reduce_pq"),
            final_init: pipe("final_init"),
            final_alpha: pipe("final_alpha"),
            update_ur: pipe("update_ur"),
            final_beta: pipe("final_beta"),
            update_p: pipe("update_p"),
        };
        let ndof = 2 * (nx + 1) * (ny + 1);
        if fixed.len() != ndof {
            return Err("fixed mask has the wrong length".into());
        }
        let npart = (ndof as u32).div_ceil(WG);
        let params = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("params"),
            contents: bytemuck::bytes_of(&Params { nx: nx as u32, ny: ny as u32, ndof: ndof as u32, npart }),
            usage: wgpu::BufferUsages::UNIFORM,
        });
        let flat: Vec<f64> = ke.iter().flatten().copied().collect();
        let ke = device.create_buffer_init(&wgpu::util::BufferInitDescriptor { label: Some("ke"), contents: f64_bytes(&flat), usage: wgpu::BufferUsages::STORAGE });
        let mask: Vec<u32> = fixed.iter().map(|b| *b as u32).collect();
        let fixed = device.create_buffer_init(&wgpu::util::BufferInitDescriptor { label: Some("fixed"), contents: bytemuck::cast_slice(&mask), usage: wgpu::BufferUsages::STORAGE });
        Ok(Self { device, queue, adapter_name, nx, ny, ndof, npart, pipes, layout, params, ke, fixed })
    }

    /// Solves K u = f with K = sum over cells of modulus[e] * ke, warm-started from `start`.
    /// `inv_diag` is 1 / diag(K) (as `mg2d::Level::diagonal`).
    pub fn solve(&self, modulus: &[f64], inv_diag: &[f64], f: &[f64], start: Option<&[f64]>, tol: f64, max: usize) -> Result<GpuSolve, String> {
        let n = self.ndof;
        if modulus.len() != self.nx * self.ny || inv_diag.len() != n || f.len() != n {
            return Err("GPU solve: input sizes do not match the grid".into());
        }
        let d = &self.device;
        let st = wgpu::BufferUsages::STORAGE;
        let init = |label: &str, data: &[f64], extra: wgpu::BufferUsages| d.create_buffer_init(&wgpu::util::BufferInitDescriptor { label: Some(label), contents: f64_bytes(data), usage: st | extra });
        let zero = vec![0.0f64; n];
        let modulus_b = init("modulus", modulus, wgpu::BufferUsages::empty());
        let invd = init("invd", inv_diag, wgpu::BufferUsages::empty());
        let fb = init("f", f, wgpu::BufferUsages::empty());
        let u = init("u", start.filter(|s| s.len() == n).unwrap_or(&zero), wgpu::BufferUsages::COPY_SRC);
        let r = init("r", &zero, wgpu::BufferUsages::empty());
        let z = init("z", &zero, wgpu::BufferUsages::empty());
        let p = init("p", &zero, wgpu::BufferUsages::empty());
        let q = init("q", &zero, wgpu::BufferUsages::empty());
        let part = init("part", &vec![0.0; 2 * self.npart as usize], wgpu::BufferUsages::empty());
        let sc = init("sc", &[0.0; 8], wgpu::BufferUsages::COPY_SRC);
        let bind = d.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("pcg"),
            layout: &self.layout,
            entries: &[&self.params, &self.ke, &modulus_b, &self.fixed, &invd, &fb, &u, &r, &z, &p, &q, &part, &sc]
                .iter()
                .enumerate()
                .map(|(i, b)| wgpu::BindGroupEntry { binding: i as u32, resource: b.as_entire_binding() })
                .collect::<Vec<_>>(),
        });
        let nodes = ((self.nx + 1) * (self.ny + 1)) as u32;
        let groups = |count: u32| count.div_ceil(WG);
        let run = |enc: &mut wgpu::CommandEncoder, pipe: &wgpu::ComputePipeline, x: u32| {
            let mut pass = enc.begin_compute_pass(&wgpu::ComputePassDescriptor { label: None, timestamp_writes: None });
            pass.set_pipeline(pipe);
            pass.set_bind_group(0, &bind, &[]);
            pass.dispatch_workgroups(x, 1, 1);
        };
        let staging = d.create_buffer(&wgpu::BufferDescriptor { label: Some("staging"), size: 64, usage: wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST, mapped_at_creation: false });
        let read_sc = |enc: wgpu::CommandEncoder| -> Result<[f64; 8], String> {
            let mut enc = enc;
            enc.copy_buffer_to_buffer(&sc, 0, &staging, 0, 64);
            self.queue.submit([enc.finish()]);
            let slice = staging.slice(..);
            slice.map_async(wgpu::MapMode::Read, |_| {});
            d.poll(wgpu::PollType::wait_indefinitely()).map_err(|e| format!("GPU poll: {e}"))?;
            let mut out = [0.0; 8];
            {
                let view = slice.get_mapped_range().map_err(|e| format!("GPU map: {e:?}"))?;
                out.copy_from_slice(bytemuck::cast_slice(&view));
            }
            staging.unmap();
            Ok(out)
        };
        let fnorm = f.iter().map(|v| v * v).sum::<f64>().sqrt().max(f64::MIN_POSITIVE);
        let (dofs, npart) = (groups(n as u32), self.npart);
        let mut enc = d.create_command_encoder(&Default::default());
        run(&mut enc, &self.pipes.copy_u_p, dofs);
        run(&mut enc, &self.pipes.apply, groups(nodes));
        run(&mut enc, &self.pipes.init_r, dofs);
        run(&mut enc, &self.pipes.reduce_rz, npart);
        run(&mut enc, &self.pipes.final_init, 1);
        let mut s = read_sc(enc)?;
        let mut it = 0;
        let mut converged = s[3].sqrt() <= tol * fnorm;
        while !converged && it < max && s[4] == 0.0 {
            let mut enc = d.create_command_encoder(&Default::default());
            let steps = CHECK_EVERY.min(max - it);
            for _ in 0..steps {
                run(&mut enc, &self.pipes.apply, groups(nodes));
                run(&mut enc, &self.pipes.reduce_pq, npart);
                run(&mut enc, &self.pipes.final_alpha, 1);
                run(&mut enc, &self.pipes.update_ur, dofs);
                run(&mut enc, &self.pipes.reduce_rz, npart);
                run(&mut enc, &self.pipes.final_beta, 1);
                run(&mut enc, &self.pipes.update_p, dofs);
            }
            it += steps;
            s = read_sc(enc)?;
            converged = s[3].sqrt() <= tol * fnorm;
        }
        // Read u back.
        let out = d.create_buffer(&wgpu::BufferDescriptor { label: Some("u out"), size: (n * 8) as u64, usage: wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST, mapped_at_creation: false });
        let mut enc = d.create_command_encoder(&Default::default());
        enc.copy_buffer_to_buffer(&u, 0, &out, 0, (n * 8) as u64);
        self.queue.submit([enc.finish()]);
        let slice = out.slice(..);
        slice.map_async(wgpu::MapMode::Read, |_| {});
        d.poll(wgpu::PollType::wait_indefinitely()).map_err(|e| format!("GPU poll: {e}"))?;
        let u_host: Vec<f64> = {
            let view = slice.get_mapped_range().map_err(|e| format!("GPU map: {e:?}"))?;
            bytemuck::cast_slice(&view).to_vec()
        };
        out.unmap();
        Ok(GpuSolve { u: u_host, iterations: it, converged, levels: 1 })
    }
}

// ---------------------------------------------------------------------------------------------
// Multigrid-preconditioned CG: the V-cycle's levels on the GPU, the coarsest solve on the CPU.

use crate::mg2d::{Cells, Hierarchy, Level, Settings};

#[repr(C)]
#[derive(Clone, Copy, bytemuck::Pod, bytemuck::Zeroable)]
struct LevelParams {
    nx: u32,
    ny: u32,
    ndof: u32,
    scaled: u32,
    npart: u32,
    omega: f32,
    pad0: u32,
    pad1: u32,
}

#[repr(C)]
#[derive(Clone, Copy, bytemuck::Pod, bytemuck::Zeroable)]
struct CoarseParams {
    nx: u32,
    ny: u32,
    ndof: u32,
    pad: u32,
}

struct MgPipes {
    apply_t: wgpu::ComputePipeline,
    apply_q: wgpu::ComputePipeline,
    zero_z: wgpu::ComputePipeline,
    smooth: wgpu::ComputePipeline,
    resid: wgpu::ComputePipeline,
    restrict_r: wgpu::ComputePipeline,
    prolong_add: wgpu::ComputePipeline,
    copy_u_p: wgpu::ComputePipeline,
    init_r: wgpu::ComputePipeline,
    set_p_z: wgpu::ComputePipeline,
    reduce_rz: wgpu::ComputePipeline,
    reduce_pq: wgpu::ComputePipeline,
    final_init: wgpu::ComputePipeline,
    final_alpha: wgpu::ComputePipeline,
    update_ur: wgpu::ComputePipeline,
    final_beta: wgpu::ComputePipeline,
    update_p: wgpu::ComputePipeline,
}

/// Multigrid-preconditioned conjugate gradients with the V-cycle on the GPU.
pub struct GpuMg {
    device: wgpu::Device,
    queue: wgpu::Queue,
    pub adapter_name: String,
    pipes: MgPipes,
    g0: wgpu::BindGroupLayout,
    g1: wgpu::BindGroupLayout,
    g2: wgpu::BindGroupLayout,
}

/// Timing split of one GPU multigrid solve, seconds.
#[derive(Default, Clone, Copy, Debug)]
pub struct MgTimes {
    pub hierarchy: f64,
    pub upload: f64,
    pub iterate: f64,
    pub coarse_cpu: f64,
    /// GPU time of the last step's two halves (up + p update, then q/alpha/update + down), seconds; 0 without timestamp support.
    pub gpu_up: f64,
    pub gpu_down: f64,
}

struct GpuLevel {
    nx: usize,
    ny: usize,
    ndof: usize,
    r: wgpu::Buffer,
    z: wgpu::Buffer,
    g0: wgpu::BindGroup,
    g1: wgpu::BindGroup,
}

impl GpuMg {
    pub fn new() -> Result<Self, String> {
        let (device, queue, adapter_name) = open_device()?;
        let module = device.create_shader_module(wgpu::ShaderModuleDescriptor { label: Some("mg"), source: wgpu::ShaderSource::Wgsl(include_str!("gpu_mg.wgsl").into()) });
        let entry = |b: u32, ty: wgpu::BufferBindingType| wgpu::BindGroupLayoutEntry {
            binding: b,
            visibility: wgpu::ShaderStages::COMPUTE,
            ty: wgpu::BindingType::Buffer { ty, has_dynamic_offset: false, min_binding_size: None },
            count: None,
        };
        let uni = wgpu::BufferBindingType::Uniform;
        let ro = wgpu::BufferBindingType::Storage { read_only: true };
        let rw = wgpu::BufferBindingType::Storage { read_only: false };
        let g0 = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
            label: Some("level"),
            entries: &[entry(0, uni), entry(1, ro), entry(2, ro), entry(3, ro), entry(4, ro), entry(5, ro), entry(6, rw), entry(7, rw), entry(8, rw)],
        });
        let g1 = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor { label: Some("coarse"), entries: &[entry(0, uni), entry(1, ro), entry(2, rw), entry(3, rw)] });
        let g2 = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
            label: Some("cg"),
            entries: &[entry(0, ro), entry(1, rw), entry(2, rw), entry(3, rw), entry(4, rw), entry(5, rw)],
        });
        let pl = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor { label: Some("mg"), bind_group_layouts: &[Some(&g0), Some(&g1), Some(&g2)], immediate_size: 0 });
        let pipe = |name: &str| {
            device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor { label: Some(name), layout: Some(&pl), module: &module, entry_point: Some(name), compilation_options: Default::default(), cache: None })
        };
        let pipes = MgPipes {
            apply_t: pipe("apply_t"),
            apply_q: pipe("apply_q"),
            zero_z: pipe("zero_z"),
            smooth: pipe("jacobi_smooth"),
            resid: pipe("resid"),
            restrict_r: pipe("restrict_r"),
            prolong_add: pipe("prolong_add"),
            copy_u_p: pipe("copy_u_p"),
            init_r: pipe("init_r"),
            set_p_z: pipe("set_p_z"),
            reduce_rz: pipe("reduce_rz"),
            reduce_pq: pipe("reduce_pq"),
            final_init: pipe("final_init"),
            final_alpha: pipe("final_alpha"),
            update_ur: pipe("update_ur"),
            final_beta: pipe("final_beta"),
            update_p: pipe("update_p"),
        };
        Ok(Self { device, queue, adapter_name, pipes, g0, g1, g2 })
    }

    /// Solves K u = f for the grid (`modulus` per cell times the unit element matrix for `nu`),
    /// like `topo2d::solve` with the multigrid solver and `settings`.
    #[allow(clippy::too_many_arguments)]
    pub fn solve(&self, nx: usize, ny: usize, nu: f64, modulus: &[f64], fixed: &[bool], f: &[f64], start: Option<&[f64]>, settings: Settings, tol: f64, max: usize) -> Result<(GpuSolve, MgTimes), String> {
        let mut times = MgTimes::default();
        let t0 = std::time::Instant::now();
        let ke = crate::topo2d::element_stiffness(nu);
        let fine = Level { nx, ny, cells: Cells::Scaled { ke: Box::new(ke), modulus: modulus.to_vec() }, fixed: fixed.to_vec() };
        let n = fine.ndof();
        if f.len() != n || modulus.len() != nx * ny {
            return Err("GPU solve: input sizes do not match the grid".into());
        }
        let h = Hierarchy::new(fine, settings);
        times.hierarchy = t0.elapsed().as_secs_f64();
        let t1 = std::time::Instant::now();
        let d = &self.device;
        let st = wgpu::BufferUsages::STORAGE;
        let cp = wgpu::BufferUsages::COPY_SRC | wgpu::BufferUsages::COPY_DST;
        let buf = |label: &str, bytes: &[u8], usage: wgpu::BufferUsages| d.create_buffer_init(&wgpu::util::BufferInitDescriptor { label: Some(label), contents: bytes, usage });
        let dummy = buf("dummy", f64_bytes(&[0.0]), st);
        let flat_ke: Vec<f64> = ke.iter().flatten().copied().collect();
        let ke_b = buf("ke", f64_bytes(&flat_ke), st);
        let npart = (n as u32).div_ceil(WG);
        let depth = h.levels.len();
        // Per level: params, cells, fixed, inverse diagonal, r, z, t.
        struct Raw {
            params: wgpu::Buffer,
            cparams: wgpu::Buffer,
            cells: wgpu::Buffer,
            fixed: wgpu::Buffer,
            invd: wgpu::Buffer,
            r: wgpu::Buffer,
            z: wgpu::Buffer,
            t: wgpu::Buffer,
            scaled: bool,
            nx: usize,
            ny: usize,
            ndof: usize,
        }
        let mut raws = Vec::with_capacity(depth);
        for (l, level) in h.levels.iter().enumerate() {
            let ndof = level.ndof();
            let (scaled, cells) = match &level.cells {
                Cells::Scaled { modulus, .. } => (true, buf("modulus", f64_bytes(modulus), st)),
                Cells::Full(all) => {
                    let flat: Vec<f64> = all.iter().flat_map(|k| k.iter().flatten().copied()).collect();
                    (false, buf("mats", f64_bytes(&flat), st))
                }
            };
            let mask: Vec<u32> = level.fixed.iter().map(|b| *b as u32).collect();
            let invd = h.inv_diag(l).ok_or("missing diagonal")?;
            let zero = vec![0.0f64; ndof];
            raws.push(Raw {
                params: buf(
                    "lp",
                    bytemuck::bytes_of(&LevelParams { nx: level.nx as u32, ny: level.ny as u32, ndof: ndof as u32, scaled: scaled as u32, npart, omega: settings.omega as f32, pad0: 0, pad1: 0 }),
                    wgpu::BufferUsages::UNIFORM,
                ),
                cparams: buf("cp", bytemuck::bytes_of(&CoarseParams { nx: level.nx as u32, ny: level.ny as u32, ndof: ndof as u32, pad: 0 }), wgpu::BufferUsages::UNIFORM),
                cells,
                fixed: buf("fixed", bytemuck::cast_slice(&mask), st),
                invd: buf("invd", f64_bytes(invd), st),
                r: buf("r", f64_bytes(&zero), st | cp),
                z: buf("z", f64_bytes(&zero), st | cp),
                t: buf("t", f64_bytes(&zero), st),
                scaled,
                nx: level.nx,
                ny: level.ny,
                ndof,
            });
        }
        let zero_n = vec![0.0f64; n];
        let fb = buf("f", f64_bytes(f), st);
        let u = buf("u", f64_bytes(start.filter(|s| s.len() == n).unwrap_or(&zero_n)), st | cp);
        let p = buf("p", f64_bytes(&zero_n), st);
        let q = buf("q", f64_bytes(&zero_n), st);
        let part = buf("part", f64_bytes(&vec![0.0; 2 * npart as usize]), st);
        let sc = buf("sc", f64_bytes(&[0.0; 8]), st | cp);
        let bg = |layout: &wgpu::BindGroupLayout, bufs: &[&wgpu::Buffer]| {
            d.create_bind_group(&wgpu::BindGroupDescriptor {
                label: None,
                layout,
                entries: &bufs.iter().enumerate().map(|(i, b)| wgpu::BindGroupEntry { binding: i as u32, resource: b.as_entire_binding() }).collect::<Vec<_>>(),
            })
        };
        let levels: Vec<GpuLevel> = raws
            .iter()
            .map(|w| GpuLevel {
                nx: w.nx,
                ny: w.ny,
                ndof: w.ndof,
                r: w.r.clone(),
                z: w.z.clone(),
                g0: bg(&self.g0, &[&w.params, &ke_b, if w.scaled { &w.cells } else { &dummy }, if w.scaled { &dummy } else { &w.cells }, &w.fixed, &w.invd, &w.r, &w.z, &w.t]),
                g1: bg(&self.g1, &[&w.cparams, &w.fixed, &w.r, &w.z]),
            })
            .collect();
        let g2 = bg(&self.g2, &[&fb, &u, &p, &q, &part, &sc]);
        times.upload = t1.elapsed().as_secs_f64();
        let t2 = std::time::Instant::now();

        let last = depth - 1;
        let coarse = levels.get(last).ok_or("no levels")?;
        let nodes = |l: &GpuLevel| ((l.nx + 1) * (l.ny + 1)) as u32;
        let groups = |count: u32| count.div_ceil(WG);
        // All dispatches of one half-step go into a single compute pass (wgpu orders them and
        // inserts the barriers); a pass per dispatch cost more than the kernels on the fine grid.
        let run = |pass: &mut wgpu::ComputePass, pipe: &wgpu::ComputePipeline, l: usize, x: u32| {
            let lv = &levels[l];
            let cg1 = &levels[(l + 1).min(last)].g1;
            pass.set_pipeline(pipe);
            pass.set_bind_group(0, &lv.g0, &[]);
            pass.set_bind_group(1, cg1, &[]);
            pass.set_bind_group(2, &g2, &[]);
            pass.dispatch_workgroups(x, 1, 1);
        };
        let sweeps = settings.sweeps;
        let down = |pass: &mut wgpu::ComputePass| {
            for l in 0..last {
                let lv = &levels[l];
                run(pass, &self.pipes.zero_z, l, groups(lv.ndof as u32));
                for _ in 0..sweeps {
                    run(pass, &self.pipes.apply_t, l, groups(nodes(lv)));
                    run(pass, &self.pipes.smooth, l, groups(lv.ndof as u32));
                }
                run(pass, &self.pipes.apply_t, l, groups(nodes(lv)));
                run(pass, &self.pipes.resid, l, groups(lv.ndof as u32));
                run(pass, &self.pipes.restrict_r, l, groups(nodes(&levels[l + 1])));
            }
        };
        let up = |pass: &mut wgpu::ComputePass| {
            for l in (0..last).rev() {
                let lv = &levels[l];
                run(pass, &self.pipes.prolong_add, l, groups(nodes(lv)));
                for _ in 0..sweeps {
                    run(pass, &self.pipes.apply_t, l, groups(nodes(lv)));
                    run(pass, &self.pipes.smooth, l, groups(lv.ndof as u32));
                }
            }
        };
        let ts = d.features().contains(wgpu::Features::TIMESTAMP_QUERY);
        let qs = ts.then(|| d.create_query_set(&wgpu::QuerySetDescriptor { label: Some("ts"), ty: wgpu::QueryType::Timestamp, count: 4 }));
        let encode_at = |f: &dyn Fn(&mut wgpu::ComputePass), slot: Option<u32>| {
            let mut enc = d.create_command_encoder(&Default::default());
            {
                let tw = match (&qs, slot) {
                    (Some(q), Some(k)) => Some(wgpu::ComputePassTimestampWrites { query_set: q, beginning_of_pass_write_index: Some(2 * k), end_of_pass_write_index: Some(2 * k + 1) }),
                    _ => None,
                };
                let mut pass = enc.begin_compute_pass(&wgpu::ComputePassDescriptor { label: None, timestamp_writes: tw });
                f(&mut pass);
            }
            enc
        };
        let encode = |f: &dyn Fn(&mut wgpu::ComputePass)| encode_at(f, None);
        let rc_bytes = (coarse.ndof * 8) as u64;
        let st_r = d.create_buffer(&wgpu::BufferDescriptor { label: Some("rc out"), size: rc_bytes, usage: wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST, mapped_at_creation: false });
        let st_sc = d.create_buffer(&wgpu::BufferDescriptor { label: Some("sc out"), size: 64, usage: wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST, mapped_at_creation: false });
        // Submit `enc` after copying the coarsest r and the scalars out; wait; solve the coarsest
        // level on the CPU and write its z back. Returns the scalars.
        let mut coarse_cpu = 0.0;
        let mut sync = |mut enc: wgpu::CommandEncoder| -> Result<[f64; 8], String> {
            enc.copy_buffer_to_buffer(&coarse.r, 0, &st_r, 0, rc_bytes);
            enc.copy_buffer_to_buffer(&sc, 0, &st_sc, 0, 64);
            self.queue.submit([enc.finish()]);
            st_r.slice(..).map_async(wgpu::MapMode::Read, |_| {});
            st_sc.slice(..).map_async(wgpu::MapMode::Read, |_| {});
            d.poll(wgpu::PollType::wait_indefinitely()).map_err(|e| format!("GPU poll: {e}"))?;
            let tc = std::time::Instant::now();
            let rc: Vec<f64> = bytemuck::cast_slice(&st_r.slice(..).get_mapped_range().map_err(|e| format!("GPU map: {e:?}"))?).to_vec();
            let mut s = [0.0; 8];
            s.copy_from_slice(bytemuck::cast_slice(&st_sc.slice(..).get_mapped_range().map_err(|e| format!("GPU map: {e:?}"))?));
            st_r.unmap();
            st_sc.unmap();
            let zc = h.coarse_solve(&rc);
            self.queue.write_buffer(&coarse.z, 0, f64_bytes(&zc));
            coarse_cpu += tc.elapsed().as_secs_f64();
            Ok(s)
        };
        let fnorm = f.iter().map(|v| v * v).sum::<f64>().sqrt().max(f64::MIN_POSITIVE);
        let dofs = groups(n as u32);
        let fine_nodes = groups(nodes(&levels[0]));
        // Start: r = f - K u, z = M r, p = z.
        sync(encode(&|pass| {
            run(pass, &self.pipes.copy_u_p, 0, dofs);
            run(pass, &self.pipes.apply_q, 0, fine_nodes);
            run(pass, &self.pipes.init_r, 0, dofs);
            down(pass);
        }))?;
        let finish_up = |first: bool| {
            encode_at(&|pass| {
                up(pass);
                if first {
                    run(pass, &self.pipes.set_p_z, 0, dofs);
                    run(pass, &self.pipes.reduce_rz, 0, npart);
                    run(pass, &self.pipes.final_init, 0, 1);
                } else {
                    run(pass, &self.pipes.reduce_rz, 0, npart);
                    run(pass, &self.pipes.final_beta, 0, 1);
                    run(pass, &self.pipes.update_p, 0, dofs);
                }
            }, Some(0))
        };
        // Each step: [the previous step's up half + p update] then [q = K p, alpha, u and r
        // update, down half], one submit, then the coarsest solve on the CPU. The scalars read at
        // a sync describe r before this step's update.
        let mut first = true;
        let mut it = 0;
        let converged;
        loop {
            let before = finish_up(first);
            first = false;
            self.queue.submit([before.finish()]);
            let s_prev = sync(encode_at(&|pass| {
                run(pass, &self.pipes.apply_q, 0, fine_nodes);
                run(pass, &self.pipes.reduce_pq, 0, npart);
                run(pass, &self.pipes.final_alpha, 0, 1);
                run(pass, &self.pipes.update_ur, 0, dofs);
                down(pass);
            }, Some(1)))?;
            if s_prev[3].sqrt() <= tol * fnorm || s_prev[4] != 0.0 || it >= max {
                // This step's update has already run; it only lowers the residual further.
                converged = s_prev[3].sqrt() <= tol * fnorm;
                break;
            }
            it += 1;
        }
        times.iterate = t2.elapsed().as_secs_f64();
        if let Some(q) = &qs {
            let resolve = d.create_buffer(&wgpu::BufferDescriptor { label: None, size: 32, usage: wgpu::BufferUsages::QUERY_RESOLVE | wgpu::BufferUsages::COPY_SRC, mapped_at_creation: false });
            let read = d.create_buffer(&wgpu::BufferDescriptor { label: None, size: 32, usage: wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST, mapped_at_creation: false });
            let mut enc = d.create_command_encoder(&Default::default());
            enc.resolve_query_set(q, 0..4, &resolve, 0);
            enc.copy_buffer_to_buffer(&resolve, 0, &read, 0, 32);
            self.queue.submit([enc.finish()]);
            read.slice(..).map_async(wgpu::MapMode::Read, |_| {});
            d.poll(wgpu::PollType::wait_indefinitely()).map_err(|e| format!("GPU poll: {e}"))?;
            let t: Vec<u64> = bytemuck::cast_slice(&read.slice(..).get_mapped_range().map_err(|e| format!("GPU map: {e:?}"))?).to_vec();
            read.unmap();
            let period = self.queue.get_timestamp_period() as f64 * 1e-9;
            times.gpu_up = t[1].wrapping_sub(t[0]) as f64 * period;
            times.gpu_down = t[3].wrapping_sub(t[2]) as f64 * period;
        }
        times.coarse_cpu = coarse_cpu;
        // Read u back.
        let out = d.create_buffer(&wgpu::BufferDescriptor { label: Some("u out"), size: (n * 8) as u64, usage: wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST, mapped_at_creation: false });
        let mut enc = d.create_command_encoder(&Default::default());
        enc.copy_buffer_to_buffer(&u, 0, &out, 0, (n * 8) as u64);
        self.queue.submit([enc.finish()]);
        out.slice(..).map_async(wgpu::MapMode::Read, |_| {});
        d.poll(wgpu::PollType::wait_indefinitely()).map_err(|e| format!("GPU poll: {e}"))?;
        let u_host: Vec<f64> = bytemuck::cast_slice(&out.slice(..).get_mapped_range().map_err(|e| format!("GPU map: {e:?}"))?).to_vec();
        out.unmap();
        Ok((GpuSolve { u: u_host, iterations: it, converged, levels: depth }, times))
    }
}

/// One GPU for the whole process, opened on first use. Solves from several threads take turns.
static SHARED: std::sync::OnceLock<Option<std::sync::Mutex<GpuMg>>> = std::sync::OnceLock::new();

/// A multigrid solve on the shared GPU, or None (with one warning on stderr) when no GPU opens
/// or the solve fails, so the caller falls back to the CPU.
#[allow(clippy::too_many_arguments)]
pub fn solve_shared(nx: usize, ny: usize, nu: f64, modulus: &[f64], fixed: &[bool], f: &[f64], start: Option<&[f64]>, settings: Settings, tol: f64, max: usize) -> Option<GpuSolve> {
    let gpu = SHARED.get_or_init(|| match GpuMg::new() {
        Ok(g) => {
            eprintln!("GPU solver: {}", g.adapter_name);
            Some(std::sync::Mutex::new(g))
        }
        Err(e) => {
            eprintln!("GPU solver not available ({e}); using the CPU");
            None
        }
    });
    let guard = gpu.as_ref()?.lock().ok()?;
    match guard.solve(nx, ny, nu, modulus, fixed, f, start, settings, tol, max) {
        Ok((s, _)) => Some(s),
        Err(e) => {
            eprintln!("GPU solve failed ({e}); this solve runs on the CPU");
            None
        }
    }
}
