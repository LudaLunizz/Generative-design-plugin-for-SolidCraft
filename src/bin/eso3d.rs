//! 3D "take out what is not needed": start from a solid pyramid (frustum) on the wall and repeatedly
//! remove the voxels that carry the least stress (hard-kill ESO), with the real 3D solver in stress3d.
//!
//! eso3d [key=value ...]   (all lengths mm)
//!   reach=150 base_y=80 base_z=100 tip=20   frustum: wall face x=0, tip face x=reach, tip x tip
//!   voxel=2.5  force=50 side=15            tip load: force N down (-Z) and side N along Y (two cases)
//!   back=6 pad=8                           kept solid: wall plate thickness (screws), load pad depth
//!   min_member=5                           no member thinner than this (opening after every round; 0 = off)
//!   holes=y:z,y:z                          countersunk screw holes (4.5 mm, 8.4 mm at 90 deg) through the wall plate, in the smooth STL
//!   from=DIR blur=2                        only re-smooth DIR/mask.bin (no ESO run); blur passes
//!   target=0.12 step=0.04 sf=4             stop at this share of the frustum, remove 4 % per round, keep SF
//!   material=petg-printed out=DIR
//!   check=FILE.stl [check_voxel=1.25]      only check an STL (same supports, loads, material) on a finer grid, heat maps in out/
//! Stops at the target share or just before the safety factor (away from the load/support zones) drops
//! below sf, and writes DIR/bracket.stl (voxel faces, wall at x=0, Z up), report.json and heat maps.

#![allow(clippy::indexing_slicing, clippy::cast_precision_loss)] // command-line tool: indices come from the grid itself

use generative_design_engine::library::Library;
use generative_design_engine::stress3d::{Case, Grid, Region, Setup, analyse, png, projection, read_stl, voxelize};
use std::collections::HashMap;
use std::time::Instant;

fn stl(grid: &Grid, solid: &[bool]) -> Vec<u8> {
    let at = |i: i64, j: i64, k: i64| {
        let ok = i >= 0 && j >= 0 && k >= 0 && (i as usize) < grid.n[0] && (j as usize) < grid.n[1] && (k as usize) < grid.n[2];
        ok && solid[grid.cell_index(i as usize, j as usize, k as usize)]
    };
    // (axis, side, other two axes in the order that winds outward)
    let faces = [(0, 1, 1, 2), (0, -1, 2, 1), (1, 1, 2, 0), (1, -1, 0, 2), (2, 1, 0, 1), (2, -1, 1, 0)];
    let mut tris: Vec<[[f32; 3]; 4]> = Vec::new(); // normal + 3 corners
    for c in (0..grid.cells()).filter(|c| solid[*c]) {
        let [i, j, k] = grid.cell_ijk(c);
        let ijk = [i as i64, j as i64, k as i64];
        for (ax, side, a, b) in faces {
            let mut n = ijk;
            n[ax] += side;
            if at(n[0], n[1], n[2]) {
                continue;
            }
            let base = |du: usize, dv: usize| {
                let mut p = [ijk[0] as f64, ijk[1] as f64, ijk[2] as f64];
                p[ax] += if side > 0 { 1.0 } else { 0.0 };
                p[a] += du as f64;
                p[b] += dv as f64;
                [0, 1, 2].map(|d| (grid.origin[d] + p[d] * grid.h) as f32)
            };
            let mut nor = [0.0f32; 3];
            nor[ax] = side as f32;
            let (p00, p10, p11, p01) = (base(0, 0), base(1, 0), base(1, 1), base(0, 1));
            tris.push([nor, p00, p10, p11]);
            tris.push([nor, p00, p11, p01]);
        }
    }
    let mut out = vec![0u8; 80];
    out.extend((tris.len() as u32).to_le_bytes());
    for t in &tris {
        for v in t {
            for x in v {
                out.extend(x.to_le_bytes());
            }
        }
        out.extend([0u8; 2]);
    }
    out
}


/// The same volume with a round form: the arms (the voxels that are not kept solid) are blurred and their
/// surface cut where the whole shape gives back exactly the voxel volume; the kept wall plate and tip pad stay exact boxes (marching tetrahedra on a grid twice as fine). The wall
/// face (x = 0) and the tip face stay flat: the design is continued one voxel past them, and that skin is
/// pressed flat. Returns triangles, the level used and the volume of the result in mm3.
/// The user's standard screw: countersunk, hole 4.5 mm, countersink 8.4 mm at 90 degrees.
const HOLE_R: f64 = 2.25;
const SINK_R: f64 = 4.2;

/// Signed distance (approximate inside the cone, exact on the cylinder) to the screw holes through the wall
/// plate (x = 0 to xb, countersink on the front face x = xb); negative inside a hole. Beyond x = xb + lim
/// there is no hole, so the cone never cuts an arm.
fn hole_dist(p: [f64; 3], holes: &[[f64; 2]], xb: f64, lim: f64) -> f64 {
    if p[0] > xb + lim {
        return f64::INFINITY;
    }
    holes.iter().map(|[y, z]| {
        let rho = ((p[1] - y).powi(2) + (p[2] - z).powi(2)).sqrt();
        (rho - HOLE_R).min((rho - SINK_R + xb - p[0]) / std::f64::consts::SQRT_2)
    }).fold(f64::INFINITY, f64::min)
}

/// Volume of one countersunk hole through a plate xb thick.
fn hole_volume(xb: f64) -> f64 {
    let hc = SINK_R - HOLE_R; // 90 degrees: as deep as the radius grows
    std::f64::consts::PI * (HOLE_R * HOLE_R * (xb - hc) + hc / 3.0 * (SINK_R * SINK_R + SINK_R * HOLE_R + HOLE_R * HOLE_R))
}

fn smooth(grid: &Grid, solid: &[bool], keep: &[bool], passes: usize, holes: &[[f64; 2]], xb: f64) -> (Vec<[[f64; 3]; 3]>, f64, f64) {
    const M: usize = 4; // empty margin in voxels
    let [nx, ny, nz] = grid.n;
    let d = [nx + 2 * M, ny + 2 * M, nz + 2 * M];
    let at = |i: usize, j: usize, k: usize| d[1] * d[2] * i + d[2] * j + k;
    // x: one voxel of the design continues past each end, then nothing
    let fill = |on: &dyn Fn(usize) -> bool| {
        let mut v = vec![0.0f64; d[0] * d[1] * d[2]];
        for i in 0..d[0] {
            let si = i as i64 - M as i64;
            if si < -1 || si > nx as i64 {
                continue;
            }
            let si = si.clamp(0, nx as i64 - 1) as usize;
            for j in M..M + ny {
                for k in M..M + nz {
                    v[at(i, j, k)] = f64::from(u8::from(on(grid.cell_index(si, j - M, k - M))));
                }
            }
        }
        v
    };
    // Where an arm meets the plate or the pad, the kept cells next to it join the blur, so the arm keeps its
    // full thickness into the plate (without this it fades out one voxel short and the part falls in two).
    let arm = |c: usize| solid[c] && !keep[c];
    let root: Vec<bool> = (0..grid.cells())
        .map(|c| {
            let [i, j, k] = grid.cell_ijk(c);
            keep[c]
                && (0..27).any(|q| {
                    let (a, b, e) = (i as i64 + q / 9 - 1, j as i64 + (q / 3) % 3 - 1, k as i64 + q % 3 - 1);
                    a >= 0 && b >= 0 && e >= 0 && (a as usize) < nx && (b as usize) < ny && (e as usize) < nz && arm(grid.cell_index(a as usize, b as usize, e as usize))
                })
        })
        .collect();
    let mut f = fill(&|c| arm(c) || root[c]);
    let kp = fill(&|c| keep[c]);
    for _ in 0..passes {
        for ax in 0..3 {
            let src = f.clone();
            let stride = [d[1] * d[2], d[2], 1][ax];
            for (idx, v) in f.iter_mut().enumerate() {
                let pos = (idx / stride) % d[ax];
                let lo = if pos > 0 { src[idx - stride] } else { src[idx] };
                let hi = if pos + 1 < d[ax] { src[idx + stride] } else { src[idx] };
                *v = (lo + src[idx] + hi) / 3.0;
            }
        }
    }
    let hf = grid.h / 2.0;
    let fd = [0, 1, 2].map(|a| 2 * (grid.n[a] + 2 * (M - 1)) + 1);
    let o = [0, 1, 2].map(|a| grid.origin[a] - (M - 1) as f64 * grid.h);
    let sample = |arr: &[f64], a: usize, b: usize, c: usize| {
        let w = [a, b, c].map(|v| v as f64 * 0.5 + 0.5); // coarse index space: fine node 0 sits half a cell inside the array
        let i0 = [0, 1, 2].map(|q| (w[q].floor() as usize).min(d[q] - 2));
        let t = [0, 1, 2].map(|q| w[q] - i0[q] as f64);
        let mut acc = 0.0;
        for corner in 0..8 {
            let (dx, dy, dz) = (corner & 1, (corner >> 1) & 1, (corner >> 2) & 1);
            let wt = [dx, dy, dz].iter().zip(&t).fold(1.0, |acc, (bit, tt)| acc * if *bit == 1 { *tt } else { 1.0 - tt });
            acc += wt * arr[at(i0[0] + dx, i0[1] + dy, i0[2] + dz)];
        }
        acc
    };
    let nidx = |a: usize, b: usize, c: usize| (a * fd[1] + b) * fd[2] + c;
    // fa: the blurred arms; fp: the kept plate and pad, 0.5 exactly on their faces
    let (mut fa, mut fp) = (vec![0.0f64; fd[0] * fd[1] * fd[2]], vec![0.0f64; fd[0] * fd[1] * fd[2]]);
    for a in 0..fd[0] {
        for b in 0..fd[1] {
            for c in 0..fd[2] {
                fa[nidx(a, b, c)] = sample(&f, a, b, c);
                fp[nidx(a, b, c)] = sample(&kp, a, b, c);
            }
        }
    }
    let pos = |a: usize, b: usize, c: usize| [o[0] + a as f64 * hf, o[1] + b as f64 * hf, o[2] + c as f64 * hf];
    let mut hd = vec![f64::INFINITY; fa.len()];
    if !holes.is_empty() {
        for a in 0..fd[0] {
            for b in 0..fd[1] {
                for c in 0..fd[2] {
                    hd[nidx(a, b, c)] = hole_dist(pos(a, b, c), holes, xb, hf);
                }
            }
        }
    }
    let reach = nx as f64 * grid.h;
    let want = solid.iter().filter(|v| **v).count() as f64 * grid.h.powi(3) - holes.len() as f64 * hole_volume(xb);
    let inside_nodes = |thr: f64| {
        let mut n = 0usize;
        for a in 0..fd[0] {
            let x = o[0] + a as f64 * hf;
            if x < 0.0 || x > reach {
                continue;
            }
            // the layers on the wall and tip faces count half: a node there owns half a cell
            let r = a * fd[1] * fd[2]..(a + 1) * fd[1] * fd[2];
            let layer = fa[r.clone()].iter().zip(&fp[r.clone()]).zip(&hd[r]).filter(|((av, pv), d)| (**av >= thr || **pv >= 0.5) && **d >= 0.0).count();
            n += if x < 0.5 * hf || x > reach - 0.5 * hf { layer } else { 2 * layer };
        }
        0.5 * n as f64 * hf.powi(3)
    };
    // Counting nodes only estimates the volume of the cut surface, so the target is corrected by the real volume a few times.
    let build = |want: f64| {
    let (mut lo, mut hi) = (0.0, 1.0);
    for _ in 0..40 {
        let mid = 0.5 * (lo + hi);
        if inside_nodes(mid) > want { lo = mid } else { hi = mid }
    }
    let thr = 0.5 * (lo + hi);
    // The plate field is bent so that its 0.5 face lands exactly on the level.
    // The holes: linear in their distance with the level on the hole wall, so the cut lands on the true circle.
    let field: Vec<f64> = fa.iter().zip(&fp).zip(&hd).map(|((av, pv), d)| av.max(if *pv >= 0.5 { thr + (pv - 0.5) * 2.0 * (1.0 - thr) } else { pv * 2.0 * thr }).min(thr + d / hf)).collect();
    // Six tetrahedra per cell around the diagonal 0-7.
    let tets = [[0, 1, 3, 7], [0, 3, 2, 7], [0, 2, 6, 7], [0, 6, 4, 7], [0, 4, 5, 7], [0, 5, 1, 7]];
    let mut tris = Vec::new();
    for a in 0..fd[0] - 1 {
        for b in 0..fd[1] - 1 {
            for c in 0..fd[2] - 1 {
                let corner: [(usize, usize, usize); 8] = std::array::from_fn(|q| (a + (q & 1), b + ((q >> 1) & 1), c + ((q >> 2) & 1)));
                for tet in tets {
                    let v: [(usize, usize, usize); 4] = tet.map(|q| corner[q]);
                    let inn: Vec<bool> = v.iter().map(|p| field[nidx(p.0, p.1, p.2)] >= thr).collect();
                    let k = inn.iter().filter(|x| **x).count();
                    if k == 0 || k == 4 {
                        continue;
                    }
                    // The crossing on the edge between two nodes; the same two nodes always give the same point.
                    let edge = |p: (usize, usize, usize), q: (usize, usize, usize)| {
                        let (p, q) = if nidx(p.0, p.1, p.2) < nidx(q.0, q.1, q.2) { (p, q) } else { (q, p) };
                        let (fp, fq) = (field[nidx(p.0, p.1, p.2)], field[nidx(q.0, q.1, q.2)]);
                        let t = (fp - thr) / (fp - fq);
                        let (pp, pq) = (pos(p.0, p.1, p.2), pos(q.0, q.1, q.2));
                        [0, 1, 2].map(|x| pp[x] + t * (pq[x] - pp[x]))
                    };
                    let ins: Vec<usize> = (0..4).filter(|q| inn[*q]).collect();
                    let out: Vec<usize> = (0..4).filter(|q| !inn[*q]).collect();
                    let centroid = |ix: &[usize]| {
                        let mut c3 = [0.0; 3];
                        for q in ix {
                            let p = pos(v[*q].0, v[*q].1, v[*q].2);
                            for x in 0..3 {
                                c3[x] += p[x] / ix.len() as f64;
                            }
                        }
                        c3
                    };
                    let dir = { let (ci, co) = (centroid(&ins), centroid(&out)); [0, 1, 2].map(|x| co[x] - ci[x]) };
                    let polys: Vec<[[f64; 3]; 3]> = match k {
                        1 | 3 => {
                            let (one, rest) = if k == 1 { (ins[0], &out) } else { (out[0], &ins) };
                            vec![[edge(v[one], v[rest[0]]), edge(v[one], v[rest[1]]), edge(v[one], v[rest[2]])]]
                        }
                        _ => {
                            let (e0, e1, e2, e3) = (edge(v[ins[0]], v[out[0]]), edge(v[ins[0]], v[out[1]]), edge(v[ins[1]], v[out[1]]), edge(v[ins[1]], v[out[0]]));
                            vec![[e0, e1, e2], [e0, e2, e3]]
                        }
                    };
                    for mut t in polys {
                        let (u, w) = ([0, 1, 2].map(|x| t[1][x] - t[0][x]), [0, 1, 2].map(|x| t[2][x] - t[0][x]));
                        let n = [u[1] * w[2] - u[2] * w[1], u[2] * w[0] - u[0] * w[2], u[0] * w[1] - u[1] * w[0]];
                        if n[0] * dir[0] + n[1] * dir[1] + n[2] * dir[2] < 0.0 {
                            t.swap(1, 2);
                        }
                        for p in &mut t {
                            p[0] = p[0].clamp(0.0, reach);
                        }
                        tris.push(t);
                    }
                }
            }
        }
    }
    // Linear crossings fall a little inside a round hole (about 0.1 mm here): put those corners on the hole wall.
    for p in tris.iter_mut().flatten().filter(|p| p[0] <= xb + 1e-9) {
        for [y, z] in holes {
            let (dy, dz) = (p[1] - y, p[2] - z);
            let (rho, want) = (dy.hypot(dz), HOLE_R.max(SINK_R - (xb - p[0])));
            if rho > 1e-9 && rho < want {
                (p[1], p[2]) = (y + dy * want / rho, z + dz * want / rho);
            }
        }
    }
    // Triangles pressed flat to nothing (the skin past the wall and tip faces) are dropped.
    tris.retain(|t| {
        let (u, w) = ([0, 1, 2].map(|x| t[1][x] - t[0][x]), [0, 1, 2].map(|x| t[2][x] - t[0][x]));
        let n = [u[1] * w[2] - u[2] * w[1], u[2] * w[0] - u[0] * w[2], u[0] * w[1] - u[1] * w[0]];
        n.iter().map(|v| v * v).sum::<f64>() > 1e-12
    });
    let vol = tris.iter().map(|t| (t[0][0] * (t[1][1] * t[2][2] - t[1][2] * t[2][1]) - t[0][1] * (t[1][0] * t[2][2] - t[1][2] * t[2][0]) + t[0][2] * (t[1][0] * t[2][1] - t[1][1] * t[2][0])) / 6.0).sum::<f64>();
    (tris, thr, vol)
    };
    let (mut target, mut res) = (want, build(want));
    for _ in 0..3 {
        target *= want / res.2;
        res = build(target);
    }
    res
}

fn tris_stl(tris: &[[[f64; 3]; 3]]) -> Vec<u8> {
    let mut out = vec![0u8; 80];
    out.extend((tris.len() as u32).to_le_bytes());
    for t in tris {
        let (u, w) = ([0, 1, 2].map(|x| t[1][x] - t[0][x]), [0, 1, 2].map(|x| t[2][x] - t[0][x]));
        let n = [u[1] * w[2] - u[2] * w[1], u[2] * w[0] - u[0] * w[2], u[0] * w[1] - u[1] * w[0]];
        let len = n.iter().map(|v| v * v).sum::<f64>().sqrt().max(f64::MIN_POSITIVE);
        for v in [n.map(|x| x / len), t[0], t[1], t[2]] {
            for x in v {
                out.extend((x as f32).to_le_bytes());
            }
        }
        out.extend([0u8; 2]);
    }
    out
}

/// Opening with a k x k x k block: what no such block fits into (members thinner than k voxels) goes.
fn open(grid: &Grid, solid: &[bool], k: usize) -> Vec<bool> {
    let [nx, ny, nz] = grid.n;
    let block = |c: usize, back: bool, test: &dyn Fn(usize) -> bool, all: bool| {
        let [i, j, l] = grid.cell_ijk(c);
        let mut hit = false; // all: every cell must pass; any: one must
        for a in 0..k {
            for b in 0..k {
                for d in 0..k {
                    let off = |v: usize, e: usize, n: usize| if back { v.checked_sub(e).filter(|x| *x < n) } else { Some(v + e).filter(|x| *x < n) };
                    let ok = match (off(i, a, nx), off(j, b, ny), off(l, d, nz)) {
                        (Some(x), Some(y), Some(z)) => test(grid.cell_index(x, y, z)),
                        _ => false,
                    };
                    if all && !ok {
                        return false;
                    }
                    if !all && ok {
                        hit = true;
                    }
                }
            }
        }
        all || hit
    };
    let eroded: Vec<bool> = (0..grid.cells()).map(|c| block(c, false, &|n| solid[n], true)).collect();
    (0..grid.cells()).map(|c| block(c, true, &|n| eroded[n], false)).collect()
}

/// "y:z,y:z" (mm on the wall face) -> hole centres; empty -> none.
fn parse_holes(s: &str) -> Result<Vec<[f64; 2]>, String> {
    s.split(',').filter(|t| !t.is_empty()).map(|t| {
        let (y, z) = t.split_once(':').ok_or(format!("holes: {t}: use y:z"))?;
        let n = |v: &str| v.parse::<f64>().ok().filter(|x| x.is_finite()).ok_or(format!("holes: {t}: numbers"));
        Ok([n(y)?, n(z)?])
    }).collect()
}

fn main() {
    if let Err(e) = run() {
        eprintln!("eso3d: {e}");
        std::process::exit(2);
    }
}

fn run() -> Result<(), String> {
    let mut o: HashMap<String, String> = HashMap::new();
    for a in std::env::args().skip(1) {
        let (k, v) = a.split_once('=').ok_or(format!("{a}: use key=value"))?;
        o.insert(k.into(), v.into());
    }
    let num = |k: &str, d: f64| -> Result<f64, String> { o.get(k).map_or(Ok(d), |v| v.parse().map_err(|_| format!("{k}: a number"))) };
    let (reach, by, bz, tip, h) = (num("reach", 150.0)?, num("base_y", 80.0)?, num("base_z", 100.0)?, num("tip", 20.0)?, num("voxel", 2.5)?);
    let (force, side, back, pad) = (num("force", 50.0)?, num("side", 15.0)?, num("back", 6.0)?, num("pad", 8.0)?);
    let (target, step, sf) = (num("target", 0.12)?, num("step", 0.04)?, num("sf", 4.0)?);
    let kmin = (num("min_member", 5.0)? / h).round().max(1.0) as usize;
    if !(0.5..=10.0).contains(&h) || !(0.0..1.0).contains(&target) || !(0.005..=0.3).contains(&step) || tip > by.min(bz) {
        return Err("check voxel (0.5-10), target (0-1), step (0.005-0.3), tip <= base".into());
    }
    let out = o.get("out").cloned().unwrap_or_else(|| "eso3d-out".into());
    let m = Library::builtin()?.find(o.get("material").map_or("petg-printed", String::as_str)).ok_or("material not in the library")?.clone();
    let s = m.strengths();
    let (e_gpa, nu) = (s.young_modulus_gpa.ok_or("no stiffness")?, s.poisson_ratio.ok_or("no Poisson ratio")?);
    let strength = s.yield_mpa.or(s.tensile_mpa).ok_or("no strength")?;

    let n = [(reach / h).round() as usize, (by / h).round() as usize, (bz / h).round() as usize];
    let grid = Grid { origin: [0.0, -0.5 * n[1] as f64 * h, -0.5 * n[2] as f64 * h], h, n };
    let mut solid = vec![false; grid.cells()];
    let mut keep = vec![false; grid.cells()];
    for c in 0..grid.cells() {
        let p = grid.cell_centre(c);
        let t = p[0] / (n[0] as f64 * h); // 0 at the wall, 1 at the tip
        let (hy, hz) = ((by + (tip - by) * t) / 2.0, (bz + (tip - bz) * t) / 2.0);
        solid[c] = p[1].abs() <= hy && p[2].abs() <= hz;
        keep[c] = solid[c] && (p[0] < back || p[0] > n[0] as f64 * h - pad);
    }
    let start = solid.iter().filter(|v| **v).count();
    let cell_cm3 = h.powi(3) / 1000.0;
    eprintln!("frustum {}x{}x{} voxels of {h} mm: {start} cells = {:.0} cm3, kept solid {}", n[0], n[1], n[2], start as f64 * cell_cm3, keep.iter().filter(|v| **v).count());

    let fixed = Region::Side { axis: 0, max: false, depth: 0.5 * h };
    let load = Region::Side { axis: 0, max: true, depth: 0.5 * h };
    let cases = [
        Case { name: "down".into(), force: [0.0, 0.0, -force], moment: [0.0; 3] },
        Case { name: "side".into(), force: [0.0, side, 0.0], moment: [0.0; 3] },
    ];
    if let Some(file) = o.get("check") {
        let cv = num("check_voxel", 1.25)?;
        let (pos, tri) = read_stl(&std::fs::read(file).map_err(|e| format!("{file}: {e}"))?)?;
        let (g, sol) = voxelize(&pos, &tri, cv)?;
        let cells = sol.iter().filter(|v| **v).count();
        eprintln!("{file}: {}x{}x{} voxels of {cv} mm, {cells} cells = {:.1} cm3", g.n[0], g.n[1], g.n[2], cells as f64 * cv.powi(3) / 1000.0);
        {
            // how much of the material is joined to the wall (cells touching by a corner count as joined)
            let mut seen = vec![false; g.cells()];
            let mut stack: Vec<usize> = (0..g.cells()).filter(|c| sol[*c] && g.cell_ijk(*c)[0] == 0).collect();
            for c in &stack {
                seen[*c] = true;
            }
            let (mut reached, mut tip) = (0usize, 0usize);
            while let Some(c) = stack.pop() {
                reached += 1;
                let [i, j, k] = g.cell_ijk(c);
                if i + 1 == g.n[0] {
                    tip += 1;
                }
                for (di, dj, dk) in (0..27).map(|q| (q / 9 - 1i64, (q / 3) % 3 - 1i64, q % 3 - 1i64)) {
                    let (a, b, d) = (i as i64 + di, j as i64 + dj, k as i64 + dk);
                    if a < 0 || b < 0 || d < 0 || a as usize >= g.n[0] || b as usize >= g.n[1] || d as usize >= g.n[2] {
                        continue;
                    }
                    let n = g.cell_index(a as usize, b as usize, d as usize);
                    if sol[n] && !seen[n] {
                        seen[n] = true;
                        stack.push(n);
                    }
                }
            }
            eprintln!("joined to the wall: {reached} of {cells} cells, {tip} cells of the tip face");
        }
        let setup = Setup { grid: &g, solid: &sol, e_mpa: e_gpa * 1000.0, nu, fixed: &fixed, load: &Region::Side { axis: 0, max: true, depth: 0.5 * cv }, cg_tol: 1e-8, cg_max: 50_000, exclude_layers: 2, cell_scale: None, multigrid: true };
        let fin = analyse(&setup, &cases, |_| {})?;
        std::fs::create_dir_all(&out).map_err(|e| format!("{out}: {e}"))?;
        for c in &fin.cases {
            for (ax, name) in [(0, "x"), (1, "y"), (2, "z")] {
                let (w, hh, rgb) = projection(&g, &c.vm, ax, c.max_vm_away, 4);
                std::fs::write(format!("{out}/{}-along-{name}.png", c.name), png(w, hh, &rgb)).map_err(|e| e.to_string())?;
            }
            println!("{}: peak {:.2} MPa, away from load/support {:.2} MPa (SF {:.1} on {strength} MPa), p99 {:.2}, tip moves {:.3} mm, CG {} converged {}", c.name, c.max_vm, c.max_vm_away, strength / c.max_vm_away.max(1e-12), c.p99_vm, c.max_disp, c.cg_iterations, c.converged);
        }
        for w in &fin.warnings {
            println!("warning: {w}");
        }
        return Ok(());
    }
    let t0 = Instant::now();
    let mut prev = solid.clone(); // the design of the round before, for the fall-back
    let mut round = 0;
    let mut last_count = usize::MAX;
    let from = o.get("from").cloned();
    let reason = if let Some(dir) = &from {
        let bytes = std::fs::read(format!("{dir}/mask.bin")).map_err(|e| format!("{dir}/mask.bin: {e}"))?;
        if bytes.len() != grid.cells() {
            return Err("mask.bin does not fit this grid (same reach/base/tip/voxel needed)".into());
        }
        solid = bytes.iter().map(|b| *b != 0).collect();
        "smoothing only".to_string()
    } else { loop {
        let setup = Setup { grid: &grid, solid: &solid, e_mpa: e_gpa * 1000.0, nu, fixed: &fixed, load: &load, cg_tol: 1e-6, cg_max: 20_000, exclude_layers: 2, cell_scale: None, multigrid: true };
        let r = analyse(&setup, &cases, |_| {})?;
        let away = r.cases.iter().map(|c| c.max_vm_away).fold(0.0, f64::max);
        let count = r.cases[0].vm.iter().filter(|v| !v.is_nan()).count();
        let found = if away > 0.0 { strength / away } else { f64::INFINITY };
        let disp = r.cases.iter().map(|c| c.max_disp).fold(0.0, f64::max);
        eprintln!("round {round:>3}: {:.1} cm3 ({:.0} %), away peak {away:.2} MPa (SF {found:.1}), tip {disp:.3} mm, CG {}, {:.0} s", count as f64 * cell_cm3, 100.0 * count as f64 / start as f64, r.cases[0].cg_iterations, t0.elapsed().as_secs_f64());
        if found < sf && round > 0 {
            solid = prev;
            break format!("stopped: this round had SF {found:.1} < {sf}, kept the one before");
        }
        // What the solver kept (loose voxels are NaN) is the new design; this round's result is its analysis.
        for (c, v) in solid.iter_mut().enumerate() {
            *v = !r.cases[0].vm[c].is_nan();
        }
        if round > 0 && count >= last_count {
            break "no progress: every member is at the minimum size".to_string();
        }
        last_count = count;
        if count as f64 / start as f64 <= target {
            break format!("reached the target share {target}");
        }
        let score: Vec<f64> = (0..grid.cells()).map(|c| r.cases.iter().map(|k| k.vm[c]).filter(|v| !v.is_nan()).fold(0.0, f64::max)).collect();
        let mut cand: Vec<usize> = (0..grid.cells()).filter(|c| solid[*c] && !keep[*c]).collect();
        cand.sort_by(|a, b| score[*a].total_cmp(&score[*b]));
        let remove = ((count as f64 * step).ceil() as usize).min(((count as f64 - target * start as f64).max(1.0)) as usize).min(cand.len());
        if remove == 0 {
            break "nothing left to remove".into();
        }
        prev = solid.clone();
        for c in &cand[..remove] {
            solid[*c] = false;
        }
        if kmin > 1 {
            solid = open(&grid, &solid, kmin);
            for (v, kp) in solid.iter_mut().zip(&keep) {
                *v |= *kp;
            }
        }
        round += 1;
        if round > 200 {
            break "round limit".into();
        }
    } };
    let count = solid.iter().filter(|v| **v).count();
    std::fs::create_dir_all(&out).map_err(|e| format!("{out}: {e}"))?;
    std::fs::write(format!("{out}/bracket.stl"), stl(&grid, &solid)).map_err(|e| e.to_string())?;
    std::fs::write(format!("{out}/mask.bin"), solid.iter().map(|v| u8::from(*v)).collect::<Vec<u8>>()).map_err(|e| e.to_string())?;
    // the kept plate is whole voxels: every column whose centre lies within `back`
    let xb = (back / h - 0.5).ceil().max(0.0) * h;
    let holes = parse_holes(o.get("holes").map_or("", String::as_str))?;
    for [y, z] in &holes {
        if xb < SINK_R - HOLE_R + 1.0 || y.abs() + SINK_R > 0.5 * by || z.abs() + SINK_R > 0.5 * bz {
            return Err(format!("hole {y}:{z}: needs a plate of at least {:.2} mm (it is {xb}) and the countersink inside the plate", SINK_R - HOLE_R + 1.0));
        }
        // the screwdriver goes in along x: no arm within the countersink radius plus 2 mm
        let blocked = (0..grid.cells()).filter(|c| solid[*c] && !keep[*c]).any(|c| {
            let p = grid.cell_centre(c);
            (p[1] - y).hypot(p[2] - z) < SINK_R + 2.0 + 0.5 * h
        });
        if blocked {
            eprintln!("warning: hole {y}:{z}: an arm is in the screwdriver's way");
        }
    }
    let (tris, thr, vol) = smooth(&grid, &solid, &keep, num("blur", 2.0)? as usize, &holes, xb);
    std::fs::write(format!("{out}/bracket-smooth.stl"), tris_stl(&tris)).map_err(|e| e.to_string())?;
    println!("smooth: {} triangles, level {thr:.3}, {:.1} cm3 against {:.1} cm3 as voxels", tris.len(), vol / 1000.0, count as f64 * cell_cm3);
    if from.is_some() {
        return Ok(());
    }
    // final analysis of exactly what was written
    let setup = Setup { grid: &grid, solid: &solid, e_mpa: e_gpa * 1000.0, nu, fixed: &fixed, load: &load, cg_tol: 1e-8, cg_max: 50_000, exclude_layers: 2, cell_scale: None, multigrid: true };
    let fin = analyse(&setup, &cases, |_| {})?;
    let mut rep = String::new();
    for c in &fin.cases {
        for (ax, name) in [(0, "x"), (1, "y"), (2, "z")] {
            let (w, hh, rgb) = projection(&grid, &c.vm, ax, c.max_vm_away, 4);
            std::fs::write(format!("{out}/{}-along-{name}.png", c.name), png(w, hh, &rgb)).map_err(|e| e.to_string())?;
        }
        rep += &format!("{}: away peak {:.2} MPa (SF {:.1} on {strength} MPa), p99 {:.2}, tip moves {:.3} mm\n", c.name, c.max_vm_away, strength / c.max_vm_away.max(1e-12), c.p99_vm, c.max_disp);
    }
    rep += &format!("{} cells = {:.1} cm3 = {:.1} g ({:.0} % of the frustum), {reason}, {:.0} s\n", count, count as f64 * cell_cm3, count as f64 * cell_cm3 * m.strengths().density_kg_m3.unwrap_or(1270.0) / 1000.0, 100.0 * count as f64 / start as f64, t0.elapsed().as_secs_f64());
    std::fs::write(format!("{out}/report.txt"), &rep).map_err(|e| e.to_string())?;
    print!("{rep}");
    Ok(())
}

#[cfg(test)]
mod tests {
    #![allow(clippy::unwrap_used, clippy::indexing_slicing)]
    use super::*;

    #[test]
    fn one_cell_has_twelve_triangles() {
        let g = Grid { origin: [0.0; 3], h: 1.0, n: [1, 1, 1] };
        let b = stl(&g, &[true]);
        assert_eq!(u32::from_le_bytes(b[80..84].try_into().unwrap()), 12);
        assert_eq!(b.len(), 84 + 12 * 50);
    }

    #[test]
    fn opening_removes_thin_members_and_keeps_thick_ones() {
        let g = Grid { origin: [0.0; 3], h: 1.0, n: [8, 4, 4] };
        let mut s = vec![false; 128];
        for i in 0..8 {
            s[g.cell_index(i, 0, 0)] = true; // a 1-voxel line
        }
        for i in 0..2 {
            for j in 2..4 {
                for k in 2..4 {
                    s[g.cell_index(i, j, k)] = true; // a 2 x 2 x 2 block
                }
            }
        }
        let o = open(&g, &s, 2);
        assert!(!o[g.cell_index(3, 0, 0)] && o[g.cell_index(1, 3, 3)]);
        assert_eq!(o.iter().filter(|v| **v).count(), 8);
    }

    #[test]
    fn smooth_keeps_the_volume_and_the_flat_faces() {
        let g = Grid { origin: [0.0, -3.0, -3.0], h: 1.0, n: [6, 6, 6] };
        let (tris, _, vol) = smooth(&g, &vec![true; 216], &vec![false; 216], 2, &[], 0.0);
        assert!((vol - 216.0).abs() < 0.05 * 216.0, "volume {vol}");
        assert!(tris.iter().flatten().all(|p| (0.0..=6.0).contains(&p[0])));
        // a kept block stays an exact block
        let (tris, _, vol) = smooth(&g, &vec![true; 216], &vec![true; 216], 2, &[], 0.0);
        assert!((vol - 216.0).abs() < 0.03 * 216.0, "kept volume {vol}"); // the edges lose a little
        assert!(tris.iter().flatten().all(|p| p[1].abs() <= 3.001 && p[2].abs() <= 3.001));
    }

    #[test]
    fn screw_hole_is_round_and_takes_its_own_volume() {
        // a 5 mm wall plate (kept) on a 10 mm grid: its front face is free, as on the bracket
        let g = Grid { origin: [0.0, -6.0, -6.0], h: 1.0, n: [10, 12, 12] };
        let plate: Vec<bool> = (0..g.cells()).map(|c| g.cell_ijk(c)[0] < 5).collect();
        let (_, _, plain) = smooth(&g, &plate, &plate, 2, &[], 5.0);
        let (tris, _, holed) = smooth(&g, &plate, &plate, 2, &[[0.0, 0.0]], 5.0);
        let hv = hole_volume(5.0);
        assert!(((plain - holed) - hv).abs() < 0.08 * hv, "removed {} against {hv}", plain - holed);
        // on the straight part of the hole every corner lies on the 4.5 mm circle
        let wall: Vec<f64> = tris.iter().flatten().filter(|p| p[0] > 0.5 && p[0] < 2.5).map(|p| p[1].hypot(p[2])).filter(|r| *r < 3.5).collect();
        assert!(!wall.is_empty() && wall.iter().all(|r| (r - HOLE_R).abs() < 0.05), "{wall:?}");
        // and no corner anywhere inside the hole, the wall face included
        let inside: Vec<&[f64; 3]> = tris.iter().flatten().filter(|p| p[1].hypot(p[2]) < HOLE_R - 1e-9).collect();
        assert!(inside.is_empty(), "{inside:?}");
        assert_eq!(parse_holes("-6.5:41,-6.5:-41").unwrap(), vec![[-6.5, 41.0], [-6.5, -41.0]]);
        assert!(parse_holes("1:x").is_err() && parse_holes("3").is_err() && parse_holes("").unwrap().is_empty());
    }
}
