// Multigrid-preconditioned conjugate gradients on the GPU, f64. Mirrors mg2d: bilinear
// prolongation, restriction = its transpose, damped Jacobi smoothing, Galerkin coarse cells.
// The coarsest level is solved on the CPU between the down and up halves of each V-cycle.
// DOF layout as mg2d::cell_dofs: node (i, j) -> dofs 2 * (i * (ny + 1) + j) + {0, 1}.

struct LP { nx: u32, ny: u32, ndof: u32, scaled: u32, npart: u32, omega: f32, pad0: u32, pad1: u32 }
@group(0) @binding(0) var<uniform> L: LP;
@group(0) @binding(1) var<storage, read> ke: array<f64, 64>;
@group(0) @binding(2) var<storage, read> modulus: array<f64>;
@group(0) @binding(3) var<storage, read> mats: array<f64>;
@group(0) @binding(4) var<storage, read> fixedm: array<u32>;
@group(0) @binding(5) var<storage, read> invd: array<f64>;
@group(0) @binding(6) var<storage, read_write> r: array<f64>;
@group(0) @binding(7) var<storage, read_write> z: array<f64>;
@group(0) @binding(8) var<storage, read_write> t: array<f64>;

struct CP { nx: u32, ny: u32, ndof: u32, pad: u32 }
@group(1) @binding(0) var<uniform> C: CP;
@group(1) @binding(1) var<storage, read> cfixed: array<u32>;
@group(1) @binding(2) var<storage, read_write> rc: array<f64>;
@group(1) @binding(3) var<storage, read_write> zc: array<f64>;

@group(2) @binding(0) var<storage, read> f: array<f64>;
@group(2) @binding(1) var<storage, read_write> u: array<f64>;
@group(2) @binding(2) var<storage, read_write> p: array<f64>;
@group(2) @binding(3) var<storage, read_write> q: array<f64>;
@group(2) @binding(4) var<storage, read_write> part: array<f64>;
// sc: 0 rz, 1 alpha, 2 beta, 3 rr, 4 breakdown flag
@group(2) @binding(5) var<storage, read_write> sc: array<f64, 8>;

var<workgroup> sh: array<f64, 256>;
var<workgroup> sh2: array<f64, 256>;

fn dof(a: u32, b: u32) -> u32 { return 2u * (a * (L.ny + 1u) + b); }

fn d4(a: vec4<f64>, b: vec4<f64>) -> f64 { return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w; }

// Rows 2 ln and 2 ln + 1 of cell e's matrix times its 8 DOF values (lo = DOFs 0-3, hi = 4-7).
// Written with fixed vectors: a dynamically indexed local array spills to scratch memory on AMD.
fn rows(e: u32, ln: u32, lo: vec4<f64>, hi: vec4<f64>) -> vec2<f64> {
    let r0 = 2u * ln;
    if (L.scaled == 1u) {
        let m = modulus[e];
        let k = r0 * 8u;
        let a0 = vec4<f64>(ke[k], ke[k + 1u], ke[k + 2u], ke[k + 3u]);
        let b0 = vec4<f64>(ke[k + 4u], ke[k + 5u], ke[k + 6u], ke[k + 7u]);
        let a1 = vec4<f64>(ke[k + 8u], ke[k + 9u], ke[k + 10u], ke[k + 11u]);
        let b1 = vec4<f64>(ke[k + 12u], ke[k + 13u], ke[k + 14u], ke[k + 15u]);
        return m * vec2<f64>(d4(a0, lo) + d4(b0, hi), d4(a1, lo) + d4(b1, hi));
    }
    let k = e * 64u + r0 * 8u;
    let a0 = vec4<f64>(mats[k], mats[k + 1u], mats[k + 2u], mats[k + 3u]);
    let b0 = vec4<f64>(mats[k + 4u], mats[k + 5u], mats[k + 6u], mats[k + 7u]);
    let a1 = vec4<f64>(mats[k + 8u], mats[k + 9u], mats[k + 10u], mats[k + 11u]);
    let b1 = vec4<f64>(mats[k + 12u], mats[k + 13u], mats[k + 14u], mats[k + 15u]);
    return vec2<f64>(d4(a0, lo) + d4(b0, hi), d4(a1, lo) + d4(b1, hi));
}

// The cells around node n: (cell x, cell y, local node), invalid ones flagged by x < 0.
fn cell_of(n: u32, c: u32) -> vec3<i32> {
    let i = i32(n / (L.ny + 1u));
    let j = i32(n % (L.ny + 1u));
    var ex = i;
    var ey = j;
    if (c == 1u) { ex = i - 1; }
    if (c == 2u) { ex = i - 1; ey = j - 1; }
    if (c == 3u) { ey = j - 1; }
    if (ex < 0 || ey < 0 || ex >= i32(L.nx) || ey >= i32(L.ny)) { return vec3<i32>(-1, 0, 0); }
    return vec3<i32>(ex, ey, i32(c));
}

// t = K z on this level (zero on fixed DOFs).
@compute @workgroup_size(256)
fn apply_t(@builtin(global_invocation_id) gid: vec3<u32>) {
    let n = gid.x;
    if (n >= (L.nx + 1u) * (L.ny + 1u)) { return; }
    var acc = vec2<f64>(f64(0.0), f64(0.0));
    for (var c = 0u; c < 4u; c++) {
        let cl = cell_of(n, c);
        if (cl.x < 0) { continue; }
        let x = u32(cl.x);
        let y = u32(cl.y);
        let d0 = dof(x, y); let d1 = dof(x + 1u, y); let d2 = dof(x + 1u, y + 1u); let d3 = dof(x, y + 1u);
        let lo = vec4<f64>(z[d0], z[d0 + 1u], z[d1], z[d1 + 1u]);
        let hi = vec4<f64>(z[d2], z[d2 + 1u], z[d3], z[d3 + 1u]);
        acc += rows(x * L.ny + y, u32(cl.z), lo, hi);
    }
    let d = 2u * n;
    t[d] = select(acc.x, f64(0.0), fixedm[d] != 0u);
    t[d + 1u] = select(acc.y, f64(0.0), fixedm[d + 1u] != 0u);
}

// q = K p on the fine level (CG).
@compute @workgroup_size(256)
fn apply_q(@builtin(global_invocation_id) gid: vec3<u32>) {
    let n = gid.x;
    if (n >= (L.nx + 1u) * (L.ny + 1u)) { return; }
    var acc = vec2<f64>(f64(0.0), f64(0.0));
    for (var c = 0u; c < 4u; c++) {
        let cl = cell_of(n, c);
        if (cl.x < 0) { continue; }
        let x = u32(cl.x);
        let y = u32(cl.y);
        let d0 = dof(x, y); let d1 = dof(x + 1u, y); let d2 = dof(x + 1u, y + 1u); let d3 = dof(x, y + 1u);
        let lo = vec4<f64>(p[d0], p[d0 + 1u], p[d1], p[d1 + 1u]);
        let hi = vec4<f64>(p[d2], p[d2 + 1u], p[d3], p[d3 + 1u]);
        acc += rows(x * L.ny + y, u32(cl.z), lo, hi);
    }
    let d = 2u * n;
    q[d] = select(acc.x, f64(0.0), fixedm[d] != 0u);
    q[d + 1u] = select(acc.y, f64(0.0), fixedm[d + 1u] != 0u);
}

@compute @workgroup_size(256)
fn zero_z(@builtin(global_invocation_id) gid: vec3<u32>) {
    if (gid.x < L.ndof) { z[gid.x] = f64(0.0); }
}

// z += omega D^-1 (r - t)
@compute @workgroup_size(256)
fn jacobi_smooth(@builtin(global_invocation_id) gid: vec3<u32>) {
    let k = gid.x;
    if (k >= L.ndof) { return; }
    z[k] = select(z[k] + f64(L.omega) * invd[k] * (r[k] - t[k]), f64(0.0), fixedm[k] != 0u);
}

// t = r - t
@compute @workgroup_size(256)
fn resid(@builtin(global_invocation_id) gid: vec3<u32>) {
    let k = gid.x;
    if (k >= L.ndof) { return; }
    t[k] = select(r[k] - t[k], f64(0.0), fixedm[k] != 0u);
}

fn w_of(fine: i32, coarse: i32) -> f64 {
    // weight with which fine node index `fine` takes coarse index `coarse` (mg2d::parents)
    if (fine == 2 * coarse) { return f64(1.0); }
    if (fine == 2 * coarse + 1 || fine == 2 * coarse - 1) { return f64(0.5); }
    return f64(0.0);
}

// rc = R t: one thread per coarse node (C dims; L is the fine level).
@compute @workgroup_size(256)
fn restrict_r(@builtin(global_invocation_id) gid: vec3<u32>) {
    let n = gid.x;
    if (n >= (C.nx + 1u) * (C.ny + 1u)) { return; }
    let ci = i32(n / (C.ny + 1u));
    let cj = i32(n % (C.ny + 1u));
    var a = vec2<f64>(f64(0.0), f64(0.0));
    for (var di = -1; di <= 1; di++) {
        let i = 2 * ci + di;
        if (i < 0 || i > i32(L.nx)) { continue; }
        for (var dj = -1; dj <= 1; dj++) {
            let j = 2 * cj + dj;
            if (j < 0 || j > i32(L.ny)) { continue; }
            let w = w_of(i, ci) * w_of(j, cj);
            let d = dof(u32(i), u32(j));
            if (fixedm[d] == 0u) { a.x += w * t[d]; }
            if (fixedm[d + 1u] == 0u) { a.y += w * t[d + 1u]; }
        }
    }
    let dc = 2u * n;
    rc[dc] = select(a.x, f64(0.0), cfixed[dc] != 0u);
    rc[dc + 1u] = select(a.y, f64(0.0), cfixed[dc + 1u] != 0u);
}

// z += P zc: one thread per fine node.
@compute @workgroup_size(256)
fn prolong_add(@builtin(global_invocation_id) gid: vec3<u32>) {
    let n = gid.x;
    if (n >= (L.nx + 1u) * (L.ny + 1u)) { return; }
    let i = i32(n / (L.ny + 1u));
    let j = i32(n % (L.ny + 1u));
    var a = vec2<f64>(f64(0.0), f64(0.0));
    for (var ci = i / 2; ci <= (i + 1) / 2; ci++) {
        if (ci > i32(C.nx)) { continue; }
        for (var cj = j / 2; cj <= (j + 1) / 2; cj++) {
            if (cj > i32(C.ny)) { continue; }
            let w = w_of(i, ci) * w_of(j, cj);
            let dc = 2u * (u32(ci) * (C.ny + 1u) + u32(cj));
            a.x += w * zc[dc];
            a.y += w * zc[dc + 1u];
        }
    }
    let d = 2u * n;
    if (fixedm[d] == 0u) { z[d] += a.x; }
    if (fixedm[d + 1u] == 0u) { z[d + 1u] += a.y; }
}

// ---- CG on the fine level (group 0 = level 0: r, z; group 2: f, u, p, q) ----

@compute @workgroup_size(256)
fn copy_u_p(@builtin(global_invocation_id) gid: vec3<u32>) {
    let k = gid.x;
    if (k >= L.ndof) { return; }
    let v = select(u[k], f64(0.0), fixedm[k] != 0u);
    u[k] = v;
    p[k] = v;
}

@compute @workgroup_size(256)
fn init_r(@builtin(global_invocation_id) gid: vec3<u32>) {
    let k = gid.x;
    if (k >= L.ndof) { return; }
    r[k] = select(f[k] - q[k], f64(0.0), fixedm[k] != 0u);
}

@compute @workgroup_size(256)
fn set_p_z(@builtin(global_invocation_id) gid: vec3<u32>) {
    if (gid.x < L.ndof) { p[gid.x] = z[gid.x]; }
}

fn tree(lid: u32) {
    for (var s = 128u; s > 0u; s = s >> 1u) {
        workgroupBarrier();
        if (lid < s) { sh[lid] += sh[lid + s]; sh2[lid] += sh2[lid + s]; }
    }
    workgroupBarrier();
}

@compute @workgroup_size(256)
fn reduce_rz(@builtin(global_invocation_id) gid: vec3<u32>, @builtin(local_invocation_id) lid: vec3<u32>, @builtin(workgroup_id) wid: vec3<u32>) {
    let k = gid.x;
    var a = f64(0.0);
    var b = f64(0.0);
    if (k < L.ndof) { a = r[k] * z[k]; b = r[k] * r[k]; }
    sh[lid.x] = a;
    sh2[lid.x] = b;
    tree(lid.x);
    if (lid.x == 0u) { part[wid.x] = sh[0]; part[L.npart + wid.x] = sh2[0]; }
}

@compute @workgroup_size(256)
fn reduce_pq(@builtin(global_invocation_id) gid: vec3<u32>, @builtin(local_invocation_id) lid: vec3<u32>, @builtin(workgroup_id) wid: vec3<u32>) {
    let k = gid.x;
    var a = f64(0.0);
    if (k < L.ndof) { a = p[k] * q[k]; }
    sh[lid.x] = a;
    sh2[lid.x] = f64(0.0);
    tree(lid.x);
    if (lid.x == 0u) { part[wid.x] = sh[0]; }
}

fn final_sums(lid: u32) {
    var a = f64(0.0);
    var b = f64(0.0);
    for (var k = lid; k < L.npart; k += 256u) { a += part[k]; b += part[L.npart + k]; }
    sh[lid] = a;
    sh2[lid] = b;
    tree(lid);
}

@compute @workgroup_size(256)
fn final_init(@builtin(local_invocation_id) lid: vec3<u32>) {
    final_sums(lid.x);
    if (lid.x == 0u) { sc[0] = sh[0]; sc[3] = sh2[0]; sc[4] = f64(0.0); }
}

@compute @workgroup_size(256)
fn final_alpha(@builtin(local_invocation_id) lid: vec3<u32>) {
    final_sums(lid.x);
    if (lid.x == 0u) {
        let pq = sh[0];
        if (pq > f64(0.0)) { sc[1] = sc[0] / pq; } else { sc[1] = f64(0.0); sc[4] = f64(1.0); }
    }
}

// u += alpha p, r -= alpha q
@compute @workgroup_size(256)
fn update_ur(@builtin(global_invocation_id) gid: vec3<u32>) {
    let k = gid.x;
    if (k >= L.ndof) { return; }
    let al = sc[1];
    u[k] += al * p[k];
    r[k] = select(r[k] - al * q[k], f64(0.0), fixedm[k] != 0u);
}

@compute @workgroup_size(256)
fn final_beta(@builtin(local_invocation_id) lid: vec3<u32>) {
    final_sums(lid.x);
    if (lid.x == 0u) {
        let rzn = sh[0];
        let old = sc[0];
        sc[2] = select(rzn / old, f64(0.0), old == f64(0.0));
        sc[0] = rzn;
        sc[3] = sh2[0];
    }
}

@compute @workgroup_size(256)
fn update_p(@builtin(global_invocation_id) gid: vec3<u32>) {
    let k = gid.x;
    if (k >= L.ndof) { return; }
    p[k] = z[k] + sc[2] * p[k];
}
