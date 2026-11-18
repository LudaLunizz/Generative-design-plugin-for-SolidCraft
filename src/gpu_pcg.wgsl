// Matrix-free Jacobi-preconditioned conjugate gradients for the 2D plane-stress grid, f64.
// DOF layout as mg2d::cell_dofs: node (i, j) -> dofs 2 * (i * (ny + 1) + j) + {0, 1}.
struct Params { nx: u32, ny: u32, ndof: u32, npart: u32 }
@group(0) @binding(0) var<uniform> P: Params;
@group(0) @binding(1) var<storage, read> ke: array<f64, 64>;
@group(0) @binding(2) var<storage, read> modulus: array<f64>;
@group(0) @binding(3) var<storage, read> fixedm: array<u32>;
@group(0) @binding(4) var<storage, read> invd: array<f64>;
@group(0) @binding(5) var<storage, read> f: array<f64>;
@group(0) @binding(6) var<storage, read_write> u: array<f64>;
@group(0) @binding(7) var<storage, read_write> r: array<f64>;
@group(0) @binding(8) var<storage, read_write> z: array<f64>;
@group(0) @binding(9) var<storage, read_write> p: array<f64>;
@group(0) @binding(10) var<storage, read_write> q: array<f64>;
@group(0) @binding(11) var<storage, read_write> part: array<f64>;
// sc: 0 rz, 1 alpha, 2 beta, 3 rr, 4 breakdown flag
@group(0) @binding(12) var<storage, read_write> sc: array<f64, 8>;

var<workgroup> sh: array<f64, 256>;
var<workgroup> sh2: array<f64, 256>;

fn ndof_of(a: u32, b: u32) -> u32 { return 2u * (a * (P.ny + 1u) + b); }

// q = K p, zero on fixed DOFs. One thread per node, gathering from its (up to) 4 cells.
@compute @workgroup_size(256)
fn apply_pq(@builtin(global_invocation_id) gid: vec3<u32>) {
    let nn = (P.nx + 1u) * (P.ny + 1u);
    let n = gid.x;
    if (n >= nn) { return; }
    let i = i32(n / (P.ny + 1u));
    let j = i32(n % (P.ny + 1u));
    var a0 = f64(0.0);
    var a1 = f64(0.0);
    // (cell x, cell y, local node of this node in that cell)
    let cx = array<i32, 4>(i, i - 1, i - 1, i);
    let cy = array<i32, 4>(j, j, j - 1, j - 1);
    for (var c = 0u; c < 4u; c++) {
        let ex = cx[c];
        let ey = cy[c];
        if (ex < 0 || ey < 0 || ex >= i32(P.nx) || ey >= i32(P.ny)) { continue; }
        let x = u32(ex);
        let y = u32(ey);
        let m = modulus[x * P.ny + y];
        let d0 = ndof_of(x, y);
        let d1 = ndof_of(x + 1u, y);
        let d2 = ndof_of(x + 1u, y + 1u);
        let d3 = ndof_of(x, y + 1u);
        let ue = array<f64, 8>(p[d0], p[d0 + 1u], p[d1], p[d1 + 1u], p[d2], p[d2 + 1u], p[d3], p[d3 + 1u]);
        let row0 = 2u * c * 8u;
        var s0 = f64(0.0);
        var s1 = f64(0.0);
        for (var k = 0u; k < 8u; k++) {
            s0 += ke[row0 + k] * ue[k];
            s1 += ke[row0 + 8u + k] * ue[k];
        }
        a0 += m * s0;
        a1 += m * s1;
    }
    let d = 2u * n;
    q[d] = select(a0, f64(0.0), fixedm[d] != 0u);
    q[d + 1u] = select(a1, f64(0.0), fixedm[d + 1u] != 0u);
}

@compute @workgroup_size(256)
fn copy_u_p(@builtin(global_invocation_id) gid: vec3<u32>) {
    let k = gid.x;
    if (k >= P.ndof) { return; }
    let v = select(u[k], f64(0.0), fixedm[k] != 0u);
    u[k] = v;
    p[k] = v;
}

// After q = K u: r = f - q, z = D^-1 r, p = z.
@compute @workgroup_size(256)
fn init_r(@builtin(global_invocation_id) gid: vec3<u32>) {
    let k = gid.x;
    if (k >= P.ndof) { return; }
    let rv = select(f[k] - q[k], f64(0.0), fixedm[k] != 0u);
    r[k] = rv;
    let zv = invd[k] * rv;
    z[k] = zv;
    p[k] = zv;
}

fn tree(lid: u32) {
    for (var s = 128u; s > 0u; s = s >> 1u) {
        workgroupBarrier();
        if (lid < s) {
            sh[lid] += sh[lid + s];
            sh2[lid] += sh2[lid + s];
        }
    }
    workgroupBarrier();
}

// Partial sums per workgroup: part[w] = sum r.z, part[npart + w] = sum r.r
@compute @workgroup_size(256)
fn reduce_rz(@builtin(global_invocation_id) gid: vec3<u32>, @builtin(local_invocation_id) lid: vec3<u32>, @builtin(workgroup_id) wid: vec3<u32>) {
    let k = gid.x;
    var a = f64(0.0);
    var b = f64(0.0);
    if (k < P.ndof) { a = r[k] * z[k]; b = r[k] * r[k]; }
    sh[lid.x] = a;
    sh2[lid.x] = b;
    tree(lid.x);
    if (lid.x == 0u) { part[wid.x] = sh[0]; part[P.npart + wid.x] = sh2[0]; }
}

// part[w] = sum p.q
@compute @workgroup_size(256)
fn reduce_pq(@builtin(global_invocation_id) gid: vec3<u32>, @builtin(local_invocation_id) lid: vec3<u32>, @builtin(workgroup_id) wid: vec3<u32>) {
    let k = gid.x;
    var a = f64(0.0);
    if (k < P.ndof) { a = p[k] * q[k]; }
    sh[lid.x] = a;
    sh2[lid.x] = f64(0.0);
    tree(lid.x);
    if (lid.x == 0u) { part[wid.x] = sh[0]; }
}

fn final_sums(lid: u32) {
    var a = f64(0.0);
    var b = f64(0.0);
    for (var k = lid; k < P.npart; k += 256u) { a += part[k]; b += part[P.npart + k]; }
    sh[lid] = a;
    sh2[lid] = b;
    tree(lid);
}

// One workgroup: rz and rr at the start.
@compute @workgroup_size(256)
fn final_init(@builtin(local_invocation_id) lid: vec3<u32>) {
    final_sums(lid.x);
    if (lid.x == 0u) { sc[0] = sh[0]; sc[3] = sh2[0]; sc[4] = f64(0.0); }
}

// One workgroup: alpha = rz / (p.q).
@compute @workgroup_size(256)
fn final_alpha(@builtin(local_invocation_id) lid: vec3<u32>) {
    final_sums(lid.x);
    if (lid.x == 0u) {
        let pq = sh[0];
        if (pq > f64(0.0)) { sc[1] = sc[0] / pq; } else { sc[1] = f64(0.0); sc[4] = f64(1.0); }
    }
}

// u += alpha p, r -= alpha q, z = D^-1 r
@compute @workgroup_size(256)
fn update_ur(@builtin(global_invocation_id) gid: vec3<u32>) {
    let k = gid.x;
    if (k >= P.ndof) { return; }
    let al = sc[1];
    u[k] += al * p[k];
    let rv = select(r[k] - al * q[k], f64(0.0), fixedm[k] != 0u);
    r[k] = rv;
    z[k] = invd[k] * rv;
}

// One workgroup: beta = rz_new / rz, rz = rz_new, rr.
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

// p = z + beta p
@compute @workgroup_size(256)
fn update_p(@builtin(global_invocation_id) gid: vec3<u32>) {
    let k = gid.x;
    if (k >= P.ndof) { return; }
    p[k] = z[k] + sc[2] * p[k];
}
