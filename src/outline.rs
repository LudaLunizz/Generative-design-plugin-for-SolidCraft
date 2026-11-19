//! From a density grid to closed outlines that SolidCraft can draw as sketch lines.
//!
//! The optimiser gives one density per square cell. Steps: threshold to solid/void, keep only what
//! is connected to the fixed edge, close diagonal-only contacts (they would draw as outlines that
//! touch at a point), trace the boundary between solid and void cells as closed loops of grid
//! corners, cut the staircase corners (edge midpoints) and simplify to straight segments.
//!
//! Loops run counter-clockwise around solid material, so a hole is a clockwise loop. No indexing
//! or unwrap here: every lookup is checked.

use std::collections::{BTreeMap, VecDeque};

/// A grid of solid/void cells, column-major like the optimiser (`i * ny + j`, row `j` upward).
#[derive(Clone, Debug, PartialEq)]
pub struct Mask {
    pub nx: usize,
    pub ny: usize,
    solid: Vec<bool>,
}

/// A grid corner in units of half a cell (so edge midpoints stay whole numbers).
type Half = (i64, i64);

impl Mask {
    /// Cells with density above `threshold` are solid. `None` when the sizes do not fit.
    pub fn from_density(density: &[f64], nx: usize, ny: usize, threshold: f64) -> Option<Self> {
        if nx == 0 || ny == 0 || density.len() != nx.checked_mul(ny)? {
            return None;
        }
        Some(Self { nx, ny, solid: density.iter().map(|&d| d.is_finite() && d > threshold).collect() })
    }

    pub fn get(&self, i: i64, j: i64) -> bool {
        if i < 0 || j < 0 {
            return false;
        }
        let (Ok(i), Ok(j)) = (usize::try_from(i), usize::try_from(j)) else { return false };
        if i >= self.nx || j >= self.ny {
            return false;
        }
        self.solid.get(i * self.ny + j).copied().unwrap_or(false)
    }

    fn set(&mut self, i: usize, j: usize, value: bool) {
        if i < self.nx && j < self.ny
            && let Some(cell) = self.solid.get_mut(i * self.ny + j)
        {
            *cell = value;
        }
    }

    pub fn count(&self) -> usize {
        self.solid.iter().filter(|&&s| s).count()
    }

    /// Keeps only the cells joined (edge to edge) to a solid cell of the left column, the fixed
    /// edge. Returns how many cells were dropped.
    pub fn keep_connected_to_left(&mut self) -> usize {
        let before = self.count();
        let mut keep = vec![false; self.solid.len()];
        let mut queue = VecDeque::new();
        for j in 0..self.ny {
            if self.get(0, j as i64)
                && let Some(k) = keep.get_mut(j)
            {
                *k = true;
                queue.push_back((0usize, j));
            }
        }
        while let Some((i, j)) = queue.pop_front() {
            let (ii, jj) = (i as i64, j as i64);
            for (di, dj) in [(1i64, 0i64), (-1, 0), (0, 1), (0, -1)] {
                let (a, b) = (ii + di, jj + dj);
                if !self.get(a, b) {
                    continue;
                }
                let (Ok(a), Ok(b)) = (usize::try_from(a), usize::try_from(b)) else { continue };
                if let Some(k) = keep.get_mut(a * self.ny + b)
                    && !*k
                {
                    *k = true;
                    queue.push_back((a, b));
                }
            }
        }
        self.solid = keep;
        before - self.count()
    }

    /// Where two solid cells (or two void cells) meet only at a corner, fills one void cell so
    /// the contact becomes an edge. Returns the number of cells filled; `Err` if it does not settle.
    pub fn close_diagonal_contacts(&mut self) -> Result<usize, String> {
        let mut filled = 0;
        for _ in 0..64 {
            let mut changed = false;
            for i in 0..self.nx.saturating_sub(1) {
                for j in 0..self.ny.saturating_sub(1) {
                    let (ii, jj) = (i as i64, j as i64);
                    let (a, b, c, d) = (self.get(ii, jj), self.get(ii + 1, jj), self.get(ii, jj + 1), self.get(ii + 1, jj + 1));
                    if a == d && b == c && a != b {
                        let (fi, fj) = if !a { (i, j) } else { (i + 1, j) };
                        self.set(fi, fj, true);
                        filled += 1;
                        changed = true;
                    }
                }
            }
            if !changed {
                return Ok(filled);
            }
        }
        Err("the outline has diagonal contacts that do not settle".into())
    }

    /// Traces every boundary as a closed loop of corners (grid units, whole numbers) with solid on
    /// the left. Call [`Mask::close_diagonal_contacts`] first.
    pub fn trace(&self) -> Result<Vec<Vec<(i64, i64)>>, String> {
        let mut next: BTreeMap<(i64, i64), Vec<(i64, i64)>> = BTreeMap::new();
        let mut edges = 0usize;
        for i in 0..self.nx as i64 {
            for j in 0..self.ny as i64 {
                if !self.get(i, j) {
                    continue;
                }
                let mut add = |from: (i64, i64), to: (i64, i64)| {
                    next.entry(from).or_default().push(to);
                    edges += 1;
                };
                if !self.get(i, j - 1) {
                    add((i, j), (i + 1, j));
                }
                if !self.get(i + 1, j) {
                    add((i + 1, j), (i + 1, j + 1));
                }
                if !self.get(i, j + 1) {
                    add((i + 1, j + 1), (i, j + 1));
                }
                if !self.get(i - 1, j) {
                    add((i, j + 1), (i, j));
                }
            }
        }
        let mut loops = Vec::new();
        let mut used = 0usize;
        while let Some(&start) = next.keys().next() {
            let mut points = vec![start];
            let mut cur = start;
            loop {
                let to = next.get_mut(&cur).and_then(Vec::pop).ok_or("the outline is not closed")?;
                if next.get(&cur).is_some_and(Vec::is_empty) {
                    next.remove(&cur);
                }
                used += 1;
                if to == start {
                    break;
                }
                points.push(to);
                cur = to;
                if points.len() > edges {
                    return Err("the outline does not close".into());
                }
            }
            loops.push(points);
        }
        if used != edges {
            return Err("the outline lost edges".into());
        }
        Ok(loops)
    }
}

/// The unit-step loop with its staircase cut: every edge contributes its midpoint, and a corner
/// is kept only where both straight runs meeting there are at least [`KEEP_CORNER_RUN`] cells long
/// (a real corner of the part). Points on a straight run are removed. Half-cell whole numbers.
const KEEP_CORNER_RUN: usize = 3;

fn midpoints(corners: &[(i64, i64)]) -> Vec<Half> {
    let n = corners.len();
    let dir = |k: usize| -> (i64, i64) {
        match (corners.get(k % n), corners.get((k + 1) % n)) {
            (Some(a), Some(b)) => (b.0 - a.0, b.1 - a.1),
            _ => (0, 0),
        }
    };
    // Length (in unit edges) of the straight run each edge belongs to.
    let Some(start) = (0..n).find(|&k| dir(k) != dir(k + n - 1)) else { return Vec::new() };
    let mut run_len = vec![0usize; n];
    let mut k = 0;
    while k < n {
        let first = start + k;
        let mut len = 1;
        while k + len < n && dir(first + len) == dir(first) {
            len += 1;
        }
        for m in 0..len {
            if let Some(slot) = run_len.get_mut((first + m) % n) {
                *slot = len;
            }
        }
        k += len;
    }
    let mut pts: Vec<Half> = Vec::with_capacity(2 * n);
    for k in 0..n {
        let (Some(a), Some(b)) = (corners.get(k), corners.get((k + 1) % n)) else { continue };
        pts.push((a.0 + b.0, a.1 + b.1));
        // The corner after edge k.
        let next = (k + 1) % n;
        let turn = dir(k) != dir(next);
        let long = |e: usize| run_len.get(e).copied().unwrap_or(0) >= KEEP_CORNER_RUN;
        if turn && long(k) && long(next) {
            pts.push((2 * b.0, 2 * b.1));
        }
    }
    remove_collinear(&pts)
}

fn remove_collinear(points: &[Half]) -> Vec<Half> {
    let n = points.len();
    if n < 3 {
        return points.to_vec();
    }
    let mut out = Vec::with_capacity(n);
    for k in 0..n {
        let (Some(p), Some(q), Some(r)) = (points.get((k + n - 1) % n), points.get(k), points.get((k + 1) % n)) else { continue };
        let cross = (q.0 - p.0) * (r.1 - q.1) - (q.1 - p.1) * (r.0 - q.0);
        let dot = (q.0 - p.0) * (r.0 - q.0) + (q.1 - p.1) * (r.1 - q.1);
        if cross != 0 || dot < 0 {
            out.push(*q);
        }
    }
    out
}

/// Douglas-Peucker on a closed loop. `tol` is in the same units as the points.
fn simplify_closed(points: &[(f64, f64)], tol: f64) -> Vec<(f64, f64)> {
    let n = points.len();
    if n <= 3 {
        return points.to_vec();
    }
    let Some(first) = points.first() else { return Vec::new() };
    let far = points
        .iter()
        .enumerate()
        .max_by(|(_, a), (_, b)| dist2(a, first).total_cmp(&dist2(b, first)))
        .map_or(0, |(k, _)| k);
    // A linear copy that starts at point 0 and ends back on it.
    let mut line: Vec<(f64, f64)> = points.to_vec();
    line.push(*first);
    let mut keep = vec![false; line.len()];
    for k in [0, far, n] {
        if let Some(slot) = keep.get_mut(k) {
            *slot = true;
        }
    }
    let mut stack = vec![(0usize, far), (far, n)];
    while let Some((lo, hi)) = stack.pop() {
        let (Some(a), Some(b)) = (line.get(lo), line.get(hi)) else { continue };
        let mut worst = (0.0f64, lo);
        for k in lo + 1..hi {
            if let Some(p) = line.get(k) {
                let d = point_segment(p, a, b);
                if d > worst.0 {
                    worst = (d, k);
                }
            }
        }
        if worst.0 > tol {
            if let Some(slot) = keep.get_mut(worst.1) {
                *slot = true;
            }
            stack.push((lo, worst.1));
            stack.push((worst.1, hi));
        }
    }
    line.iter().take(n).zip(&keep).filter(|(_, k)| **k).map(|(p, _)| *p).collect()
}

fn dist2(a: &(f64, f64), b: &(f64, f64)) -> f64 {
    (a.0 - b.0).powi(2) + (a.1 - b.1).powi(2)
}

fn point_segment(p: &(f64, f64), a: &(f64, f64), b: &(f64, f64)) -> f64 {
    let (dx, dy) = (b.0 - a.0, b.1 - a.1);
    let len2 = dx * dx + dy * dy;
    if len2 == 0.0 {
        return dist2(p, a).sqrt();
    }
    let t = (((p.0 - a.0) * dx + (p.1 - a.1) * dy) / len2).clamp(0.0, 1.0);
    dist2(p, &(a.0 + t * dx, a.1 + t * dy)).sqrt()
}

/// Signed area of a loop (positive counter-clockwise).
pub fn signed_area(points: &[(f64, f64)]) -> f64 {
    let n = points.len();
    let mut s = 0.0;
    for k in 0..n {
        if let (Some(a), Some(b)) = (points.get(k), points.get((k + 1) % n)) {
            s += a.0 * b.1 - b.0 * a.1;
        }
    }
    s / 2.0
}

/// The outlines of a mask as simplified loops in cell units (x right, y up, origin at the corner
/// of cell (0, 0)). `tol` is the largest deviation from the staircase, in cells.
pub fn outlines(mask: &Mask, tol: f64) -> Result<Vec<Vec<(f64, f64)>>, String> {
    let mut out = Vec::new();
    for corners in mask.trace()? {
        let mids = midpoints(&corners);
        let pts: Vec<(f64, f64)> = mids.iter().map(|&(x, y)| (x as f64 / 2.0, y as f64 / 2.0)).collect();
        let simple = simplify_closed(&pts, tol);
        if simple.len() >= 3 && signed_area(&simple).abs() > 1e-9 {
            out.push(simple);
        }
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    #![allow(clippy::unwrap_used, clippy::expect_used, clippy::panic, clippy::indexing_slicing)]
    use super::*;

    fn mask(rows: &[&str]) -> Mask {
        // rows are written top row first; the mask is column-major with row 0 at the bottom.
        let ny = rows.len();
        let nx = rows[0].len();
        let mut d = vec![0.0; nx * ny];
        for (r, line) in rows.iter().enumerate() {
            let j = ny - 1 - r;
            for (i, c) in line.chars().enumerate() {
                if c == '#' {
                    d[i * ny + j] = 1.0;
                }
            }
        }
        Mask::from_density(&d, nx, ny, 0.5).unwrap()
    }

    fn total_area(m: &Mask, tol: f64) -> f64 {
        outlines(m, tol).unwrap().iter().map(|l| signed_area(l)).sum()
    }

    #[test]
    fn bad_sizes_are_refused() {
        assert!(Mask::from_density(&[1.0; 5], 2, 3, 0.5).is_none());
        assert!(Mask::from_density(&[], 0, 0, 0.5).is_none());
    }

    #[test]
    fn nan_density_is_void() {
        let m = Mask::from_density(&[f64::NAN, 1.0], 1, 2, 0.5).unwrap();
        assert_eq!(m.count(), 1);
    }

    #[test]
    fn a_rectangle_is_one_ccw_loop_of_four_points() {
        let m = mask(&["######", "######", "######", "######"]);
        let l = outlines(&m, 0.25).unwrap();
        assert_eq!(l.len(), 1);
        assert_eq!(l[0].len(), 4, "{:?}", l[0]);
        assert!((signed_area(&l[0]) - 24.0).abs() < 1e-9, "real corners stay exact");
    }

    #[test]
    fn a_thin_end_gets_its_corners_cut_a_little() {
        let m = mask(&["####", "####"]);
        let l = outlines(&m, 0.25).unwrap();
        assert_eq!(l.len(), 1);
        assert!(signed_area(&l[0]) > 0.0);
        assert!((signed_area(&l[0]) - 8.0).abs() < 1.0, "{:?}", l[0]);
    }

    #[test]
    fn a_hole_is_a_clockwise_loop() {
        let m = mask(&["######", "#....#", "#....#", "######"]);
        let l = outlines(&m, 0.25).unwrap();
        assert_eq!(l.len(), 2);
        let areas: Vec<f64> = l.iter().map(|x| signed_area(x)).collect();
        assert!(areas.iter().any(|a| *a > 0.0) && areas.iter().any(|a| *a < 0.0), "{areas:?}");
        let net: f64 = areas.iter().sum();
        assert!((net - 16.0).abs() < 1.5, "net {net}");
    }

    #[test]
    fn only_what_touches_the_fixed_edge_is_kept() {
        let mut m = mask(&["#.#", "#.#", "#.#"]);
        // The right column is separate from the left one.
        assert_eq!(m.keep_connected_to_left(), 3);
        assert_eq!(m.count(), 3);
        let mut empty = mask(&[".#", ".#"]);
        empty.keep_connected_to_left();
        assert_eq!(empty.count(), 0);
        assert!(empty.trace().unwrap().is_empty());
    }

    #[test]
    fn diagonal_contacts_are_closed_and_traced() {
        let mut m = mask(&["#.", ".#"]);
        assert_eq!(m.close_diagonal_contacts().unwrap(), 1);
        let loops = m.trace().unwrap();
        assert_eq!(loops.len(), 1, "{loops:?}");
        let mut checker = mask(&["#.#.", ".#.#", "#.#.", ".#.#"]);
        assert!(checker.close_diagonal_contacts().is_ok());
        assert!(checker.trace().is_ok());
    }

    #[test]
    fn a_slanted_bar_is_few_points_and_keeps_its_area() {
        // A 45 degree staircase band: many unit steps, a few straight segments after simplifying.
        let n = 24;
        let mut rows = Vec::new();
        for r in 0..n {
            let mut line = String::new();
            for c in 0..n {
                let d = c as i64 - (n as i64 - 1 - r as i64);
                line.push(if (0..4).contains(&d) { '#' } else { '.' });
            }
            rows.push(line);
        }
        let refs: Vec<&str> = rows.iter().map(String::as_str).collect();
        let mut m = mask(&refs);
        m.close_diagonal_contacts().unwrap();
        let cells = m.count() as f64;
        let l = outlines(&m, 0.4).unwrap();
        assert_eq!(l.len(), 1);
        assert!(l[0].len() < 14, "{} points", l[0].len());
        let a = signed_area(&l[0]);
        assert!((a - cells).abs() / cells < 0.06, "area {a} vs {cells} cells");
    }

    #[test]
    fn simplified_loops_never_cross_themselves_on_random_blobs() {
        // A deterministic pseudo-random field; every loop must stay simple (no two segments cross).
        let (nx, ny) = (30usize, 20usize);
        let mut s = 12345u64;
        let mut d = vec![0.0; nx * ny];
        for v in &mut d {
            s = s.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            *v = if (s >> 33) % 100 < 60 { 1.0 } else { 0.0 };
        }
        let mut m = Mask::from_density(&d, nx, ny, 0.5).unwrap();
        m.close_diagonal_contacts().unwrap();
        for tol in [0.0, 0.25, 0.5] {
            let loops = outlines(&m, tol).unwrap();
            for l in &loops {
                assert!(is_simple(l), "tol {tol}: loop with {} points crosses itself", l.len());
            }
        }
    }

    fn is_simple(l: &[(f64, f64)]) -> bool {
        let n = l.len();
        for a in 0..n {
            for b in a + 1..n {
                if b == a + 1 || (a == 0 && b == n - 1) {
                    continue;
                }
                if crosses(l[a], l[(a + 1) % n], l[b], l[(b + 1) % n]) {
                    return false;
                }
            }
        }
        true
    }

    fn crosses(p: (f64, f64), q: (f64, f64), r: (f64, f64), s: (f64, f64)) -> bool {
        let o = |a: (f64, f64), b: (f64, f64), c: (f64, f64)| (b.0 - a.0) * (c.1 - a.1) - (b.1 - a.1) * (c.0 - a.0);
        let (d1, d2, d3, d4) = (o(p, q, r), o(p, q, s), o(r, s, p), o(r, s, q));
        d1 * d2 < -1e-12 && d3 * d4 < -1e-12
    }

    #[test]
    fn total_area_is_close_to_the_cell_count() {
        let m = mask(&["######", "##..##", "##..##", "######"]);
        let a = total_area(&m, 0.25);
        assert!((a - m.count() as f64).abs() < 2.0, "{a}");
    }
}
