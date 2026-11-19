//! From a density grid to closed, smooth outlines that SolidCraft can draw as sketch lines and splines.
//!
//! The optimiser gives one density per square cell. Steps: threshold to solid/void, keep only what
//! is connected to the fixed edge, close diagonal-only contacts (they would draw as outlines that
//! touch at a point), trace the boundary between solid and void cells as closed loops of grid
//! corners, sample them at the edge midpoints, smooth them (Taubin; the plate's own edges stay
//! pinned) and split them into straight runs along the plate's edges and curves inside it.
//!
//! Loops run counter-clockwise around solid material, so a hole is a clockwise loop. No indexing
//! or unwrap here: every lookup is checked.

use std::collections::{BTreeMap, HashMap, VecDeque};

/// A grid of solid/void cells, column-major like the optimiser (`i * ny + j`, row `j` upward).
#[derive(Clone, Debug, PartialEq)]
pub struct Mask {
    pub nx: usize,
    pub ny: usize,
    solid: Vec<bool>,
    /// The densities the mask came from (non-finite as 0), to place the outline between cells.
    density: Vec<f64>,
}


impl Mask {
    /// Cells with density above `threshold` are solid. `None` when the sizes do not fit.
    pub fn from_density(density: &[f64], nx: usize, ny: usize, threshold: f64) -> Option<Self> {
        if nx == 0 || ny == 0 || density.len() != nx.checked_mul(ny)? {
            return None;
        }
        Some(Self {
            nx,
            ny,
            solid: density.iter().map(|&d| d.is_finite() && d > threshold).collect(),
            density: density.iter().map(|&d| if d.is_finite() { d - threshold + 0.5 } else { 0.0 }).collect(),
        })
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

    /// Density of a cell shifted so the threshold sits at 0.5; 0 outside the grid.
    fn density_at(&self, i: i64, j: i64) -> f64 {
        let (Ok(iu), Ok(ju)) = (usize::try_from(i), usize::try_from(j)) else { return 0.0 };
        if iu >= self.nx || ju >= self.ny {
            return 0.0;
        }
        self.density.get(iu * self.ny + ju).copied().unwrap_or(0.0)
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

/// A point in cell units (x right, y up, origin at the corner of cell (0, 0)).
pub type Point = (f64, f64);

/// How one closed outline is drawn.
#[derive(Clone, Debug, PartialEq)]
pub enum Shape {
    /// Straight lines all round (an outline that lies entirely on the plate's edges).
    Polygon(Vec<Point>),
    /// One smooth closed curve through these points (an outline that never touches the plate's edges).
    Spline(Vec<Point>),
    /// Straight runs along the plate's edges and smooth curves between them, joined end to end.
    Mixed(Vec<Piece>),
}

#[derive(Clone, Debug, PartialEq)]
pub enum Piece {
    /// An open polyline.
    Line(Vec<Point>),
    /// An open smooth curve through these points.
    Curve(Vec<Point>),
}

#[derive(Clone, Debug, PartialEq)]
pub struct Outline {
    pub shape: Shape,
    /// Signed area of the smoothed loop (positive counter-clockwise, so holes are negative).
    pub area: f64,
    /// The smoothed loop itself, densely sampled (what the curves follow).
    pub dense: Vec<Point>,
}

impl Outline {
    /// Points handed to the sketch.
    pub fn drawn_points(&self) -> usize {
        match &self.shape {
            Shape::Polygon(p) | Shape::Spline(p) => p.len(),
            Shape::Mixed(pieces) => pieces.iter().map(|p| match p { Piece::Line(v) | Piece::Curve(v) => v.len() }).sum(),
        }
    }
}

/// Taubin smoothing weights (Taubin 1995, "A signal processing approach to fair surface design"):
/// a shrinking step `LAMBDA` then an inflating step `MU`, so the outline loses its stair steps
/// without shrinking. The same family of methods as Blender's Laplacian Smooth (Desbrun et al. 1999),
/// written here from the papers, not from Blender's (GPL) code.
const LAMBDA: f64 = 0.5;
const MU: f64 = -0.53;
/// Smoothing passes tried in turn; a pass count whose loops cross themselves or each other falls back
/// to the next, gentler one (0 = the stair steps' edge midpoints, which never cross).
const PASSES: [usize; 4] = [12, 6, 2, 0];
/// Roughly one spline point per this many dense points (about 1.5 cells apart).
const CURVE_STEP: usize = 3;
/// SolidCraft's limit for one fit-point spline.
const MAX_CURVE_POINTS: usize = 500;

/// The outlines of a mask, smoothed: every loop is traced on the cell edges, sampled at the edge
/// midpoints, and smoothed. Points on the plate's own edges (the fixed edge, the load and support
/// edges) are pinned, so those stay straight and in place; everything inside becomes curves.
pub fn smooth_outlines(mask: &Mask) -> Result<Vec<Outline>, String> {
    let loops: Vec<Vec<(Point, bool)>> = mask.trace()?.iter().map(|c| dense(c, mask)).collect();
    for passes in PASSES {
        let smoothed: Vec<Vec<(Point, bool)>> = loops.iter().map(|l| taubin(l, passes)).collect();
        let dense_loops: Vec<Vec<Point>> = smoothed.iter().map(|l| l.iter().map(|&(p, _)| p).collect()).collect();
        if passes > 0 && any_crossing(&dense_loops) {
            continue;
        }
        return Ok(smoothed
            .iter()
            .zip(dense_loops)
            .filter(|(l, _)| l.len() >= 3)
            .map(|(l, d)| Outline { shape: shape(l), area: signed_area(&d), dense: d })
            .filter(|o| o.area.abs() > 1e-9)
            .collect());
    }
    Err("the smoothed outline crosses itself".into())
}

/// The unit-step loop as one point per cell edge, plus the turning corners that lie on the plate's
/// edge. A point on the plate's edge is pinned. The others sit where the density crosses the
/// threshold between the solid and the void cell (as in marching squares), not at the edge's middle.
fn dense(corners: &[(i64, i64)], mask: &Mask) -> Vec<(Point, bool)> {
    let (nx, ny) = (mask.nx as f64, mask.ny as f64);
    let n = corners.len();
    let on_edge = |p: Point| p.0 == 0.0 || p.1 == 0.0 || p.0 == nx || p.1 == ny;
    let mut out = Vec::with_capacity(2 * n);
    for k in 0..n {
        let (Some(a), Some(b), Some(c)) = (corners.get(k), corners.get((k + 1) % n), corners.get((k + 2) % n)) else { continue };
        let mid = ((a.0 + b.0) as f64 / 2.0, (a.1 + b.1) as f64 / 2.0);
        if on_edge(mid) {
            out.push((mid, true));
        } else {
            out.push((crossing(mid, (b.0 - a.0, b.1 - a.1), mask), false));
        }
        let corner = (b.0 as f64, b.1 as f64);
        let turns = (b.0 - a.0, b.1 - a.1) != (c.0 - b.0, c.1 - b.1);
        if turns && on_edge(corner) {
            out.push((corner, true));
        }
    }
    out
}

/// Where the density crosses the threshold on the line between the centres of the solid cell (left
/// of the edge) and the void cell (right of it), as an offset of at most half a cell from the edge's
/// middle. Where the mask was edited (dropped islands, closed diagonal contacts) the density does not
/// match and the middle is used.
fn crossing(mid: Point, dir: (i64, i64), mask: &Mask) -> Point {
    let normal = (-(dir.1 as f64), dir.0 as f64); // left of the edge, into the solid
    let cell = |s: f64| ((mid.0 + s * normal.0).floor() as i64, (mid.1 + s * normal.1).floor() as i64);
    let (solid, void) = (cell(0.5), cell(-0.5));
    let (ds, dv) = (mask.density_at(solid.0, solid.1), mask.density_at(void.0, void.1));
    if !(ds > 0.5 && dv <= 0.5) {
        return mid;
    }
    let t = ((ds - 0.5) / (ds - dv)).clamp(0.0, 1.0);
    (mid.0 + (0.5 - t) * normal.0, mid.1 + (0.5 - t) * normal.1)
}

fn taubin(points: &[(Point, bool)], passes: usize) -> Vec<(Point, bool)> {
    let n = points.len();
    let mut cur = points.to_vec();
    if n < 3 {
        return cur;
    }
    for step in 0..2 * passes {
        let f = if step % 2 == 0 { LAMBDA } else { MU };
        let prev = cur.clone();
        for (k, slot) in cur.iter_mut().enumerate() {
            let (Some(&(a, _)), Some(&(p, pinned)), Some(&(b, _))) = (prev.get((k + n - 1) % n), prev.get(k), prev.get((k + 1) % n)) else { continue };
            if !pinned {
                slot.0 = (p.0 + f * ((a.0 + b.0) / 2.0 - p.0), p.1 + f * ((a.1 + b.1) / 2.0 - p.1));
            }
        }
    }
    cur
}

/// Splits a smoothed loop into straight runs (between pinned points) and curves (through the free
/// points between two pinned ones).
fn shape(points: &[(Point, bool)]) -> Shape {
    let n = points.len();
    let Some(start) = points.iter().position(|&(_, pinned)| pinned) else {
        return Shape::Spline(resample(&points.iter().map(|&(p, _)| p).collect::<Vec<_>>(), true));
    };
    if points.iter().all(|&(_, pinned)| pinned) {
        let poly: Vec<Point> = points.iter().map(|&(p, _)| p).collect();
        return Shape::Polygon(drop_collinear(&poly, true));
    }
    // Walk once round from a pinned point, back to it.
    let at = |k: usize| points.get((start + k) % n).copied().unwrap_or(((0.0, 0.0), true));
    let mut pieces: Vec<Piece> = Vec::new();
    let mut k = 0;
    while k < n {
        let (p, _) = at(k);
        if at(k + 1).1 {
            // A straight step to the next pinned point; extend the last line if there is one.
            let q = at(k + 1).0;
            match pieces.last_mut() {
                Some(Piece::Line(v)) => v.push(q),
                _ => pieces.push(Piece::Line(vec![p, q])),
            }
            k += 1;
        } else {
            let mut run = vec![p];
            let mut m = k + 1;
            while m < n && !at(m).1 {
                run.push(at(m).0);
                m += 1;
            }
            run.push(at(m).0);
            pieces.push(Piece::Curve(resample(&run, false)));
            k = m;
        }
    }
    for piece in &mut pieces {
        if let Piece::Line(v) = piece {
            *v = drop_collinear(v, false);
        }
    }
    Shape::Mixed(pieces)
}

/// Every `CURVE_STEP`-th point (more apart when the curve would pass SolidCraft's limit); an open
/// run keeps both ends exactly.
fn resample(points: &[Point], closed: bool) -> Vec<Point> {
    let n = points.len();
    let step = CURVE_STEP.max(n.div_ceil(MAX_CURVE_POINTS - 2));
    let mut out: Vec<Point> = points.iter().step_by(step).copied().collect();
    if closed {
        if out.len() < 3 {
            return points.iter().take(MAX_CURVE_POINTS).copied().collect();
        }
    } else if let Some(&last) = points.last() {
        // Drop a sample that sits too close to the end, then end exactly on it.
        if n > 1 && (n - 1) % step != 0 && (n - 1) % step < step / 2 + 1 && out.len() > 1 {
            out.pop();
        }
        if out.last() != Some(&last) {
            out.push(last);
        }
    }
    out
}

fn drop_collinear(points: &[Point], closed: bool) -> Vec<Point> {
    let n = points.len();
    if n < 3 {
        return points.to_vec();
    }
    let mut out = Vec::with_capacity(n);
    for k in 0..n {
        let Some(&q) = points.get(k) else { continue };
        let ends = !closed && (k == 0 || k == n - 1);
        let (Some(p), Some(r)) = (points.get((k + n - 1) % n), points.get((k + 1) % n)) else { continue };
        let cross = (q.0 - p.0) * (r.1 - q.1) - (q.1 - p.1) * (r.0 - q.0);
        let dot = (q.0 - p.0) * (r.0 - q.0) + (q.1 - p.1) * (r.1 - q.1);
        if ends || cross.abs() > 1e-9 || dot < 0.0 {
            out.push(q);
        }
    }
    out
}

/// Whether any two segments of the loops cross (within a loop or between loops). Segments are at
/// most about a cell long, so a grid of 1-cell buckets finds every candidate pair.
fn any_crossing(loops: &[Vec<Point>]) -> bool {
    let mut buckets: HashMap<(i64, i64), Vec<(usize, usize)>> = HashMap::new();
    for (l, pts) in loops.iter().enumerate() {
        for k in 0..pts.len() {
            if let (Some(a), Some(b)) = (pts.get(k), pts.get((k + 1) % pts.len())) {
                let key = (((a.0 + b.0) / 2.0).floor() as i64, ((a.1 + b.1) / 2.0).floor() as i64);
                buckets.entry(key).or_default().push((l, k));
            }
        }
    }
    let seg = |(l, k): (usize, usize)| -> Option<(Point, Point)> {
        let pts = loops.get(l)?;
        Some((*pts.get(k)?, *pts.get((k + 1) % pts.len())?))
    };
    for (&(bx, by), list) in &buckets {
        for &s in list {
            for dx in -1..=1 {
                for dy in -1..=1 {
                    let Some(other) = buckets.get(&(bx + dx, by + dy)) else { continue };
                    for &t in other {
                        if t <= s {
                            continue;
                        }
                        let n = loops.get(s.0).map_or(0, Vec::len);
                        let neighbours = s.0 == t.0 && (t.1 == (s.1 + 1) % n || s.1 == (t.1 + 1) % n);
                        if neighbours {
                            continue;
                        }
                        if let (Some((p, q)), Some((r, u))) = (seg(s), seg(t))
                            && crosses(p, q, r, u)
                        {
                            return true;
                        }
                    }
                }
            }
        }
    }
    false
}

/// Proper crossing or touching of two segments (shared end points of neighbours are excluded by the caller).
fn crosses(p: Point, q: Point, r: Point, s: Point) -> bool {
    let o = |a: Point, b: Point, c: Point| (b.0 - a.0) * (c.1 - a.1) - (b.1 - a.1) * (c.0 - a.0);
    let (d1, d2, d3, d4) = (o(p, q, r), o(p, q, s), o(r, s, p), o(r, s, q));
    d1 * d2 < 1e-12 && d3 * d4 < 1e-12 && {
        // Collinear segments: only overlapping ones count.
        let (minx, maxx) = (p.0.min(q.0), p.0.max(q.0));
        let (miny, maxy) = (p.1.min(q.1), p.1.max(q.1));
        let (rminx, rmaxx) = (r.0.min(s.0), r.0.max(s.0));
        let (rminy, rmaxy) = (r.1.min(s.1), r.1.max(s.1));
        minx <= rmaxx + 1e-12 && rminx <= maxx + 1e-12 && miny <= rmaxy + 1e-12 && rminy <= maxy + 1e-12
    }
}

/// Signed area of a loop (positive counter-clockwise).
pub fn signed_area(points: &[Point]) -> f64 {
    let n = points.len();
    let mut s = 0.0;
    for k in 0..n {
        if let (Some(a), Some(b)) = (points.get(k), points.get((k + 1) % n)) {
            s += a.0 * b.1 - b.0 * a.1;
        }
    }
    s / 2.0
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

    fn total_area(m: &Mask) -> f64 {
        smooth_outlines(m).unwrap().iter().map(|o| o.area).sum()
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
    fn a_full_plate_is_one_ccw_polygon_of_four_points() {
        let m = mask(&["######", "######", "######", "######"]);
        let l = smooth_outlines(&m).unwrap();
        assert_eq!(l.len(), 1);
        let Shape::Polygon(p) = &l[0].shape else { panic!("{:?}", l[0].shape) };
        assert_eq!(p.len(), 4, "{p:?}");
        assert!((l[0].area - 24.0).abs() < 1e-9, "the plate's corners stay exact");
    }

    #[test]
    fn a_hole_is_a_clockwise_spline() {
        let m = mask(&["########", "##....##", "##....##", "##....##", "########"]);
        let l = smooth_outlines(&m).unwrap();
        assert_eq!(l.len(), 2);
        let hole = l.iter().find(|o| o.area < 0.0).unwrap();
        assert!(matches!(hole.shape, Shape::Spline(_)), "{:?}", hole.shape);
        let net = total_area(&m);
        assert!((net - m.count() as f64).abs() < 1.5, "net {net} vs {} cells", m.count());
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

    /// A slanted bar from the left edge into the plate: the part on the left edge stays a straight
    /// line exactly on x = 0, the rest becomes curves, and the area stays close to the cell count.
    #[test]
    fn a_slanted_bar_is_straight_on_the_edge_and_curved_inside() {
        let n = 24;
        let mut rows = Vec::new();
        for r in 0..n {
            let mut line = String::new();
            for c in 0..n {
                let d = c as i64 - (n as i64 - 1 - r as i64) / 2;
                line.push(if (0..4).contains(&d) || (c == 0 && r >= n / 2) { '#' } else { '.' });
            }
            rows.push(line);
        }
        let refs: Vec<&str> = rows.iter().map(String::as_str).collect();
        let mut m = mask(&refs);
        m.keep_connected_to_left();
        m.close_diagonal_contacts().unwrap();
        let cells = m.count() as f64;
        let l = smooth_outlines(&m).unwrap();
        assert_eq!(l.len(), 1);
        let Shape::Mixed(pieces) = &l[0].shape else { panic!("{:?}", l[0].shape) };
        assert!(pieces.iter().any(|p| matches!(p, Piece::Curve(_))));
        let on_left = pieces.iter().any(|p| matches!(p, Piece::Line(v) if v.windows(2).any(|w| w[0].0 == 0.0 && w[1].0 == 0.0)));
        assert!(on_left, "{pieces:?}");
        // Pieces join end to end and close.
        let ends: Vec<(Point, Point)> = pieces.iter().map(|p| match p { Piece::Line(v) | Piece::Curve(v) => (v[0], *v.last().unwrap()) }).collect();
        for k in 0..ends.len() {
            assert_eq!(ends[k].1, ends[(k + 1) % ends.len()].0, "piece {k} does not join the next");
        }
        assert!(l[0].drawn_points() < l[0].dense.len() / 2);
        assert!((l[0].area - cells).abs() / cells < 0.06, "area {} vs {cells} cells", l[0].area);
    }

    /// A disc inside the plate: one closed spline, much rounder than the stair steps.
    #[test]
    fn a_disc_becomes_round() {
        let (nx, ny, r) = (40usize, 40usize, 12.0);
        let mut d = vec![0.0; nx * ny];
        for i in 0..nx {
            for j in 0..ny {
                let (x, y) = (i as f64 + 0.5 - 20.0, j as f64 + 0.5 - 20.0);
                if x.hypot(y) <= r {
                    d[i * ny + j] = 1.0;
                }
            }
        }
        let m = Mask::from_density(&d, nx, ny, 0.5).unwrap();
        let l = smooth_outlines(&m).unwrap();
        assert_eq!(l.len(), 1);
        assert!(matches!(l[0].shape, Shape::Spline(_)));
        let spread = |pts: &[Point]| {
            let radii: Vec<f64> = pts.iter().map(|p| (p.0 - 20.0).hypot(p.1 - 20.0)).collect();
            radii.iter().cloned().fold(f64::MIN, f64::max) - radii.iter().cloned().fold(f64::MAX, f64::min)
        };
        let raw: Vec<Point> = dense(&m.trace().unwrap()[0], &m).iter().map(|&(p, _)| p).collect();
        let (before, after) = (spread(&raw), spread(&l[0].dense));
        assert!(after < 0.7 * before, "radius spread {after} (stair steps {before})");
        assert!((l[0].area - m.count() as f64).abs() / (m.count() as f64) < 0.02, "Taubin keeps the area");

        // With a graded density, as the optimiser's filter gives, the outline follows the 0.5 level
        // between the cells and comes out nearly a true circle.
        let mut g = vec![0.0; nx * ny];
        for i in 0..nx {
            for j in 0..ny {
                let (x, y) = (i as f64 + 0.5 - 20.0, j as f64 + 0.5 - 20.0);
                g[i * ny + j] = (0.5 + (r - x.hypot(y)) / 3.0).clamp(0.0, 1.0);
            }
        }
        let m = Mask::from_density(&g, nx, ny, 0.5).unwrap();
        let l = smooth_outlines(&m).unwrap();
        let s = spread(&l[0].dense);
        assert!(s < 0.1, "radius spread {s} with a graded density");
        let true_area = std::f64::consts::PI * r * r;
        assert!((l[0].area - true_area).abs() / true_area < 0.01, "area {} vs {true_area}", l[0].area);
    }

    #[test]
    fn smoothed_loops_never_cross_on_random_blobs() {
        // A deterministic pseudo-random field: thin features everywhere, so the gentler fallbacks run too.
        let (nx, ny) = (30usize, 20usize);
        let mut s = 12345u64;
        let mut d = vec![0.0; nx * ny];
        for v in &mut d {
            s = s.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            *v = if (s >> 33) % 100 < 60 { 1.0 } else { 0.0 };
        }
        let mut m = Mask::from_density(&d, nx, ny, 0.5).unwrap();
        m.close_diagonal_contacts().unwrap();
        let loops = smooth_outlines(&m).unwrap();
        let dense: Vec<Vec<Point>> = loops.iter().map(|o| o.dense.clone()).collect();
        assert!(!any_crossing(&dense));
        for o in &loops {
            if let Shape::Spline(p) | Shape::Polygon(p) = &o.shape {
                assert!(p.len() <= MAX_CURVE_POINTS);
            }
        }
    }

    #[test]
    fn crossing_finds_a_figure_eight() {
        let eight = vec![(0.0, 0.0), (1.0, 1.0), (1.0, 0.0), (0.0, 1.0)];
        assert!(any_crossing(&[eight]));
        let square = vec![(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)];
        assert!(!any_crossing(&[square]));
    }

    #[test]
    fn long_curves_stay_under_the_spline_limit() {
        let run: Vec<Point> = (0..3000).map(|k| (k as f64 * 0.5, 0.0)).collect();
        let r = resample(&run, false);
        assert!(r.len() <= MAX_CURVE_POINTS, "{}", r.len());
        assert_eq!(r[0], run[0]);
        assert_eq!(*r.last().unwrap(), *run.last().unwrap());
    }

    #[test]
    fn total_area_is_close_to_the_cell_count() {
        let m = mask(&["######", "##..##", "##..##", "######"]);
        let a = total_area(&m);
        assert!((a - m.count() as f64).abs() < 2.0, "{a}");
    }
}
