#!/usr/bin/env python3
"""One-layer connected Fermat prototype for strict continuous slicing.

This is intentionally standalone and dependency-free. It is not the production
OrcaSlicer implementation. The goal is to make the geometry easy to debug:

* polygon signed-distance field
* marching-squares offset contours
* one continuous extrusion polyline
* JSON metrics and SVG preview output

The production slicer should replace the sampled SDF/marching-squares contour
generator with Orca's Clipper offset geometry, while keeping the same invariants:
one path, exact endpoint continuity, and no travel-like gaps.
"""

from __future__ import annotations

import argparse
import struct
import heapq
import json
import math
import sys
import time
import zlib
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence


Point = tuple[float, float]
EPS = 1e-9


def add(a: Point, b: Point) -> Point:
    return (a[0] + b[0], a[1] + b[1])


def sub(a: Point, b: Point) -> Point:
    return (a[0] - b[0], a[1] - b[1])


def mul(a: Point, s: float) -> Point:
    return (a[0] * s, a[1] * s)


def dot(a: Point, b: Point) -> float:
    return a[0] * b[0] + a[1] * b[1]


def dist2(a: Point, b: Point) -> float:
    dx = a[0] - b[0]
    dy = a[1] - b[1]
    return dx * dx + dy * dy


def dist(a: Point, b: Point) -> float:
    return math.sqrt(dist2(a, b))


def lerp(a: Point, b: Point, t: float) -> Point:
    return (a[0] + (b[0] - a[0]) * t, a[1] + (b[1] - a[1]) * t)


def point_segment_distance(p: Point, a: Point, b: Point) -> float:
    ab = sub(b, a)
    denom = dot(ab, ab)
    if denom <= EPS:
        return dist(p, a)
    t = max(0.0, min(1.0, dot(sub(p, a), ab) / denom))
    q = lerp(a, b, t)
    return dist(p, q)


def orient(a: Point, b: Point, c: Point) -> float:
    return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])


def on_segment(a: Point, b: Point, p: Point, eps: float = 1e-7) -> bool:
    return (
        min(a[0], b[0]) - eps <= p[0] <= max(a[0], b[0]) + eps
        and min(a[1], b[1]) - eps <= p[1] <= max(a[1], b[1]) + eps
        and abs(orient(a, b, p)) <= eps
    )


def segments_intersect(a: Point, b: Point, c: Point, d: Point, eps: float = 1e-7) -> bool:
    o1 = orient(a, b, c)
    o2 = orient(a, b, d)
    o3 = orient(c, d, a)
    o4 = orient(c, d, b)

    if ((o1 > eps and o2 < -eps) or (o1 < -eps and o2 > eps)) and (
        (o3 > eps and o4 < -eps) or (o3 < -eps and o4 > eps)
    ):
        return True

    return on_segment(a, b, c, eps) or on_segment(a, b, d, eps) or on_segment(c, d, a, eps) or on_segment(c, d, b, eps)


def segment_distance(a: Point, b: Point, c: Point, d: Point) -> float:
    if segments_intersect(a, b, c, d):
        return 0.0
    return min(
        point_segment_distance(a, c, d),
        point_segment_distance(b, c, d),
        point_segment_distance(c, a, b),
        point_segment_distance(d, a, b),
    )


def point_in_polygon(p: Point, poly: Sequence[Point]) -> bool:
    x, y = p
    inside = False
    n = len(poly)
    if n < 3:
        return False
    j = n - 1
    for i in range(n):
        xi, yi = poly[i]
        xj, yj = poly[j]
        if (yi > y) != (yj > y):
            x_at_y = (xj - xi) * (y - yi) / (yj - yi + EPS) + xi
            if x < x_at_y:
                inside = not inside
        j = i
    return inside


def signed_area(poly: Sequence[Point]) -> float:
    area = 0.0
    for i, p in enumerate(poly):
        q = poly[(i + 1) % len(poly)]
        area += p[0] * q[1] - q[0] * p[1]
    return 0.5 * area


def polygon_centroid(poly: Sequence[Point]) -> Point:
    area2 = 0.0
    cx = 0.0
    cy = 0.0
    for i, p in enumerate(poly):
        q = poly[(i + 1) % len(poly)]
        cross = p[0] * q[1] - q[0] * p[1]
        area2 += cross
        cx += (p[0] + q[0]) * cross
        cy += (p[1] + q[1]) * cross
    if abs(area2) <= EPS:
        return (
            sum(p[0] for p in poly) / max(1, len(poly)),
            sum(p[1] for p in poly) / max(1, len(poly)),
        )
    return (cx / (3.0 * area2), cy / (3.0 * area2))


def polyline_length(points: Sequence[Point], closed: bool = False) -> float:
    if len(points) < 2:
        return 0.0
    total = 0.0
    for i in range(len(points) - 1):
        total += dist(points[i], points[i + 1])
    if closed:
        total += dist(points[-1], points[0])
    return total


def point_at_closed_fraction(points: Sequence[Point], fraction: float) -> Point:
    if not points:
        return (0.0, 0.0)
    if len(points) == 1:
        return points[0]
    fraction = fraction % 1.0
    total = polyline_length(points, closed=True)
    if total <= EPS:
        return points[0]
    target = total * fraction
    walked = 0.0
    for i, p in enumerate(points):
        q = points[(i + 1) % len(points)]
        length = dist(p, q)
        if walked + length >= target:
            return lerp(p, q, (target - walked) / max(length, EPS))
        walked += length
    return points[-1]


def distance_to_loop(p: Point, loop: Sequence[Point]) -> float:
    best = float("inf")
    for i, a in enumerate(loop):
        b = loop[(i + 1) % len(loop)]
        best = min(best, point_segment_distance(p, a, b))
    return best


def nearest_index(points: Sequence[Point], p: Point) -> int:
    best_idx = 0
    best_dist = float("inf")
    for i, q in enumerate(points):
        d = dist2(p, q)
        if d < best_dist:
            best_idx = i
            best_dist = d
    return best_idx


def circle_polygon(radius: float, segments: int, center: Point = (0.0, 0.0)) -> list[Point]:
    cx, cy = center
    return [
        (cx + math.cos(2.0 * math.pi * i / segments) * radius, cy + math.sin(2.0 * math.pi * i / segments) * radius)
        for i in range(segments)
    ]


def star_polygon(outer_radius: float, inner_radius: float, points: int) -> list[Point]:
    out: list[Point] = []
    for i in range(points * 2):
        r = outer_radius if i % 2 == 0 else inner_radius
        a = -math.pi / 2.0 + math.pi * i / points
        out.append((math.cos(a) * r, math.sin(a) * r))
    return out


@dataclass(frozen=True)
class PolygonModel:
    name: str
    outer: list[Point]
    holes: list[list[Point]]

    def bounds(self, margin: float = 0.0) -> tuple[float, float, float, float]:
        pts = list(self.outer)
        for hole in self.holes:
            pts.extend(hole)
        min_x = min(p[0] for p in pts) - margin
        min_y = min(p[1] for p in pts) - margin
        max_x = max(p[0] for p in pts) + margin
        max_y = max(p[1] for p in pts) + margin
        return (min_x, min_y, max_x, max_y)

    def contains(self, p: Point) -> bool:
        return point_in_polygon(p, self.outer) and not any(point_in_polygon(p, hole) for hole in self.holes)

    def sdf(self, p: Point) -> float:
        boundary_distance = distance_to_loop(p, self.outer)
        for hole in self.holes:
            boundary_distance = min(boundary_distance, distance_to_loop(p, hole))
        return boundary_distance if self.contains(p) else -boundary_distance


def available_shapes() -> dict[str, PolygonModel]:
    return {
        "rectangle": PolygonModel(
            "rectangle",
            [(-45.0, -28.0), (45.0, -28.0), (45.0, 28.0), (-45.0, 28.0)],
            [],
        ),
        "dumbbell": PolygonModel(
            "dumbbell",
            [
                (-52.0, -25.0),
                (-20.0, -25.0),
                (-20.0, -8.0),
                (20.0, -8.0),
                (20.0, -25.0),
                (52.0, -25.0),
                (52.0, 25.0),
                (20.0, 25.0),
                (20.0, 8.0),
                (-20.0, 8.0),
                (-20.0, 25.0),
                (-52.0, 25.0),
            ],
            [],
        ),
        "c_shape": PolygonModel(
            "c_shape",
            [
                (-48.0, -32.0),
                (48.0, -32.0),
                (48.0, -16.0),
                (-18.0, -16.0),
                (-18.0, 16.0),
                (48.0, 16.0),
                (48.0, 32.0),
                (-48.0, 32.0),
            ],
            [],
        ),
        "annulus": PolygonModel(
            "annulus",
            circle_polygon(42.0, 160),
            [circle_polygon(16.0, 96)],
        ),
        "square_hole": PolygonModel(
            "square_hole",
            [(-45.0, -36.0), (45.0, -36.0), (45.0, 36.0), (-45.0, 36.0)],
            [[(-12.0, -12.0), (-12.0, 12.0), (12.0, 12.0), (12.0, -12.0)]],
        ),
        "two_holes": PolygonModel(
            "two_holes",
            [(-56.0, -34.0), (56.0, -34.0), (56.0, 34.0), (-56.0, 34.0)],
            [circle_polygon(10.0, 64, (-21.0, 0.0)), circle_polygon(10.0, 64, (21.0, 0.0))],
        ),
        "star": PolygonModel("star", star_polygon(45.0, 22.0, 7), []),
    }


@dataclass
class SDFGrid:
    model: PolygonModel
    min_x: float
    min_y: float
    max_x: float
    max_y: float
    nx: int
    ny: int
    values: list[float]

    @property
    def dx(self) -> float:
        return (self.max_x - self.min_x) / self.nx

    @property
    def dy(self) -> float:
        return (self.max_y - self.min_y) / self.ny

    @property
    def step(self) -> float:
        return max(self.dx, self.dy)

    def value(self, i: int, j: int) -> float:
        return self.values[j * (self.nx + 1) + i]

    def point(self, i: int, j: int) -> Point:
        return (self.min_x + i * self.dx, self.min_y + j * self.dy)

    def max_sdf(self) -> float:
        return max(self.values)


def build_sdf_grid(model: PolygonModel, grid_cells: int, margin: float) -> SDFGrid:
    min_x, min_y, max_x, max_y = model.bounds(margin)
    width = max_x - min_x
    height = max_y - min_y
    if width >= height:
        nx = max(24, grid_cells)
        ny = max(24, int(round(grid_cells * height / width)))
    else:
        ny = max(24, grid_cells)
        nx = max(24, int(round(grid_cells * width / height)))

    values: list[float] = []
    for j in range(ny + 1):
        y = min_y + j * height / ny
        for i in range(nx + 1):
            x = min_x + i * width / nx
            values.append(model.sdf((x, y)))

    return SDFGrid(model, min_x, min_y, max_x, max_y, nx, ny, values)


EDGE_CORNERS = {
    0: (0, 1),  # bottom
    1: (1, 2),  # right
    2: (2, 3),  # top
    3: (3, 0),  # left
}

CASE_SEGMENTS: dict[int, list[tuple[int, int]]] = {
    0: [],
    1: [(3, 0)],
    2: [(0, 1)],
    3: [(3, 1)],
    4: [(1, 2)],
    6: [(0, 2)],
    7: [(3, 2)],
    8: [(2, 3)],
    9: [(0, 2)],
    11: [(1, 2)],
    12: [(1, 3)],
    13: [(0, 1)],
    14: [(3, 0)],
    15: [],
}


@dataclass
class ContourLoop:
    level_index: int
    offset: float
    points: list[Point]
    area: float
    length: float
    centroid: Point


def interpolate_edge(corners: Sequence[Point], values: Sequence[float], level: float, edge: int) -> Point:
    a_idx, b_idx = EDGE_CORNERS[edge]
    va = values[a_idx]
    vb = values[b_idx]
    if abs(vb - va) <= EPS:
        t = 0.5
    else:
        t = (level - va) / (vb - va)
    t = max(0.0, min(1.0, t))
    return lerp(corners[a_idx], corners[b_idx], t)


def marching_cell_segments(grid: SDFGrid, i: int, j: int, level: float) -> list[tuple[Point, Point]]:
    corners = [
        grid.point(i, j),
        grid.point(i + 1, j),
        grid.point(i + 1, j + 1),
        grid.point(i, j + 1),
    ]
    values = [
        grid.value(i, j),
        grid.value(i + 1, j),
        grid.value(i + 1, j + 1),
        grid.value(i, j + 1),
    ]

    case = 0
    for bit, value in enumerate(values):
        if value >= level:
            case |= 1 << bit

    if case == 5:
        center = sum(values) / 4.0
        pairs = [(3, 2), (0, 1)] if center >= level else [(3, 0), (1, 2)]
    elif case == 10:
        center = sum(values) / 4.0
        pairs = [(0, 3), (1, 2)] if center >= level else [(0, 1), (2, 3)]
    else:
        pairs = CASE_SEGMENTS.get(case, [])

    return [(interpolate_edge(corners, values, level, a), interpolate_edge(corners, values, level, b)) for a, b in pairs]


def assemble_loops(segments: list[tuple[Point, Point]], quant: float, min_length: float) -> list[list[Point]]:
    def key(p: Point) -> tuple[int, int]:
        return (int(round(p[0] / quant)), int(round(p[1] / quant)))

    adjacency: dict[tuple[int, int], list[int]] = {}
    for idx, (a, b) in enumerate(segments):
        adjacency.setdefault(key(a), []).append(idx)
        adjacency.setdefault(key(b), []).append(idx)

    used = [False] * len(segments)
    loops: list[list[Point]] = []

    for start_idx in range(len(segments)):
        if used[start_idx]:
            continue

        a, b = segments[start_idx]
        used[start_idx] = True
        start_key = key(a)
        current_key = key(b)
        current_point = b
        path = [a, b]

        guard = 0
        while current_key != start_key and guard < len(segments) + 4:
            guard += 1
            next_idx = None
            for candidate in adjacency.get(current_key, []):
                if not used[candidate]:
                    next_idx = candidate
                    break
            if next_idx is None:
                break

            p, q = segments[next_idx]
            used[next_idx] = True
            if key(p) == current_key:
                current_point = q
            else:
                current_point = p
            path.append(current_point)
            current_key = key(current_point)

        if len(path) >= 4 and dist(path[-1], path[0]) <= quant * 4.0:
            if dist(path[-1], path[0]) > EPS:
                path[-1] = path[0]
            open_loop = path[:-1]
            length = polyline_length(open_loop, closed=True)
            if length >= min_length and abs(signed_area(open_loop)) >= quant * quant:
                loops.append(open_loop)

    return loops


def resample_closed_loop(points: Sequence[Point], target_step: float) -> list[Point]:
    length = polyline_length(points, closed=True)
    if length <= EPS:
        return list(points)

    sample_count = max(16, int(math.ceil(length / max(target_step, EPS))))
    out: list[Point] = []
    segment_lengths: list[float] = []
    for i, p in enumerate(points):
        segment_lengths.append(dist(p, points[(i + 1) % len(points)]))

    seg_idx = 0
    seg_start_distance = 0.0
    for sample_idx in range(sample_count):
        target = length * sample_idx / sample_count
        while seg_start_distance + segment_lengths[seg_idx] < target and seg_idx < len(points) - 1:
            seg_start_distance += segment_lengths[seg_idx]
            seg_idx += 1
        seg_len = max(segment_lengths[seg_idx], EPS)
        t = (target - seg_start_distance) / seg_len
        out.append(lerp(points[seg_idx], points[(seg_idx + 1) % len(points)], t))

    return out


def extract_contours(grid: SDFGrid, level_index: int, offset: float, spacing: float) -> list[ContourLoop]:
    segments: list[tuple[Point, Point]] = []
    for j in range(grid.ny):
        for i in range(grid.nx):
            segments.extend(marching_cell_segments(grid, i, j, offset))

    raw_loops = assemble_loops(segments, quant=grid.step * 0.01, min_length=spacing * 2.5)
    loops: list[ContourLoop] = []
    for raw in raw_loops:
        points = resample_closed_loop(raw, target_step=max(spacing * 0.45, grid.step * 0.75))
        area = signed_area(points)
        if abs(area) < spacing * spacing:
            continue
        length = polyline_length(points, closed=True)
        loops.append(
            ContourLoop(
                level_index=level_index,
                offset=offset,
                points=points if area >= 0.0 else list(reversed(points)),
                area=abs(area),
                length=length,
                centroid=polygon_centroid(points),
            )
        )

    loops.sort(key=lambda loop: loop.area, reverse=True)
    return loops


def generate_offset_contours(
    grid: SDFGrid,
    line_width: float,
    spacing: float,
    max_levels: int,
) -> list[ContourLoop]:
    contours: list[ContourLoop] = []
    max_sdf = grid.max_sdf()
    first_offset = line_width * 0.5
    if max_sdf < first_offset:
        return contours

    offset = first_offset
    level_index = 0
    while offset <= max_sdf - spacing * 0.15 and level_index < max_levels:
        loops = extract_contours(grid, level_index, offset, spacing)
        contours.extend(loops)
        offset += spacing
        level_index += 1

    return contours


def append_point(path: list[Point], p: Point) -> None:
    if not path or dist(path[-1], p) > EPS:
        path.append(p)


def append_points(path: list[Point], points: Iterable[Point]) -> None:
    for p in points:
        append_point(path, p)


def closed_arc(points: Sequence[Point], start: int, end: int, direction: int) -> list[Point]:
    if not points:
        return []
    n = len(points)
    out = [points[start % n]]
    idx = start % n
    guard = 0
    while idx != end % n and guard <= n + 2:
        idx = (idx + direction) % n
        out.append(points[idx])
        guard += 1
    return out


def almost_full_loop(points: Sequence[Point], start: int, direction: int, gap_points: int) -> list[Point]:
    n = len(points)
    gap_points = max(1, min(gap_points, max(1, n // 3)))
    end = (start - direction * gap_points) % n
    return closed_arc(points, start, end, direction)


def full_loop(points: Sequence[Point], start: int, direction: int) -> list[Point]:
    n = len(points)
    out = [points[start % n]]
    idx = start % n
    for _ in range(n):
        idx = (idx + direction) % n
        out.append(points[idx])
    return out


def loop_arc_between_params(points: Sequence[Point], start: float, end: float, direction: int) -> list[Point]:
    if not points:
        return []

    n = len(points)
    start %= n
    end %= n
    out = [loop_point_at_index(points, start)]

    if direction >= 0:
        i = math.ceil(start) % n
        if abs(float(i) - start) <= 1e-6:
            i = (i + 1) % n
        guard = 0
        while subset_cycle_param(start, end, float(i), False, False) and guard <= n + 2:
            append_point(out, points[i])
            i = (i + 1) % n
            guard += 1
    else:
        i = math.floor(start) % n
        if abs(float(i) - start) <= 1e-6:
            i = (i + n - 1) % n
        guard = 0
        while subset_cycle_param(end, start, float(i), False, False) and guard <= n + 2:
            append_point(out, points[i])
            i = (i + n - 1) % n
            guard += 1

    append_loop_param(out, points, end)
    return out


def loops_by_level(contours: Sequence[ContourLoop]) -> dict[int, list[ContourLoop]]:
    grouped: dict[int, list[ContourLoop]] = {}
    for loop in contours:
        grouped.setdefault(loop.level_index, []).append(loop)
    return grouped


def filter_printable_contours(contours: Sequence[ContourLoop], spacing: float) -> list[ContourLoop]:
    min_area = spacing * spacing * 8.0
    min_length = spacing * 5.0
    return [loop for loop in contours if loop.area >= min_area and loop.length >= min_length]


def centroid_spread(loops: Sequence[ContourLoop]) -> float:
    spread = 0.0
    for i, a in enumerate(loops):
        for b in loops[i + 1 :]:
            spread = max(spread, dist(a.centroid, b.centroid))
    return spread


def filter_topology_stable_levels(contours: Sequence[ContourLoop], spacing: float) -> list[ContourLoop]:
    grouped = loops_by_level(contours)
    if not grouped:
        return []

    first_level = min(grouped)
    initial_count = len(grouped[first_level])
    initial_spread = centroid_spread(grouped[first_level])
    max_allowed_spread_jump = spacing * 8.0
    stable: list[ContourLoop] = []

    for level in sorted(grouped):
        loops = grouped[level]
        spread = centroid_spread(loops)
        excessive_fragmentation = len(loops) > max(initial_count * 3, initial_count + 4)
        medial_split = level != first_level and len(loops) <= initial_count and spread > initial_spread + max_allowed_spread_jump
        if excessive_fragmentation or medial_split:
            break
        stable.extend(loops)

    return stable


def sampled_loop_gap(a: ContourLoop, b: ContourLoop, max_samples: int = 128) -> float:
    """Estimate the minimum distance between two closed contours."""

    best = float("inf")
    stride_a = max(1, len(a.points) // max_samples)
    stride_b = max(1, len(b.points) // max_samples)

    for p in a.points[::stride_a]:
        best = min(best, distance_to_loop(p, b.points))
    for p in b.points[::stride_b]:
        best = min(best, distance_to_loop(p, a.points))

    return best


def filter_ring_medial_overlap(contours: Sequence[ContourLoop], spacing: float) -> list[ContourLoop]:
    """Drop one-hole ring levels whose opposing fronts are closer than spacing.

    Marching squares may emit one extra contour pair near the medial axis of a
    ring. Those two centerlines are less than one line spacing apart, so no
    Fermat ordering can route both without violating the equal-spacing check.
    """

    grouped = loops_by_level(contours)
    stable: list[ContourLoop] = []

    for level in sorted(grouped):
        loops = grouped[level]
        if len(loops) == 2:
            outer, inner = sorted(loops, key=lambda loop: loop.area, reverse=True)
            if sampled_loop_gap(outer, inner) < spacing * 0.98:
                break
        stable.extend(loops)

    return stable if stable else list(contours)


def build_even_odd_fermat_chain(level_loops: Sequence[ContourLoop], anchor: Point, spacing: float) -> list[Point]:
    """Build a boundary-to-boundary Fermat-like contour spiral for one-loop levels.

    Even contour levels form the inward branch. Odd contour levels form the
    outward branch. This mirrors the useful CFS property that the path can start
    and end near the boundary instead of stopping at the medial axis.
    """

    path: list[Point] = []
    even = [i for i in range(len(level_loops)) if i % 2 == 0]
    odd = [i for i in range(len(level_loops)) if i % 2 == 1]

    direction = 1
    for level_pos, loop_idx in enumerate(even):
        loop = level_loops[loop_idx]
        start = nearest_index(loop.points, anchor if not path else path[-1])
        gap = max(2, int(round(len(loop.points) * 0.025)))
        segment = almost_full_loop(loop.points, start, direction, gap)
        append_points(path, segment)
        direction *= -1

    if odd:
        direction *= -1
        for loop_idx in reversed(odd):
            loop = level_loops[loop_idx]
            start = nearest_index(loop.points, path[-1] if path else anchor)
            gap = max(2, int(round(len(loop.points) * 0.025)))
            segment = almost_full_loop(loop.points, start, direction, gap)
            append_points(path, segment)
            direction *= -1

    return path


def seam_gap_points(loop: ContourLoop, spacing: float) -> int:
    if not loop.points:
        return 1
    average_step = max(loop.length / len(loop.points), EPS)
    return max(2, min(max(2, len(loop.points) // 4), int(math.ceil(spacing * 1.15 / average_step))))


def seam_gap_points_for_width(loop: ContourLoop, spacing: float, width_factor: float) -> int:
    if not loop.points:
        return 1
    average_step = max(loop.length / len(loop.points), EPS)
    return max(2, min(max(2, len(loop.points) // 3), int(math.ceil(spacing * width_factor / average_step))))


def build_open_contour_spiral(
    model: PolygonModel,
    grid: SDFGrid,
    ordered_loops: Sequence[ContourLoop],
    anchor: Point,
    line_width: float,
    spacing: float,
) -> tuple[list[Point], int]:
    """Build a non-crossing contour-parallel spiral by opening every loop.

    A closed contour loop is not printable as part of a single no-travel path
    without either retracing or crossing into the next loop. This routine leaves
    one seam gap in each contour and connects through that gap to the next
    contour. For simply nested contours this creates an ordinary contour spiral.
    """

    path: list[Point] = []
    current_point = anchor
    direction = 1
    unsafe_connectors = 0
    required_clearance = line_width * 0.5 - spacing * 0.15

    for loop_idx, loop in enumerate(ordered_loops):
        if not loop.points:
            continue

        if path:
            start_idx, _, safe = sampled_connection_to_loop(
                model,
                current_point,
                loop,
                required_clearance=required_clearance,
                spacing=spacing,
            )
            start_point = loop.points[start_idx]
            if not safe:
                connector = grid_astar_connector(grid, current_point, start_point, required_clearance=required_clearance)
                if connector is None:
                    unsafe_connectors += 1
                else:
                    append_points(path, connector)
        else:
            start_idx = nearest_index(loop.points, current_point)

        gap = seam_gap_points(loop, spacing)
        append_points(path, almost_full_loop(loop.points, start_idx, direction, gap))
        current_point = path[-1]

    return path, unsafe_connectors


def build_fixed_seam_open_contour_spiral(
    model: PolygonModel,
    grid: SDFGrid,
    ordered_loops: Sequence[ContourLoop],
    seam_target: Point,
    line_width: float,
    spacing: float,
) -> tuple[list[Point], int]:
    path: list[Point] = []
    unsafe_connectors = 0
    direction = 1
    required_clearance = line_width * 0.5 - spacing * 0.15

    for loop in ordered_loops:
        entry_idx = nearest_index(loop.points, seam_target)
        entry = loop.points[entry_idx]
        if path and not segment_is_inside(model, path[-1], entry, required_clearance, spacing * 0.5):
            connector = grid_astar_connector(grid, path[-1], entry, required_clearance=required_clearance)
            if connector is None:
                unsafe_connectors += 1
            else:
                append_points(path, connector)
        append_points(path, almost_full_loop(loop.points, entry_idx, direction, seam_gap_points_for_width(loop, spacing, 2.5)))

    return path, unsafe_connectors


def build_two_slot_contour_weave(
    model: PolygonModel,
    grid: SDFGrid,
    ordered_loops: Sequence[ContourLoop],
    start_anchor: Point,
    exit_anchor: Point,
    line_width: float,
    spacing: float,
) -> tuple[list[Point], int]:
    """Build a closed cyclic spiral from contour arcs.

    Every contour is split at two boundary-controlled slots. The inward pass
    consumes alternating arcs, then the outward pass consumes the complementary
    arcs. This produces a cyclic path whose first and last point are on the
    outer contour, so a later layer-to-layer solver can rotate/cut the cycle.
    """

    if not ordered_loops:
        return [], 1

    slots: list[tuple[int, int, int, int]] = []
    for loop in ordered_loops:
        a = nearest_index(loop.points, start_anchor)
        b = nearest_index(loop.points, exit_anchor)
        if abs(a - b) < max(3, len(loop.points) // 20):
            b = (a + len(loop.points) // 2) % len(loop.points)
        gap = seam_gap_points_for_width(loop, spacing, 0.65)
        slots.append(((a + gap) % len(loop.points), (b - gap) % len(loop.points), (a - gap) % len(loop.points), (b + gap) % len(loop.points)))

    path: list[Point] = []
    unsafe_connectors = 0
    required_clearance = line_width * 0.5 - spacing * 0.15
    direction = 1
    current_slot = 0

    def connect_to(point: Point) -> None:
        nonlocal unsafe_connectors
        if not path:
            append_point(path, point)
            return
        if dist(path[-1], point) <= EPS:
            return
        if segment_is_inside(model, path[-1], point, required_clearance, spacing * 0.5):
            append_point(path, point)
            return
        connector = grid_astar_connector(grid, path[-1], point, required_clearance=required_clearance)
        if connector is None:
            unsafe_connectors += 1
        else:
            append_points(path, connector)

    for i, loop in enumerate(ordered_loops):
        inward_slots = (slots[i][0], slots[i][1])
        start_idx = inward_slots[current_slot]
        end_idx = inward_slots[1 - current_slot]
        connect_to(loop.points[start_idx])
        append_points(path, closed_arc(loop.points, start_idx, end_idx, direction))
        current_slot = 1 - current_slot

    for i in reversed(range(len(ordered_loops))):
        loop = ordered_loops[i]
        outward_slots = (slots[i][2], slots[i][3])
        start_idx = outward_slots[current_slot]
        end_idx = outward_slots[1 - current_slot]
        connect_to(loop.points[start_idx])
        append_points(path, closed_arc(loop.points, start_idx, end_idx, direction))
        current_slot = 1 - current_slot

    return path, unsafe_connectors


def ordered_loops_for_spiral(contours: Sequence[ContourLoop]) -> list[ContourLoop]:
    grouped = loops_by_level(contours)
    if not grouped:
        return []

    first_count = len(grouped[min(grouped)])
    if first_count == 1:
        ordered: list[ContourLoop] = []
        for level in sorted(grouped):
            loops = grouped[level]
            if len(loops) == 1:
                ordered.append(loops[0])
            else:
                ordered.extend(sorted(loops, key=lambda loop: (loop.centroid[0], loop.centroid[1])))
        return ordered

    if first_count == 2 and all(len(grouped[level]) == 2 for level in grouped):
        outer_family: list[ContourLoop] = []
        hole_family: list[ContourLoop] = []
        for level in sorted(grouped):
            larger, smaller = sorted(grouped[level], key=lambda loop: loop.area, reverse=True)
            outer_family.append(larger)
            hole_family.append(smaller)
        return outer_family + list(reversed(hole_family))

    ordered = []
    for level in sorted(grouped):
        ordered.extend(sorted(grouped[level], key=lambda loop: -loop.area))
    return ordered


def subset_cycle_param(left: float, right: float, query: float, close_left: bool = False, close_right: bool = False) -> bool:
    if abs(query - left) <= 1e-9:
        return close_left
    if abs(query - right) <= 1e-9:
        return close_right
    if abs(left - right) <= 1e-9:
        return False
    if left < right:
        return left < query < right
    return query > left or query < right


def loop_point_at_index(points: Sequence[Point], index: float) -> Point:
    n = len(points)
    if n == 0:
        return (0.0, 0.0)
    index %= n
    base = math.floor(index)
    alpha = index - base
    return lerp(points[base % n], points[(base + 1) % n], alpha)


def loop_segment_length(points: Sequence[Point], index: int) -> float:
    n = len(points)
    return dist(points[index % n], points[(index + 1) % n])


def loop_nearest_index(points: Sequence[Point], p: Point) -> float:
    if not points:
        return 0.0

    best_idx = 0.0
    best_dist = float("inf")
    n = len(points)
    for i, a in enumerate(points):
        b = points[(i + 1) % n]
        ab = sub(b, a)
        denom = dot(ab, ab)
        t = 0.0 if denom <= EPS else max(0.0, min(1.0, dot(sub(p, a), ab) / denom))
        q = lerp(a, b, t)
        d = dist2(p, q)
        if d < best_dist:
            best_dist = d
            best_idx = i + t
    return best_idx % n


def loop_furthest_vertex_index(points: Sequence[Point], p: Point) -> float:
    best_idx = 0
    best_dist = -1.0
    for i, q in enumerate(points):
        d = dist2(p, q)
        if d > best_dist:
            best_dist = d
            best_idx = i
    return float(best_idx)


def loop_forward_by_distance(points: Sequence[Point], index: float, distance_along: float) -> float:
    if distance_along < 0.0:
        return loop_back_by_distance(points, index, -distance_along)

    n = len(points)
    index %= n
    p = loop_point_at_index(points, index)
    ceil_idx = math.ceil(index) % n
    walked = dist(p, points[ceil_idx])
    if walked >= distance_along and walked > EPS:
        alpha = distance_along / walked
        return (math.ceil(index) * alpha + index * (1.0 - alpha)) % n

    i = ceil_idx
    for _ in range(n + 2):
        segment = loop_segment_length(points, i)
        if walked + segment <= distance_along:
            walked += segment
            i = (i + 1) % n
        else:
            alpha = (distance_along - walked) / max(segment, EPS)
            return (i + alpha) % n
    return index


def loop_back_by_distance(points: Sequence[Point], index: float, distance_along: float) -> float:
    if distance_along < 0.0:
        return loop_forward_by_distance(points, index, -distance_along)

    n = len(points)
    index %= n
    p = loop_point_at_index(points, index)
    floor_idx = math.floor(index)
    walked = dist(p, points[floor_idx])
    if walked >= distance_along and walked > EPS:
        alpha = distance_along / walked
        return (floor_idx * alpha + index * (1.0 - alpha)) % n

    i = (floor_idx + n - 1) % n
    for _ in range(n + 2):
        segment = loop_segment_length(points, i)
        if walked + segment <= distance_along:
            walked += segment
            i = (i + n - 1) % n
        else:
            alpha = (distance_along - walked) / max(segment, EPS)
            return (i + 1.0 - alpha) % n
    return index


def loop_length_between(points: Sequence[Point], start: float, end: float) -> float:
    n = len(points)
    start %= n
    end %= n
    start_floor = math.floor(start)
    end_floor = math.floor(end)
    if start_floor == end_floor and end > start:
        return (end - start) * loop_segment_length(points, start_floor)

    length = (1.0 - start + start_floor) * loop_segment_length(points, start_floor)
    length += (end - end_floor) * loop_segment_length(points, end_floor)
    i = (start_floor + 1) % n
    while i != end_floor:
        length += loop_segment_length(points, i)
        i = (i + 1) % n
    return length


def append_loop_param(path: list[Point], points: Sequence[Point], index: float) -> None:
    append_point(path, loop_point_at_index(points, index))


def append_loop_vertices_forward(
    path: list[Point],
    points: Sequence[Point],
    start_index: int,
    stop_index: int,
    guard_extra: int = 2,
) -> None:
    n = len(points)
    i = start_index % n
    guard = 0
    while i != stop_index % n and guard <= n + guard_extra:
        append_point(path, points[i])
        i = (i + 1) % n
        guard += 1


def append_loop_vertices_backward_until(
    path: list[Point],
    points: Sequence[Point],
    start_index: int,
    left: float,
    right: float,
) -> None:
    n = len(points)
    i = start_index % n
    guard = 0
    while subset_cycle_param(left, float(i), right, False, False) and guard <= n + 2:
        append_point(path, points[i])
        i = (i + n - 1) % n
        guard += 1


def build_single_minimum_connected_fermat(
    ordered_loops: Sequence[ContourLoop],
    start_anchor: Point,
    spacing: float,
    port_spacing: float | None = None,
    exit_anchor: Point | None = None,
    preserve_medial_pockets: bool = True,
) -> tuple[list[Point], int]:
    """Build a boundary-to-boundary connected Fermat spiral for one contour chain.

    This is a direct contour-index version of the CFS single-minimum connector:
    each contour owns an in point and an out point. The in branch walks inward,
    the out branch walks inward in parallel, and the final path appends the out
    branch in reverse. The result starts and ends on the outer contour without
    the slot self-touches produced by the earlier two-slot weave.
    """

    loops = [loop for loop in ordered_loops if len(loop.points) >= 4]
    if not loops:
        return [], 0
    port_spacing = spacing if port_spacing is None else port_spacing
    if len(loops) == 1:
        start = loop_nearest_index(loops[0].points, start_anchor)
        return full_loop(loops[0].points, int(round(start)), 1), 0

    # Keep the small medial contours for concave pocket shapes. One-hole rings
    # are different: opposing fronts can become closer than one line spacing,
    # so their final medial caps are not printable without spacing violations.
    if not preserve_medial_pockets:
        while len(loops) > 2 and loops[-1].length < spacing * 20.0:
            loops.pop()

    in_index = loop_nearest_index(loops[0].points, start_anchor)
    out_index = loop_nearest_index(loops[0].points, exit_anchor) if exit_anchor is not None else in_index
    out_forward_in = True
    in_run = False
    first_circle = True
    in_branch: list[Point] = []
    out_branch: list[Point] = []
    loop_index = 0

    while True:
        loop = loops[loop_index]
        points = loop.points
        n = len(points)
        circle_small = loop.length < port_spacing * 2.0

        if circle_small:
            near_index = (0.5 * (in_index + out_index)) % n
            near = loop_point_at_index(points, near_index)
            far_index = loop_furthest_vertex_index(points, near)
            far = loop_point_at_index(points, far_index)

            in_index = float(math.ceil(near_index) % n)

            def small_score(idx: float) -> float:
                p = points[int(idx) % n]
                return min(dist(near, p), dist(far, p))

            best = small_score(in_index)
            i = (int(in_index) + 1) % n
            while subset_cycle_param(near_index, far_index, float(i)):
                score = small_score(float(i))
                if score > best:
                    best = score
                    in_index = float(i)
                i = (i + 1) % n

            out_index = float(math.ceil(far_index) % n)
            best = small_score(out_index)
            i = (int(out_index) + 1) % n
            while subset_cycle_param(far_index, near_index, float(i)):
                score = small_score(float(i))
                if score > best:
                    best = score
                    out_index = float(i)
                i = (i + 1) % n

            out_forward_in = loop_length_between(points, in_index, out_index) <= loop_length_between(points, out_index, in_index)
        else:
            if abs(in_index - out_index) < 1e-6:
                if first_circle:
                    out_forward = loop_forward_by_distance(points, in_index, port_spacing)
                    out_backward = loop_back_by_distance(points, in_index, port_spacing)
                    next_loop = loops[min(loop_index + 1, len(loops) - 1)]
                    if distance_to_loop(loop_point_at_index(points, out_forward), next_loop.points) < distance_to_loop(
                        loop_point_at_index(points, out_backward), next_loop.points
                    ):
                        out_index = out_backward
                        out_forward_in = False
                    else:
                        out_index = out_forward
                        out_forward_in = True
                else:
                    in0 = loop_point_at_index(points, in_index)
                    out_forward = loop_forward_by_distance(points, in_index, port_spacing)
                    out_backward = loop_back_by_distance(points, in_index, port_spacing)
                    p_forward = loop_point_at_index(points, out_forward)
                    p_backward = loop_point_at_index(points, out_backward)
                    prev_in = in_branch[-1]
                    prev_out = out_branch[-1]

                    forward_score = min(dist(in0, prev_in) + dist(p_forward, prev_out), dist(p_forward, prev_in) + dist(in0, prev_out))
                    backward_score = min(
                        dist(in0, prev_in) + dist(p_backward, prev_out), dist(p_backward, prev_in) + dist(in0, prev_out)
                    )
                    if forward_score < backward_score:
                        if dist(in0, prev_in) + dist(p_forward, prev_out) <= dist(p_forward, prev_in) + dist(in0, prev_out):
                            out_index = out_forward
                            out_forward_in = True
                        else:
                            out_index = in_index
                            in_index = out_forward
                            out_forward_in = False
                    else:
                        if dist(in0, prev_in) + dist(p_backward, prev_out) <= dist(p_backward, prev_in) + dist(in0, prev_out):
                            out_index = out_backward
                            out_forward_in = False
                        else:
                            out_index = in_index
                            in_index = out_backward
                            out_forward_in = True
            else:
                length_io = loop_length_between(points, in_index, out_index)
                length_oi = loop_length_between(points, out_index, in_index)
                out_forward_in = length_io <= length_oi
                arc_length = min(length_io, length_oi)
                if first_circle:
                    if arc_length < port_spacing:
                        if out_forward_in:
                            out_index = loop_forward_by_distance(points, in_index, port_spacing)
                        else:
                            out_index = loop_back_by_distance(points, in_index, port_spacing)
                    out_forward_in = loop_length_between(points, in_index, out_index) <= loop_length_between(
                        points, out_index, in_index
                    )
                    append_loop_param(in_branch, points, in_index)
                    append_loop_param(out_branch, points, out_index)
                    if loop_index + 1 < len(loops):
                        if out_forward_in:
                            in_back = loop_back_by_distance(points, in_index, port_spacing)
                            i = math.ceil(out_index) % n
                            guard = 0
                            while subset_cycle_param(float(i), in_index, in_back, False, False) and guard <= n + 2:
                                append_point(out_branch, points[i])
                                i = (i + 1) % n
                                guard += 1
                            append_loop_param(out_branch, points, in_back)
                        else:
                            in_far = loop_forward_by_distance(points, in_index, port_spacing)
                            i = math.floor(out_index) % n
                            append_loop_vertices_backward_until(out_branch, points, i, in_index, in_far)
                            append_loop_param(out_branch, points, in_far)

                        in_run = True
                        loop_index += 1
                        child_points = loops[loop_index].points
                        in_index = loop_nearest_index(child_points, in_branch[-1])
                        out_index = loop_nearest_index(child_points, out_branch[-1])
                        first_circle = False
                        continue
                in0 = loop_point_at_index(points, in_index)
                out0 = loop_point_at_index(points, out_index)

                if out_forward_in:
                    in1_index = loop_forward_by_distance(points, in_index, arc_length - port_spacing)
                    out1_index = loop_back_by_distance(points, out_index, arc_length - port_spacing)
                else:
                    in1_index = loop_back_by_distance(points, in_index, arc_length - port_spacing)
                    out1_index = loop_forward_by_distance(points, out_index, arc_length - port_spacing)

                in1 = loop_point_at_index(points, in1_index)
                out1 = loop_point_at_index(points, out1_index)
                prev_in = in_branch[-1]
                prev_out = out_branch[-1]

                keep_in_score = min(dist(in0, prev_in) + dist(out1, prev_out), dist(in0, prev_out) + dist(out1, prev_in))
                keep_out_score = min(dist(in1, prev_in) + dist(out0, prev_out), dist(in1, prev_out) + dist(out0, prev_in))
                if keep_in_score <= keep_out_score:
                    if dist(in0, prev_in) + dist(out1, prev_out) < dist(in0, prev_out) + dist(out1, prev_in):
                        out_index = out1_index
                    else:
                        out_index = in_index
                        in_index = out1_index
                        out_forward_in = not out_forward_in
                else:
                    if dist(in1, prev_in) + dist(out0, prev_out) < dist(in1, prev_out) + dist(out0, prev_in):
                        in_index = in1_index
                    else:
                        in_index = out_index
                        out_index = in1_index
                        out_forward_in = not out_forward_in

        if not first_circle and segments_intersect(
            loop_point_at_index(points, in_index),
            in_branch[-1],
            loop_point_at_index(points, out_index),
            out_branch[-1],
        ):
            in_index, out_index = out_index, in_index
            out_forward_in = not out_forward_in

        append_loop_param(in_branch, points, in_index)
        append_loop_param(out_branch, points, out_index)

        if loop_index + 1 < len(loops):
            if in_run:
                if out_forward_in:
                    out_far = loop_forward_by_distance(points, out_index, port_spacing)
                    i = math.floor(in_index) % n
                    append_loop_vertices_backward_until(in_branch, points, i, out_index, out_far)
                    append_loop_param(in_branch, points, out_far)
                else:
                    out_back = loop_back_by_distance(points, out_index, port_spacing)
                    i = math.ceil(in_index) % n
                    guard = 0
                    while subset_cycle_param(float(i), out_index, out_back, False, False) and guard <= n + 2:
                        append_point(in_branch, points[i])
                        i = (i + 1) % n
                        guard += 1
                    append_loop_param(in_branch, points, out_back)
            else:
                if out_forward_in:
                    in_back = loop_back_by_distance(points, in_index, port_spacing)
                    i = math.ceil(out_index) % n
                    guard = 0
                    while subset_cycle_param(float(i), in_index, in_back, False, False) and guard <= n + 2:
                        append_point(out_branch, points[i])
                        i = (i + 1) % n
                        guard += 1
                    append_loop_param(out_branch, points, in_back)
                else:
                    in_far = loop_forward_by_distance(points, in_index, port_spacing)
                    i = math.floor(out_index) % n
                    append_loop_vertices_backward_until(out_branch, points, i, in_index, in_far)
                    append_loop_param(out_branch, points, in_far)

            in_run = not in_run
            loop_index += 1
            child_points = loops[loop_index].points
            in_index = loop_nearest_index(child_points, in_branch[-1])
            out_index = loop_nearest_index(child_points, out_branch[-1])
        else:
            prev_in = in_branch[-1]
            prev_out = out_branch[-1]
            medial_index = max(range(n), key=lambda i: min(dist(points[i], prev_in), dist(points[i], prev_out)))
            if subset_cycle_param(in_index, out_index, float(medial_index), True, False):
                append_loop_vertices_forward(in_branch, points, math.ceil(in_index) % n, math.floor(out_index) % n)
            else:
                append_loop_vertices_forward(out_branch, points, math.ceil(out_index) % n, math.floor(in_index) % n)
            break

        first_circle = False

    for p in reversed(out_branch):
        append_point(in_branch, p)
    cleaned: list[Point] = []
    append_points(cleaned, in_branch)
    return cleaned, 0


def project_open_polyline(points: Sequence[Point], p: Point) -> tuple[float, Point, float]:
    if len(points) < 2:
        q = points[0] if points else (0.0, 0.0)
        return 0.0, q, dist(p, q)

    best_index = 0.0
    best_point = points[0]
    best_distance2 = float("inf")
    for i, (a, b) in enumerate(zip(points, points[1:])):
        ab = sub(b, a)
        denom = dot(ab, ab)
        t = 0.0 if denom <= EPS else max(0.0, min(1.0, dot(sub(p, a), ab) / denom))
        q = lerp(a, b, t)
        d = dist2(p, q)
        if d < best_distance2:
            best_distance2 = d
            best_index = i + t
            best_point = q
    return best_index, best_point, math.sqrt(best_distance2)


def cut_open_polyline(points: Sequence[Point], index: float) -> tuple[list[Point], list[Point], Point]:
    if not points:
        return [], [], (0.0, 0.0)
    if len(points) == 1:
        return [points[0]], [points[0]], points[0]

    index = max(0.0, min(index, len(points) - 1.0))
    base = min(math.floor(index), len(points) - 2)
    alpha = index - base
    cut = lerp(points[base], points[base + 1], alpha)
    before = list(points[: base + 1])
    append_point(before, cut)
    after = [cut]
    after.extend(points[base + 1 :])
    return before, after, cut


def merge_child_spiral(parent_path: Sequence[Point], child_path: Sequence[Point], spacing: float) -> list[Point]:
    if not parent_path:
        return list(child_path)
    if not child_path:
        return list(parent_path)

    start_index, _, _ = project_open_polyline(parent_path, child_path[0])
    end_index, _, _ = project_open_polyline(parent_path, child_path[-1])
    child = list(child_path)
    if start_index > end_index:
        child.reverse()
        start_index, _, _ = project_open_polyline(parent_path, child[0])
        end_index, _, _ = project_open_polyline(parent_path, child[-1])

    parent_start, _, _ = cut_open_polyline(parent_path, start_index)
    _, parent_end, _ = cut_open_polyline(parent_path, end_index)
    merged: list[Point] = []
    append_points(merged, parent_start)
    append_points(merged, child)
    append_points(merged, parent_end)
    return simplify_polyline(merged, tolerance=spacing * 0.04)


def complete_outer_boundary_cycle(
    path: Sequence[Point],
    outer_loop: ContourLoop,
    spacing: float,
    spacing_tolerance: float,
) -> list[Point]:
    """Append the missing outer-contour arc so the CFS can be cut anywhere.

    The CFS routines produce a boundary-to-boundary curve. For layer-to-layer
    continuity we need the layer path to behave as a closed cycle with only a
    tiny boundary cut gap. Try both possible outer-boundary arcs and keep the
    one that preserves the non-crossing/equal-spacing validator best.
    """

    if len(path) < 2 or len(outer_loop.points) < 4:
        return list(path)

    target_gap = spacing * 0.12
    if dist(path[0], path[-1]) <= target_gap * 1.5:
        return list(path)

    start_index = loop_nearest_index(outer_loop.points, path[0])
    end_index = loop_nearest_index(outer_loop.points, path[-1])
    candidates: list[tuple[tuple[int, int, float, float], list[Point]]] = []

    for direction in (1, -1):
        if direction > 0:
            stop_index = loop_back_by_distance(outer_loop.points, start_index, target_gap)
        else:
            stop_index = loop_forward_by_distance(outer_loop.points, start_index, target_gap)

        arc = loop_arc_between_params(outer_loop.points, end_index, stop_index, direction)
        if len(arc) < 2:
            continue
        arc[0] = path[-1]

        candidate: list[Point] = []
        append_points(candidate, path)
        append_points(candidate, arc[1:])
        if dist(candidate[0], candidate[-1]) > spacing * 0.20:
            continue

        crossings, close_pairs, min_spacing = path_pair_metrics(
            candidate,
            spacing=spacing,
            spacing_tolerance=spacing_tolerance,
        )
        arc_length = polyline_length(arc, closed=False)
        candidates.append(((crossings, close_pairs, -min_spacing, arc_length), candidate))

    if not candidates:
        return list(path)

    candidates.sort(key=lambda item: item[0])
    return candidates[0][1]


def open_path_prefix_lengths(path: Sequence[Point]) -> list[float]:
    prefix = [0.0]
    for a, b in zip(path, path[1:]):
        prefix.append(prefix[-1] + dist(a, b))
    return prefix


def open_path_point_at_index(path: Sequence[Point], index: float) -> Point:
    if not path:
        return (0.0, 0.0)
    if len(path) == 1:
        return path[0]
    index = max(0.0, min(index, len(path) - 1.0))
    base = min(math.floor(index), len(path) - 2)
    return lerp(path[base], path[base + 1], index - base)


def open_path_length_at_index(path: Sequence[Point], prefix: Sequence[float], index: float) -> float:
    if len(path) < 2:
        return 0.0
    index = max(0.0, min(index, len(path) - 1.0))
    base = min(math.floor(index), len(path) - 2)
    return prefix[base] + (index - base) * dist(path[base], path[base + 1])


def open_path_index_at_length(path: Sequence[Point], prefix: Sequence[float], target: float) -> float:
    if len(path) < 2:
        return 0.0

    target = max(0.0, min(target, prefix[-1]))
    lo = 0
    hi = len(prefix) - 1
    while lo < hi:
        mid = (lo + hi) // 2
        if prefix[mid] < target:
            lo = mid + 1
        else:
            hi = mid

    idx = max(0, lo - 1)
    segment_length = max(dist(path[idx], path[idx + 1]), EPS)
    return idx + (target - prefix[idx]) / segment_length


def count_containment_violations(model: PolygonModel, path: Sequence[Point], spacing: float) -> int:
    violations = 0
    outside_tolerance = -spacing * 0.05
    for a, b in zip(path, path[1:]):
        length = dist(a, b)
        samples = max(2, int(math.ceil(length / max(spacing * 0.25, EPS))))
        for i in range(samples + 1):
            if model.sdf(lerp(a, b, i / samples)) < outside_tolerance:
                violations += 1
    return violations


def loop_uncovered_fraction(loop: ContourLoop, path: Sequence[Point], spacing: float) -> float:
    if not path:
        return 1.0

    stride = max(1, len(loop.points) // 32)
    sampled = loop.points[::stride]
    if not sampled:
        return 0.0

    uncovered = 0
    for p in sampled:
        _, _, distance = project_open_polyline(path, p)
        if distance > spacing * 0.55:
            uncovered += 1
    return uncovered / len(sampled)


def uncovered_pocket_groups(contours: Sequence[ContourLoop], path: Sequence[Point], spacing: float) -> list[list[ContourLoop]]:
    uncovered = [
        loop
        for loop in contours
        if loop.level_index > 0 and loop.area < spacing * spacing * 90.0 and loop_uncovered_fraction(loop, path, spacing) > 0.65
    ]

    groups: list[list[ContourLoop]] = []
    for loop in sorted(uncovered, key=lambda item: -item.area):
        for group in groups:
            if dist(loop.centroid, group[0].centroid) < spacing * 7.0:
                group.append(loop)
                break
        else:
            groups.append([loop])
    return groups


def unique_port_candidates(candidates: Sequence[tuple[float, float]], limit: int) -> list[float]:
    selected: list[float] = []
    for _, index in sorted(candidates, key=lambda item: item[0]):
        if all(abs(index - existing) > 2.0 for existing in selected):
            selected.append(index)
            if len(selected) >= limit:
                break
    return selected


def try_insert_pocket_group(
    model: PolygonModel,
    parent_path: Sequence[Point],
    group: Sequence[ContourLoop],
    spacing: float,
    spacing_tolerance: float,
) -> list[Point] | None:
    if len(parent_path) < 4 or not group:
        return None

    prefix = open_path_prefix_lengths(parent_path)
    centroid = (
        sum(loop.centroid[0] for loop in group) / len(group),
        sum(loop.centroid[1] for loop in group) / len(group),
    )
    port_candidates: list[tuple[float, float]] = []
    index, _, distance = project_open_polyline(parent_path, centroid)
    port_candidates.append((distance, index))
    for loop in group:
        stride = max(1, len(loop.points) // 16)
        for p in loop.points[::stride]:
            index, _, distance = project_open_polyline(parent_path, p)
            port_candidates.append((distance, index))

    ordered = ordered_loops_for_spiral(group)
    best: tuple[tuple[int, int, int, float], list[Point]] | None = None
    for center_index in unique_port_candidates(port_candidates, limit=2):
        center_index = max(3.0, min(center_index, len(parent_path) - 4.0))
        center_length = open_path_length_at_index(parent_path, prefix, center_index)

        for half_width in (spacing * 1.35,):
            left_index = open_path_index_at_length(parent_path, prefix, center_length - half_width)
            right_index = open_path_index_at_length(parent_path, prefix, center_length + half_width)

            left_port = open_path_point_at_index(parent_path, left_index)
            right_port = open_path_point_at_index(parent_path, right_index)
            if dist(left_port, right_port) < spacing * 0.5:
                continue
            for start_anchor, exit_anchor in ((left_port, right_port), (right_port, left_port)):
                child_path, _ = build_single_minimum_connected_fermat(
                    ordered,
                    start_anchor=start_anchor,
                    spacing=spacing,
                    port_spacing=spacing * 2.5,
                    exit_anchor=exit_anchor,
                    preserve_medial_pockets=True,
                )
                for candidate_child in (child_path, list(reversed(child_path))):
                    before, _, _ = cut_open_polyline(parent_path, left_index)
                    _, after, _ = cut_open_polyline(parent_path, right_index)
                    candidate: list[Point] = []
                    append_points(candidate, before)
                    append_points(candidate, candidate_child)
                    append_points(candidate, after)

                    crossings, close_pairs, min_spacing = path_pair_metrics(
                        candidate,
                        spacing=spacing,
                        spacing_tolerance=spacing_tolerance,
                    )
                    containment = 0
                    if crossings == 0 and close_pairs == 0:
                        containment = count_containment_violations(model, candidate, spacing)

                    score = (crossings, close_pairs, containment, -min_spacing)
                    if best is None or score < best[0]:
                        best = (score, candidate)
                    if crossings == 0 and close_pairs == 0 and containment == 0:
                        return candidate

    if best is not None and best[0][0] == 0 and best[0][1] == 0 and best[0][2] == 0:
        return best[1]
    return None


def insert_uncovered_pocket_spirals(
    model: PolygonModel,
    contours: Sequence[ContourLoop],
    path: Sequence[Point],
    spacing: float,
    spacing_tolerance: float,
) -> tuple[list[Point], int]:
    current = list(path)
    inserted = 0
    for group in uncovered_pocket_groups(contours, current, spacing):
        candidate = try_insert_pocket_group(
            model,
            current,
            group,
            spacing=spacing,
            spacing_tolerance=spacing_tolerance,
        )
        if candidate is not None:
            current = candidate
            inserted += 1
    return current, inserted


def build_branch_connected_fermat(
    contours: Sequence[ContourLoop],
    start_anchor: Point,
    spacing: float,
) -> tuple[list[Point], int]:
    grouped = loops_by_level(contours)
    levels = sorted(grouped)
    split_level = next((level for level in levels if len(grouped[level]) > 1), None)
    if split_level is None:
        return build_single_minimum_connected_fermat(
            ordered_loops_for_spiral(contours),
            start_anchor=start_anchor,
            spacing=spacing,
            port_spacing=spacing * 2.5,
        )

    trunk = [grouped[level][0] for level in levels if level < split_level and len(grouped[level]) == 1]
    if not trunk:
        return build_single_minimum_connected_fermat(
            ordered_loops_for_spiral(contours),
            start_anchor=start_anchor,
            spacing=spacing,
            port_spacing=spacing * 2.5,
        )

    parent_path, unsafe = build_single_minimum_connected_fermat(
        trunk,
        start_anchor=start_anchor,
        spacing=spacing,
        port_spacing=spacing * 2.5,
    )

    chains = [[loop] for loop in sorted(grouped[split_level], key=lambda loop: (loop.centroid[0], loop.centroid[1]))]
    for level in (level for level in levels if level > split_level):
        remaining = list(grouped[level])
        used: set[int] = set()
        for chain in chains:
            last = chain[-1]
            best: tuple[float, int, ContourLoop] | None = None
            for idx, loop in enumerate(remaining):
                if idx in used:
                    continue
                d = dist2(last.centroid, loop.centroid)
                if best is None or d < best[0]:
                    best = (d, idx, loop)
            if best is not None:
                used.add(best[1])
                chain.append(best[2])

    for chain in chains:
        stride = max(1, len(chain[0].points) // 48)
        anchor = min((p for p in chain[0].points[::stride]), key=lambda p: project_open_polyline(parent_path, p)[2])
        child_path, child_unsafe = build_single_minimum_connected_fermat(
            chain,
            start_anchor=anchor,
            spacing=spacing,
            port_spacing=spacing * 2.5,
        )
        unsafe += child_unsafe
        parent_path = merge_child_spiral(parent_path, child_path, spacing)

    return parent_path, unsafe


def build_horizontal_capsule_loop(
    left_center: Point,
    right_center: Point,
    radius: float,
    level_index: int,
    offset: float,
    segments: int = 64,
) -> ContourLoop:
    cy = (left_center[1] + right_center[1]) * 0.5
    left_x = left_center[0]
    right_x = right_center[0]
    half_segments = max(8, segments // 2)
    points: list[Point] = [(left_x, cy + radius), (right_x, cy + radius)]

    for i in range(1, half_segments + 1):
        angle = math.pi * 0.5 - math.pi * i / half_segments
        points.append((right_x + math.cos(angle) * radius, cy + math.sin(angle) * radius))

    points.append((left_x, cy - radius))

    for i in range(1, half_segments + 1):
        angle = -math.pi * 0.5 - math.pi * i / half_segments
        points.append((left_x + math.cos(angle) * radius, cy + math.sin(angle) * radius))

    if signed_area(points) < 0.0:
        points.reverse()
    return ContourLoop(
        level_index=level_index,
        offset=offset,
        points=points,
        area=abs(signed_area(points)),
        length=polyline_length(points, closed=True),
        centroid=polygon_centroid(points),
    )


def build_merged_hole_connected_fermat(
    model: PolygonModel,
    contours: Sequence[ContourLoop],
    start_anchor: Point,
    line_width: float,
    spacing: float,
) -> tuple[list[Point], int]:
    """Conservative multi-hole fallback.

    A strict no-travel CFS cannot always preserve independent hole islands
    without a multi-slot contour tree. For the prototype, merge multiple holes
    into one synthetic capsule cavity and fill the printable area around that
    cavity. This intentionally leaves the bridge material between holes
    unprinted, which is acceptable for this continuous mode when topology would
    otherwise force crossings.
    """

    grouped = loops_by_level(contours)
    outer: list[ContourLoop] = []
    min_outer_area = spacing * spacing * 20.0
    for level in sorted(grouped):
        loops = sorted(grouped[level], key=lambda loop: loop.area, reverse=True)
        if loops and loops[0].area > min_outer_area:
            outer.append(loops[0])

    if not outer or len(model.holes) < 2:
        return build_single_minimum_connected_fermat(
            ordered_loops_for_spiral(contours),
            start_anchor=start_anchor,
            spacing=spacing,
            port_spacing=spacing * 2.5,
        )

    # Keep only the large outer-family contours; later medial fragments are
    # children of the hole topology and are intentionally replaced by the
    # conservative capsule cavity.
    stable_outer: list[ContourLoop] = []
    previous_area = None
    for loop in outer:
        if previous_area is not None and loop.area < previous_area * 0.2:
            break
        stable_outer.append(loop)
        previous_area = loop.area
    outer = stable_outer

    centers = sorted((polygon_centroid(hole) for hole in model.holes), key=lambda p: p[0])
    left_center = centers[0]
    right_center = centers[-1]
    base_radius = 0.0
    for hole, center in zip(model.holes, centers):
        base_radius = max(base_radius, max(dist(p, center) for p in hole))
    base_radius += line_width * 0.5

    _, min_y, _, max_y = model.bounds(margin=0.0)
    center_y = (left_center[1] + right_center[1]) * 0.5
    max_radius = min(center_y - min_y, max_y - center_y) - line_width

    inner: list[ContourLoop] = []
    for level_index in range(len(outer)):
        radius = base_radius + level_index * spacing
        if radius > max_radius:
            break
        inner.append(
            build_horizontal_capsule_loop(
                left_center,
                right_center,
                radius,
                level_index=level_index,
                offset=line_width * 0.5 + level_index * spacing,
            )
        )

    if not inner:
        return build_single_minimum_connected_fermat(
            ordered_loops_for_spiral(contours),
            start_anchor=start_anchor,
            spacing=spacing,
            port_spacing=spacing * 2.5,
        )

    ordered = outer[: len(inner)] + list(reversed(inner))
    candidate_anchors = [start_anchor]
    for fraction in (0.04, 0.06, 0.10, 0.25, 0.50, 0.75):
        candidate_anchors.append(point_at_closed_fraction(ordered[0].points, fraction))

    best: tuple[tuple[int, int, float], list[Point], int] | None = None
    for anchor in candidate_anchors:
        candidate_path, unsafe = build_single_minimum_connected_fermat(
            ordered,
            start_anchor=anchor,
            spacing=spacing,
            port_spacing=spacing * 2.5,
        )
        crossings, close_pairs, min_spacing = path_pair_metrics(
            candidate_path,
            spacing=spacing,
            spacing_tolerance=0.25,
        )
        score = (crossings, close_pairs, -min_spacing)
        if best is None or score < best[0]:
            best = (score, candidate_path, unsafe)
        if crossings == 0 and close_pairs == 0:
            return candidate_path, unsafe

    assert best is not None
    return best[1], best[2]


def close_path_to_cycle(
    model: PolygonModel,
    grid: SDFGrid,
    path: list[Point],
    line_width: float,
    spacing: float,
) -> int:
    if len(path) < 2:
        return 1
    if dist(path[-1], path[0]) <= EPS:
        return 0

    required_clearance = line_width * 0.5 - spacing * 0.15
    if segment_is_inside(model, path[-1], path[0], required_clearance, spacing * 0.5):
        append_point(path, path[0])
        return 0

    connector = grid_astar_connector(grid, path[-1], path[0], required_clearance=required_clearance)
    if connector is None:
        return 1
    append_points(path, connector)
    if dist(path[-1], path[0]) > EPS:
        append_point(path, path[0])
    return 0


def segment_is_inside(model: PolygonModel, a: Point, b: Point, required_clearance: float, step: float) -> bool:
    length = dist(a, b)
    samples = max(2, int(math.ceil(length / max(step, EPS))))
    for i in range(samples + 1):
        p = lerp(a, b, i / samples)
        if model.sdf(p) < required_clearance:
            return False
    return True


def sampled_connection_to_loop(
    model: PolygonModel,
    start: Point,
    loop: ContourLoop,
    required_clearance: float,
    spacing: float,
) -> tuple[int, float, bool]:
    stride = max(1, len(loop.points) // 96)
    candidates: list[tuple[float, int]] = []
    for idx in range(0, len(loop.points), stride):
        p = loop.points[idx]
        candidates.append((dist(start, p), idx))

    candidates.sort(key=lambda item: item[0])
    for candidate_dist, idx in candidates:
        if segment_is_inside(model, start, loop.points[idx], required_clearance, spacing * 0.7):
            return idx, candidate_dist, True
    best_dist, best_idx = candidates[0]
    return best_idx, best_dist, False


def nearest_valid_grid_node(grid: SDFGrid, p: Point, required_clearance: float) -> tuple[int, int] | None:
    ix = int(round((p[0] - grid.min_x) / max(grid.dx, EPS)))
    iy = int(round((p[1] - grid.min_y) / max(grid.dy, EPS)))
    ix = max(0, min(grid.nx, ix))
    iy = max(0, min(grid.ny, iy))

    best: tuple[float, int, int] | None = None
    max_radius = max(grid.nx, grid.ny)
    for radius in range(max_radius + 1):
        found_at_radius = False
        for y in range(max(0, iy - radius), min(grid.ny, iy + radius) + 1):
            for x in range(max(0, ix - radius), min(grid.nx, ix + radius) + 1):
                if abs(x - ix) != radius and abs(y - iy) != radius:
                    continue
                if grid.value(x, y) >= required_clearance:
                    d = dist2(p, grid.point(x, y))
                    if best is None or d < best[0]:
                        best = (d, x, y)
                        found_at_radius = True
        if found_at_radius:
            break
    if best is None:
        return None
    return (best[1], best[2])


def grid_astar_connector(grid: SDFGrid, start: Point, goal: Point, required_clearance: float) -> list[Point] | None:
    start_node = nearest_valid_grid_node(grid, start, required_clearance)
    goal_node = nearest_valid_grid_node(grid, goal, required_clearance)
    if start_node is None or goal_node is None:
        return None
    if start_node == goal_node:
        return [start, goal]

    def heuristic(node: tuple[int, int]) -> float:
        return dist(grid.point(node[0], node[1]), grid.point(goal_node[0], goal_node[1]))

    neighbors = [
        (-1, -1, math.sqrt(2.0)),
        (0, -1, 1.0),
        (1, -1, math.sqrt(2.0)),
        (-1, 0, 1.0),
        (1, 0, 1.0),
        (-1, 1, math.sqrt(2.0)),
        (0, 1, 1.0),
        (1, 1, math.sqrt(2.0)),
    ]

    open_heap: list[tuple[float, tuple[int, int]]] = [(heuristic(start_node), start_node)]
    came_from: dict[tuple[int, int], tuple[int, int]] = {}
    g_score: dict[tuple[int, int], float] = {start_node: 0.0}
    closed: set[tuple[int, int]] = set()

    while open_heap:
        _, current = heapq.heappop(open_heap)
        if current in closed:
            continue
        if current == goal_node:
            nodes = [current]
            while nodes[-1] in came_from:
                nodes.append(came_from[nodes[-1]])
            nodes.reverse()
            points = [start]
            points.extend(grid.point(x, y) for x, y in nodes)
            points.append(goal)
            return simplify_polyline(points, tolerance=grid.step * 0.35)

        closed.add(current)
        cx, cy = current
        for ox, oy, step_cost in neighbors:
            nx = cx + ox
            ny = cy + oy
            if nx < 0 or nx > grid.nx or ny < 0 or ny > grid.ny:
                continue
            if grid.value(nx, ny) < required_clearance:
                continue
            neighbor = (nx, ny)
            tentative = g_score[current] + step_cost * grid.step
            if tentative < g_score.get(neighbor, float("inf")):
                came_from[neighbor] = current
                g_score[neighbor] = tentative
                heapq.heappush(open_heap, (tentative + heuristic(neighbor), neighbor))

    return None


def connector_clear_of_path(
    existing_path: Sequence[Point],
    connector: Sequence[Point],
    spacing: float,
    spacing_tolerance: float,
    skip_tail: int = 10,
) -> bool:
    if len(existing_path) < 2 or len(connector) < 2:
        return True

    min_allowed = spacing * (1.0 - spacing_tolerance)
    existing_segments = list(zip(existing_path, existing_path[1:]))
    connector_segments = list(zip(connector, connector[1:]))
    compare_until = max(0, len(existing_segments) - skip_tail)

    for ca, cb in connector_segments:
        for ea, eb in existing_segments[:compare_until]:
            d = segment_distance(ca, cb, ea, eb)
            if d <= 1e-7 or d < min_allowed:
                return False
    return True


def choose_self_avoiding_entry(
    model: PolygonModel,
    path: Sequence[Point],
    current_point: Point,
    loop: ContourLoop,
    line_width: float,
    spacing: float,
    spacing_tolerance: float,
) -> tuple[int, float, bool]:
    required_clearance = line_width * 0.5 - spacing * 0.15
    stride = max(1, len(loop.points) // 80)
    candidates: list[tuple[float, int]] = [(dist(current_point, loop.points[idx]), idx) for idx in range(0, len(loop.points), stride)]
    candidates.sort(key=lambda item: item[0])

    for connector_length, idx in candidates:
        p = loop.points[idx]
        if not segment_is_inside(model, current_point, p, required_clearance, spacing * 0.7):
            continue
        if connector_clear_of_path(path, [current_point, p], spacing, spacing_tolerance):
            return idx, connector_length, True

    best_length, best_idx = candidates[0]
    return best_idx, best_length, False


def build_self_avoiding_open_spiral(
    model: PolygonModel,
    grid: SDFGrid,
    contours: Sequence[ContourLoop],
    anchor: Point,
    line_width: float,
    spacing: float,
    spacing_tolerance: float,
) -> tuple[list[Point], int]:
    if not contours:
        return [], 0

    remaining = set(range(len(contours)))
    start_idx = max(remaining, key=lambda idx: contours[idx].area)
    path: list[Point] = []
    unsafe_connectors = 0
    current_point = anchor
    direction = 1
    required_clearance = line_width * 0.5 - spacing * 0.15
    last_area = contours[start_idx].area

    while remaining:
        if not path:
            loop_idx = start_idx
            loop = contours[loop_idx]
            entry_idx = nearest_index(loop.points, current_point)
            safe = True
        else:
            best: tuple[float, int, int, bool] | None = None
            candidate_order = sorted(
                remaining,
                key=lambda idx: (
                    dist2(current_point, contours[idx].centroid),
                    abs(math.sqrt(max(contours[idx].area, 0.0)) - math.sqrt(max(last_area, 0.0))),
                ),
            )[:10]
            for candidate_idx in candidate_order:
                loop = contours[candidate_idx]
                entry_idx, connector_length, safe = choose_self_avoiding_entry(
                    model,
                    path,
                    current_point,
                    loop,
                    line_width=line_width,
                    spacing=spacing,
                    spacing_tolerance=spacing_tolerance,
                )
                if not safe:
                    continue
                area_penalty = abs(math.sqrt(max(loop.area, 0.0)) - math.sqrt(max(last_area, 0.0))) * 0.35
                level_penalty = abs(loop.level_index - contours[start_idx].level_index) * spacing * 0.02
                score = connector_length + area_penalty + level_penalty
                if best is None or score < best[0]:
                    best = (score, candidate_idx, entry_idx, safe)

            if best is None:
                loop_idx = min(remaining, key=lambda idx: dist2(current_point, contours[idx].centroid))
                loop = contours[loop_idx]
                entry_idx, _, safe = sampled_connection_to_loop(
                    model,
                    current_point,
                    loop,
                    required_clearance=required_clearance,
                    spacing=spacing,
                )
                safe = safe and connector_clear_of_path(
                    path,
                    [current_point, loop.points[entry_idx]],
                    spacing=spacing,
                    spacing_tolerance=spacing_tolerance,
                )
            else:
                _, loop_idx, entry_idx, safe = best
                loop = contours[loop_idx]

            entry_point = loop.points[entry_idx]
            if not safe:
                connector = grid_astar_connector(grid, current_point, entry_point, required_clearance=required_clearance)
                if connector is not None and connector_clear_of_path(path, connector, spacing, spacing_tolerance):
                    append_points(path, connector)
                else:
                    unsafe_connectors += 1

        gap = seam_gap_points(loop, spacing)
        append_points(path, almost_full_loop(loop.points, entry_idx, direction, gap))
        current_point = path[-1]
        last_area = loop.area
        remaining.remove(loop_idx)

    return path, unsafe_connectors


def scanline_intervals(
    model: PolygonModel,
    fixed: float,
    variable_min: float,
    variable_max: float,
    line_width: float,
    sample_step: float,
    vertical: bool,
) -> list[tuple[float, float]]:
    required = line_width * 0.5
    samples = max(8, int(math.ceil((variable_max - variable_min) / max(sample_step, EPS))))
    inside = []
    for i in range(samples + 1):
        v = variable_min + (variable_max - variable_min) * i / samples
        p = (fixed, v) if vertical else (v, fixed)
        inside.append(model.sdf(p) >= required)

    intervals: list[tuple[float, float]] = []
    start: float | None = None
    for i, is_inside in enumerate(inside):
        v = variable_min + (variable_max - variable_min) * i / samples
        if is_inside and start is None:
            start = v
        elif not is_inside and start is not None:
            prev = variable_min + (variable_max - variable_min) * (i - 1) / samples
            if prev - start >= line_width:
                intervals.append((start, prev))
            start = None
    if start is not None:
        if variable_max - start >= line_width:
            intervals.append((start, variable_max))
    return intervals


def build_monotone_scanline_path(
    model: PolygonModel,
    grid: SDFGrid,
    line_width: float,
    spacing: float,
    vertical: bool = True,
) -> tuple[list[Point], bool]:
    min_x, min_y, max_x, max_y = model.bounds(margin=0.0)
    fixed_min, fixed_max = (min_x + line_width * 0.5, max_x - line_width * 0.5) if vertical else (min_y + line_width * 0.5, max_y - line_width * 0.5)
    variable_min, variable_max = (min_y + line_width * 0.5, max_y - line_width * 0.5) if vertical else (min_x + line_width * 0.5, max_x - line_width * 0.5)
    if fixed_max <= fixed_min or variable_max <= variable_min:
        return [], False

    count = max(2, int(math.floor((fixed_max - fixed_min) / spacing)) + 1)
    path: list[Point] = []
    reverse = False
    unsupported_multi_interval = False

    for i in range(count):
        fixed = fixed_min + i * spacing
        intervals = scanline_intervals(
            model,
            fixed,
            variable_min,
            variable_max,
            line_width=line_width,
            sample_step=spacing * 0.25,
            vertical=vertical,
        )
        if not intervals:
            continue
        if len(intervals) > 1:
            unsupported_multi_interval = True
            intervals = [max(intervals, key=lambda interval: interval[1] - interval[0])]

        a, b = intervals[0]
        segment = [(fixed, b), (fixed, a)] if reverse else [(fixed, a), (fixed, b)]
        if not vertical:
            segment = [(p[1], p[0]) for p in segment]
        if path and not segment_is_inside(model, path[-1], segment[0], line_width * 0.5 - spacing * 0.15, spacing * 0.5):
            connector = grid_astar_connector(grid, path[-1], segment[0], required_clearance=line_width * 0.5 - spacing * 0.15)
            if connector is None:
                unsupported_multi_interval = True
            else:
                append_points(path, connector)
        append_points(path, segment)
        reverse = not reverse

    return path, not unsupported_multi_interval


def build_two_band_scanline_path(
    model: PolygonModel,
    grid: SDFGrid,
    line_width: float,
    spacing: float,
) -> list[Point]:
    min_x, min_y, max_x, max_y = model.bounds(margin=0.0)
    if model.holes:
        total_points = sum(len(hole) for hole in model.holes)
        y_cut = sum(p[1] for hole in model.holes for p in hole) / max(1, total_points)
    else:
        y_cut = (min_y + max_y) * 0.5

    fixed_values: list[float] = []
    x = min_x + line_width * 0.5
    while x <= max_x - line_width * 0.5 + EPS:
        fixed_values.append(x)
        x += spacing

    lower: list[tuple[float, float, float]] = []
    upper: list[tuple[float, float, float]] = []
    for fixed in fixed_values:
        intervals = scanline_intervals(
            model,
            fixed,
            min_y + line_width * 0.5,
            max_y - line_width * 0.5,
            line_width=line_width,
            sample_step=spacing * 0.25,
            vertical=True,
        )
        low_parts: list[tuple[float, float]] = []
        high_parts: list[tuple[float, float]] = []
        for a, b in intervals:
            if a < y_cut - spacing * 0.5:
                low_parts.append((a, min(b, y_cut - spacing * 0.5)))
            if b > y_cut + spacing * 0.5:
                high_parts.append((max(a, y_cut + spacing * 0.5), b))
        if low_parts:
            a = min(part[0] for part in low_parts)
            b = max(part[1] for part in low_parts)
            if b - a >= line_width:
                lower.append((fixed, a, b))
        if high_parts:
            a = min(part[0] for part in high_parts)
            b = max(part[1] for part in high_parts)
            if b - a >= line_width:
                upper.append((fixed, a, b))

    path: list[Point] = []
    reverse = False
    required_clearance = line_width * 0.5 - spacing * 0.15

    def append_segment_with_routed_connector(segment: list[Point]) -> None:
        if path and not segment_is_inside(model, path[-1], segment[0], required_clearance, spacing * 0.5):
            connector = grid_astar_connector(grid, path[-1], segment[0], required_clearance=required_clearance)
            if connector is not None:
                append_points(path, connector)
        append_points(path, segment)

    for fixed, a, b in lower:
        append_segment_with_routed_connector([(fixed, b), (fixed, a)] if reverse else [(fixed, a), (fixed, b)])
        reverse = not reverse

    reverse = False
    for fixed, a, b in reversed(upper):
        append_segment_with_routed_connector([(fixed, a), (fixed, b)] if reverse else [(fixed, b), (fixed, a)])
        reverse = not reverse

    return path


def point_line_distance(p: Point, a: Point, b: Point) -> float:
    return point_segment_distance(p, a, b)


def simplify_polyline(points: Sequence[Point], tolerance: float) -> list[Point]:
    if len(points) <= 2:
        return list(points)

    def simplify_range(start: int, end: int, out: list[Point]) -> None:
        max_dist = -1.0
        max_idx = start
        for idx in range(start + 1, end):
            d = point_line_distance(points[idx], points[start], points[end])
            if d > max_dist:
                max_dist = d
                max_idx = idx
        if max_dist > tolerance:
            simplify_range(start, max_idx, out)
            out.pop()
            simplify_range(max_idx, end, out)
        else:
            out.append(points[start])
            out.append(points[end])

    simplified: list[Point] = []
    simplify_range(0, len(points) - 1, simplified)
    deduped: list[Point] = []
    append_points(deduped, simplified)
    return deduped


def build_safe_chained_contours(
    model: PolygonModel,
    grid: SDFGrid,
    contours: Sequence[ContourLoop],
    anchor: Point,
    line_width: float,
    spacing: float,
) -> tuple[list[Point], int]:
    """Connect all contour loops into one stroke.

    This handles holes and branch-heavy cases in the prototype. It is deliberately
    conservative: each contour loop is covered as a loop, and the next loop is
    reached by the shortest sampled connector that stays in the printable domain.
    Production code should replace this with a real CFS contour-tree traversal,
    but this version is useful because it exposes containment and speed problems.
    """

    if not contours:
        return [], 0

    ordered = sorted(range(len(contours)), key=lambda idx: (contours[idx].level_index, -contours[idx].area))
    path: list[Point] = []
    unsafe_connectors = 0
    current_point = anchor
    direction = 1

    for order_pos, current_idx in enumerate(ordered):
        loop = contours[current_idx]
        start_point_idx, _, safe = sampled_connection_to_loop(
            model,
            current_point,
            loop,
            required_clearance=line_width * 0.5 - spacing * 0.15,
            spacing=spacing,
        )
        if path and not safe:
            connector = grid_astar_connector(
                grid,
                current_point,
                loop.points[start_point_idx],
                required_clearance=line_width * 0.5 - spacing * 0.15,
            )
            if connector is None:
                unsafe_connectors += 1
            else:
                append_points(path, connector)

        append_points(path, full_loop(loop.points, start_point_idx, direction))
        current_point = path[-1]
        direction *= -1

    return path, unsafe_connectors


@dataclass
class ValidationMetrics:
    path_points: int
    path_length: float
    contour_count: int
    contour_levels: int
    max_segment: float
    containment_violations: int
    clearance_warnings: int
    min_sdf: float
    self_intersections: int
    spacing_violations: int
    min_nonlocal_spacing: float
    start_outer_distance: float
    end_outer_distance: float
    start_on_outer: bool
    end_on_outer: bool
    unsafe_connectors: int
    coverage_ratio: float
    underfill_ratio: float
    overfill_ratio: float
    elapsed_seconds: float
    strategy: str
    ok: bool


def raster_coverage(
    model: PolygonModel,
    path: Sequence[Point],
    bounds: tuple[float, float, float, float],
    line_width: float,
    cells: int,
) -> tuple[float, float, float]:
    min_x, min_y, max_x, max_y = bounds
    width = max_x - min_x
    height = max_y - min_y
    if width <= EPS or height <= EPS:
        return 0.0, 1.0, 0.0
    if width >= height:
        nx = cells
        ny = max(12, int(round(cells * height / width)))
    else:
        ny = cells
        nx = max(12, int(round(cells * width / height)))
    cell = max(width / nx, height / ny)
    radius_cells = max(1, int(math.ceil((line_width * 0.5) / cell)) + 1)

    covered = bytearray(nx * ny)

    def mark(p: Point) -> None:
        ix = int((p[0] - min_x) / width * nx)
        iy = int((p[1] - min_y) / height * ny)
        if ix < -radius_cells or ix >= nx + radius_cells or iy < -radius_cells or iy >= ny + radius_cells:
            return
        for y in range(max(0, iy - radius_cells), min(ny, iy + radius_cells + 1)):
            cy = min_y + (y + 0.5) * height / ny
            for x in range(max(0, ix - radius_cells), min(nx, ix + radius_cells + 1)):
                cx = min_x + (x + 0.5) * width / nx
                if dist((cx, cy), p) <= line_width * 0.5 + cell * 0.75:
                    covered[y * nx + x] = 1

    for a, b in zip(path, path[1:]):
        length = dist(a, b)
        steps = max(1, int(math.ceil(length / max(cell * 0.45, EPS))))
        for i in range(steps + 1):
            mark(lerp(a, b, i / steps))

    inside_count = 0
    inside_covered = 0
    outside_covered = 0
    outside_count = 0
    for y in range(ny):
        cy = min_y + (y + 0.5) * height / ny
        for x in range(nx):
            cx = min_x + (x + 0.5) * width / nx
            idx = y * nx + x
            inside = model.contains((cx, cy))
            if inside:
                inside_count += 1
                if covered[idx]:
                    inside_covered += 1
            else:
                outside_count += 1
                if covered[idx]:
                    outside_covered += 1

    coverage_ratio = inside_covered / inside_count if inside_count else 0.0
    underfill_ratio = 1.0 - coverage_ratio
    overfill_ratio = outside_covered / max(1, inside_count)
    return coverage_ratio, underfill_ratio, overfill_ratio


def path_pair_metrics(
    path: Sequence[Point],
    spacing: float,
    spacing_tolerance: float,
    local_skip: int = 6,
) -> tuple[int, int, float]:
    if len(path) < 4:
        return 0, 0, float("inf")

    min_allowed = spacing * (1.0 - spacing_tolerance)
    bin_size = max(spacing * 1.5, EPS)
    bins: dict[tuple[int, int], list[int]] = {}
    segments = list(zip(path, path[1:]))
    segment_lengths = [dist(a, b) for a, b in segments]
    prefix_lengths = [0.0]
    for length in segment_lengths:
        prefix_lengths.append(prefix_lengths[-1] + length)
    path_length = prefix_lengths[-1]
    local_skip_distance = spacing * 4.0
    closed_path = dist(path[0], path[-1]) <= spacing * 0.20

    def key(x: float, y: float) -> tuple[int, int]:
        return (math.floor(x / bin_size), math.floor(y / bin_size))

    for idx, (a, b) in enumerate(segments):
        min_x = min(a[0], b[0]) - min_allowed
        min_y = min(a[1], b[1]) - min_allowed
        max_x = max(a[0], b[0]) + min_allowed
        max_y = max(a[1], b[1]) + min_allowed
        k0 = key(min_x, min_y)
        k1 = key(max_x, max_y)
        for by in range(k0[1], k1[1] + 1):
            for bx in range(k0[0], k1[0] + 1):
                bins.setdefault((bx, by), []).append(idx)

    checked: set[tuple[int, int]] = set()
    self_intersections = 0
    spacing_violations = 0
    min_nonlocal_spacing = float("inf")

    for i, (a, b) in enumerate(segments):
        min_x = min(a[0], b[0]) - min_allowed
        min_y = min(a[1], b[1]) - min_allowed
        max_x = max(a[0], b[0]) + min_allowed
        max_y = max(a[1], b[1]) + min_allowed
        k0 = key(min_x, min_y)
        k1 = key(max_x, max_y)
        candidates: set[int] = set()
        for by in range(k0[1], k1[1] + 1):
            for bx in range(k0[0], k1[0] + 1):
                candidates.update(bins.get((bx, by), []))

        for j in candidates:
            if j <= i:
                continue
            index_distance = abs(i - j)
            if closed_path:
                index_distance = min(index_distance, len(segments) - index_distance)
            if index_distance <= local_skip:
                continue
            path_gap = max(0.0, prefix_lengths[j] - prefix_lengths[i + 1])
            if closed_path:
                occupied_span = prefix_lengths[j + 1] - prefix_lengths[i]
                path_gap = min(path_gap, max(0.0, path_length - occupied_span))
            if path_gap <= local_skip_distance:
                continue
            pair = (i, j)
            if pair in checked:
                continue
            checked.add(pair)

            c, d = segments[j]
            pair_distance = segment_distance(a, b, c, d)
            min_nonlocal_spacing = min(min_nonlocal_spacing, pair_distance)

            if pair_distance <= 1e-7:
                self_intersections += 1
            elif pair_distance < min_allowed:
                spacing_violations += 1

    return self_intersections, spacing_violations, min_nonlocal_spacing


def validate_path(
    model: PolygonModel,
    contours: Sequence[ContourLoop],
    path: Sequence[Point],
    line_width: float,
    spacing: float,
    coverage_cells: int,
    coverage_threshold: float,
    spacing_tolerance: float,
    elapsed_seconds: float,
    strategy: str,
    unsafe_connectors: int,
) -> ValidationMetrics:
    max_segment = 0.0
    containment_violations = 0
    clearance_warnings = 0
    min_sdf = float("inf")
    outside_tolerance = -spacing * 0.05
    warning_clearance = line_width * 0.5 - spacing * 0.20

    for a, b in zip(path, path[1:]):
        length = dist(a, b)
        max_segment = max(max_segment, length)
        samples = max(2, int(math.ceil(length / max(spacing * 0.25, EPS))))
        for i in range(samples + 1):
            p = lerp(a, b, i / samples)
            sdf = model.sdf(p)
            min_sdf = min(min_sdf, sdf)
            if sdf < outside_tolerance:
                containment_violations += 1
            elif sdf < warning_clearance:
                clearance_warnings += 1

    coverage_ratio, underfill_ratio, overfill_ratio = raster_coverage(
        model,
        path,
        model.bounds(margin=line_width),
        line_width,
        coverage_cells,
    )
    self_intersections, spacing_violations, min_nonlocal_spacing = path_pair_metrics(
        path,
        spacing=spacing,
        spacing_tolerance=spacing_tolerance,
    )

    contour_levels = len({loop.level_index for loop in contours})
    outer_loop = max((loop for loop in contours if loop.level_index == 0), key=lambda loop: loop.area, default=None)
    if outer_loop is not None and path:
        start_outer_distance = distance_to_loop(path[0], outer_loop.points)
        end_outer_distance = distance_to_loop(path[-1], outer_loop.points)
    else:
        start_outer_distance = 0.0
        end_outer_distance = 0.0
    outer_tolerance = max(spacing * 0.20, line_width * 0.20)
    start_on_outer = bool(path) and start_outer_distance <= outer_tolerance
    end_on_outer = bool(path) and end_outer_distance <= outer_tolerance
    path_length = polyline_length(path, closed=False)
    ok = (
        bool(path)
        and len(path) >= 2
        and containment_violations == 0
        and self_intersections == 0
        and spacing_violations == 0
        and start_on_outer
        and end_on_outer
        and unsafe_connectors == 0
        and coverage_ratio >= coverage_threshold
    )

    return ValidationMetrics(
        path_points=len(path),
        path_length=path_length,
        contour_count=len(contours),
        contour_levels=contour_levels,
        max_segment=max_segment,
        containment_violations=containment_violations,
        clearance_warnings=clearance_warnings,
        min_sdf=min_sdf if min_sdf < float("inf") else 0.0,
        self_intersections=self_intersections,
        spacing_violations=spacing_violations,
        min_nonlocal_spacing=min_nonlocal_spacing if min_nonlocal_spacing < float("inf") else 0.0,
        start_outer_distance=start_outer_distance,
        end_outer_distance=end_outer_distance,
        start_on_outer=start_on_outer,
        end_on_outer=end_on_outer,
        unsafe_connectors=unsafe_connectors,
        coverage_ratio=coverage_ratio,
        underfill_ratio=underfill_ratio,
        overfill_ratio=overfill_ratio,
        elapsed_seconds=elapsed_seconds,
        strategy=strategy,
        ok=ok,
    )


@dataclass
class LayerResult:
    shape: str
    grid: SDFGrid
    contours: list[ContourLoop]
    path: list[Point]
    metrics: ValidationMetrics
    diagnostics: list[str]

    def to_json(self) -> dict[str, object]:
        return {
            "shape": self.shape,
            "grid": {
                "nx": self.grid.nx,
                "ny": self.grid.ny,
                "dx": self.grid.dx,
                "dy": self.grid.dy,
                "max_sdf": self.grid.max_sdf(),
            },
            "contours": {
                "count": len(self.contours),
                "levels": sorted({loop.level_index for loop in self.contours}),
                "per_level": {
                    str(level): sum(1 for loop in self.contours if loop.level_index == level)
                    for level in sorted({loop.level_index for loop in self.contours})
                },
            },
            "path": {
                "points": len(self.path),
                "start": list(self.path[0]) if self.path else None,
                "end": list(self.path[-1]) if self.path else None,
            },
            "metrics": self.metrics.__dict__,
            "diagnostics": self.diagnostics,
        }


def plan_one_layer(
    model: PolygonModel,
    grid_cells: int,
    line_width: float,
    spacing: float,
    max_levels: int,
    coverage_cells: int,
    coverage_threshold: float,
    spacing_tolerance: float,
    start_fraction: float,
    exit_fraction: float,
) -> LayerResult:
    started = time.perf_counter()
    grid = build_sdf_grid(model, grid_cells=grid_cells, margin=line_width * 2.0)
    contours = filter_topology_stable_levels(
        filter_printable_contours(
            generate_offset_contours(grid, line_width=line_width, spacing=spacing, max_levels=max_levels),
            spacing=spacing,
        ),
        spacing=spacing,
    )
    diagnostics: list[str] = []

    if not contours:
        elapsed = time.perf_counter() - started
        metrics = validate_path(
            model,
            [],
            [],
            line_width,
            spacing,
            coverage_cells,
            coverage_threshold,
            spacing_tolerance,
            elapsed,
            "none",
            unsafe_connectors=0,
        )
        diagnostics.append("No valid contours were generated. The shape is probably too narrow for this line width.")
        return LayerResult(model.name, grid, contours, [], metrics, diagnostics)

    outer_loop = max((loop for loop in contours if loop.level_index == 0), key=lambda loop: loop.area, default=contours[0])
    start_anchor = point_at_closed_fraction(outer_loop.points, start_fraction)
    exit_anchor = point_at_closed_fraction(outer_loop.points, exit_fraction)
    grouped = loops_by_level(contours)
    multi_loop = any(len(grouped[level]) > 1 for level in grouped)
    first_level = min(grouped) if grouped else 0
    one_hole_ring = multi_loop and len(grouped.get(first_level, [])) == 2 and centroid_spread(grouped[first_level]) < spacing * 2.0
    branch_single_island = multi_loop and len(grouped.get(first_level, [])) == 1
    if one_hole_ring:
        filtered_contours = filter_ring_medial_overlap(contours, spacing)
        if len(filtered_contours) != len(contours):
            diagnostics.append(
                f"Dropped {len(contours) - len(filtered_contours)} medial contour(s) whose opposing ring fronts were below nominal spacing."
            )
            contours = filtered_contours
            grouped = loops_by_level(contours)
            multi_loop = any(len(grouped[level]) > 1 for level in grouped)
            first_level = min(grouped) if grouped else 0
            one_hole_ring = (
                multi_loop and len(grouped.get(first_level, [])) == 2 and centroid_spread(grouped[first_level]) < spacing * 2.0
            )
            branch_single_island = multi_loop and len(grouped.get(first_level, [])) == 1

    ordered_contours = ordered_loops_for_spiral(contours)
    if branch_single_island:
        path, unsafe_connectors = build_branch_connected_fermat(
            contours,
            start_anchor=start_anchor,
            spacing=spacing,
        )
        strategy = "branch_connected_fermat"
    elif multi_loop and len(model.holes) >= 2:
        path, unsafe_connectors = build_merged_hole_connected_fermat(
            model,
            contours,
            start_anchor=start_anchor,
            line_width=line_width,
            spacing=spacing,
        )
        strategy = "merged_hole_connected_fermat"
    else:
        candidate_exits: list[tuple[float, Point]] = [(exit_fraction, exit_anchor)]
        for fraction in (0.67, 0.50, 0.33, 0.25, 0.75, 0.10, 0.90):
            if all(abs(((fraction - existing + 0.5) % 1.0) - 0.5) > 1e-3 for existing, _ in candidate_exits):
                candidate_exits.append((fraction, point_at_closed_fraction(outer_loop.points, fraction)))

        best_candidate: tuple[tuple[int, int, int, float], list[Point], int, float] | None = None
        for fraction, candidate_exit in candidate_exits:
            candidate_path, candidate_unsafe = build_single_minimum_connected_fermat(
                ordered_contours,
                start_anchor=start_anchor,
                spacing=spacing,
                port_spacing=spacing * 2.5,
                exit_anchor=candidate_exit,
                preserve_medial_pockets=not one_hole_ring,
            )
            crossings, close_pairs, min_spacing = path_pair_metrics(
                candidate_path,
                spacing=spacing,
                spacing_tolerance=spacing_tolerance,
            )
            score = (candidate_unsafe, crossings, close_pairs, -min_spacing)
            if best_candidate is None or score < best_candidate[0]:
                best_candidate = (score, candidate_path, candidate_unsafe, fraction)
            if candidate_unsafe == 0 and crossings == 0 and close_pairs == 0:
                best_candidate = (score, candidate_path, candidate_unsafe, fraction)
                break

        assert best_candidate is not None
        path = best_candidate[1]
        unsafe_connectors = best_candidate[2]
        if abs(((best_candidate[3] - exit_fraction + 0.5) % 1.0) - 0.5) > 1e-3:
            diagnostics.append(f"Selected boundary exit fraction {best_candidate[3]:.2f} to satisfy local spacing checks.")
        strategy = "single_minimum_connected_fermat"

    if one_hole_ring:
        diagnostics.append("One-hole topology is currently routed as a single contour chain; inspect contour-tree ordering.")
    elif branch_single_island:
        diagnostics.append("Branch topology used contour-tree Fermat child insertion.")
    elif multi_loop:
        diagnostics.append(
            "Multi-hole topology used a conservative merged-cavity CFS fallback; material between merged holes may be omitted."
        )

    completed_path = complete_outer_boundary_cycle(path, outer_loop, spacing, spacing_tolerance)
    if len(completed_path) != len(path) or (path and completed_path and dist(path[-1], completed_path[-1]) > EPS):
        diagnostics.append("Completed the outer-boundary cycle so start/end are adjacent at a tiny cut gap.")
    path = completed_path
    if multi_loop and len(model.holes) >= 2:
        pocket_path, inserted_pockets = insert_uncovered_pocket_spirals(
            model,
            contours,
            path,
            spacing=spacing,
            spacing_tolerance=spacing_tolerance,
        )
        if inserted_pockets:
            diagnostics.append(f"Inserted {inserted_pockets} uncovered pocket Fermat spiral(s) into the global path.")
            path = pocket_path

    elapsed = time.perf_counter() - started
    metrics = validate_path(
        model,
        contours,
        path,
        line_width,
        spacing,
        coverage_cells,
        coverage_threshold,
        spacing_tolerance,
        elapsed,
        strategy,
        unsafe_connectors,
    )

    if metrics.containment_violations:
        diagnostics.append(f"{metrics.containment_violations} sampled path points leave the printable polygon.")
    if metrics.self_intersections:
        diagnostics.append(f"{metrics.self_intersections} non-adjacent red path segment pair(s) cross or touch.")
    if metrics.spacing_violations:
        diagnostics.append(
            f"{metrics.spacing_violations} non-adjacent red path segment pair(s) are closer than "
            f"{spacing * (1.0 - spacing_tolerance):.3f}."
        )
    if not metrics.start_on_outer:
        diagnostics.append(f"Path start is {metrics.start_outer_distance:.3f} from the outer contour.")
    if not metrics.end_on_outer:
        diagnostics.append(f"Path end is {metrics.end_outer_distance:.3f} from the outer contour.")
    if metrics.clearance_warnings:
        diagnostics.append(
            f"{metrics.clearance_warnings} sampled path points are inside the polygon but below nominal half-width clearance."
        )
    if metrics.unsafe_connectors:
        diagnostics.append(f"{metrics.unsafe_connectors} connector(s) could not be proven safe by sampled clearance checks.")
    if metrics.coverage_ratio < coverage_threshold:
        diagnostics.append(
            f"Coverage {metrics.coverage_ratio:.3f} is below threshold {coverage_threshold:.3f}; inspect SVG."
        )

    return LayerResult(model.name, grid, contours, path, metrics, diagnostics)


def svg_points(points: Sequence[Point], bounds: tuple[float, float, float, float], width: int, height: int, pad: int) -> str:
    min_x, min_y, max_x, max_y = bounds
    scale = min((width - pad * 2) / max(max_x - min_x, EPS), (height - pad * 2) / max(max_y - min_y, EPS))

    def tx(p: Point) -> tuple[float, float]:
        x = pad + (p[0] - min_x) * scale
        y = height - pad - (p[1] - min_y) * scale
        return x, y

    return " ".join(f"{tx(p)[0]:.2f},{tx(p)[1]:.2f}" for p in points)


def write_svg(result: LayerResult, path: Path, draw_contours: bool) -> None:
    width = 1100
    height = 850
    pad = 32
    bounds = result.grid.model.bounds(margin=4.0)
    model = result.grid.model

    lines: list[str] = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#fbfbf8"/>',
    ]

    outer = svg_points(model.outer + [model.outer[0]], bounds, width, height, pad)
    lines.append(f'<polyline points="{outer}" fill="none" stroke="#242424" stroke-width="2.2"/>')
    for hole in model.holes:
        hole_points = svg_points(hole + [hole[0]], bounds, width, height, pad)
        lines.append(f'<polyline points="{hole_points}" fill="none" stroke="#555" stroke-width="1.8" stroke-dasharray="5 5"/>')

    if draw_contours:
        for loop in result.contours:
            color = "#9ab7d8" if loop.level_index % 2 == 0 else "#b7cda3"
            pts = svg_points(loop.points + [loop.points[0]], bounds, width, height, pad)
            lines.append(f'<polyline points="{pts}" fill="none" stroke="{color}" stroke-width="0.8" opacity="0.45"/>')

    if result.path:
        pts = svg_points(result.path, bounds, width, height, pad)
        lines.append(f'<polyline points="{pts}" fill="none" stroke="#d12f1f" stroke-width="2.0" stroke-linejoin="round" stroke-linecap="round"/>')

        def svg_point(p: Point) -> tuple[float, float]:
            min_x, min_y, max_x, max_y = bounds
            scale = min((width - pad * 2) / max(max_x - min_x, EPS), (height - pad * 2) / max(max_y - min_y, EPS))
            return (pad + (p[0] - min_x) * scale, height - pad - (p[1] - min_y) * scale)

        sx, sy = svg_point(result.path[0])
        ex, ey = svg_point(result.path[-1])
        lines.append(f'<circle cx="{sx:.2f}" cy="{sy:.2f}" r="5" fill="#12873f"><title>start</title></circle>')
        lines.append(f'<circle cx="{ex:.2f}" cy="{ey:.2f}" r="5" fill="#6e35b8"><title>end</title></circle>')

    status = "PASS" if result.metrics.ok else "FAIL"
    subtitle = (
        f"{result.shape}: {status}, strategy={result.metrics.strategy}, "
        f"coverage={result.metrics.coverage_ratio:.3f}, "
        f"cross={result.metrics.self_intersections}, close={result.metrics.spacing_violations}, "
        f"minsep={result.metrics.min_nonlocal_spacing:.3f}, "
        f"gap={dist(result.path[0], result.path[-1]) if result.path else 0.0:.3f}, "
        f"time={result.metrics.elapsed_seconds:.3f}s"
    )
    lines.append(f'<text x="32" y="32" font-family="monospace" font-size="18" fill="#222">{subtitle}</text>')
    lines.append("</svg>")
    path.write_text("\n".join(lines), encoding="utf-8")


def write_png_file(path: Path, width: int, height: int, pixels: bytearray) -> None:
    def chunk(kind: bytes, data: bytes) -> bytes:
        return struct.pack(">I", len(data)) + kind + data + struct.pack(">I", zlib.crc32(kind + data) & 0xFFFFFFFF)

    raw = bytearray()
    stride = width * 3
    for y in range(height):
        raw.append(0)
        raw.extend(pixels[y * stride : (y + 1) * stride])

    data = b"".join(
        [
            b"\x89PNG\r\n\x1a\n",
            chunk(b"IHDR", struct.pack(">IIBBBBB", width, height, 8, 2, 0, 0, 0)),
            chunk(b"IDAT", zlib.compress(bytes(raw), level=6)),
            chunk(b"IEND", b""),
        ]
    )
    path.write_bytes(data)


def write_png_preview(result: LayerResult, path: Path, draw_contours: bool) -> None:
    width = 1400
    height = 1000
    pad = 48
    bounds = result.grid.model.bounds(margin=4.0)
    min_x, min_y, max_x, max_y = bounds
    scale = min((width - pad * 2) / max(max_x - min_x, EPS), (height - pad * 2) / max(max_y - min_y, EPS))
    x_offset = (width - (max_x - min_x) * scale) * 0.5
    y_offset = (height - (max_y - min_y) * scale) * 0.5
    pixels = bytearray([251, 251, 248] * width * height)

    def tx(p: Point) -> tuple[float, float]:
        return (x_offset + (p[0] - min_x) * scale, height - y_offset - (p[1] - min_y) * scale)

    def set_pixel(x: int, y: int, color: tuple[int, int, int]) -> None:
        if 0 <= x < width and 0 <= y < height:
            idx = (y * width + x) * 3
            pixels[idx : idx + 3] = bytes(color)

    def draw_disc(x: float, y: float, radius: float, color: tuple[int, int, int]) -> None:
        r = int(math.ceil(radius))
        for yy in range(int(y) - r, int(y) + r + 1):
            for xx in range(int(x) - r, int(x) + r + 1):
                if (xx - x) * (xx - x) + (yy - y) * (yy - y) <= radius * radius:
                    set_pixel(xx, yy, color)

    def draw_polyline(points: Sequence[Point], color: tuple[int, int, int], radius: float, closed: bool = False) -> None:
        if len(points) < 2:
            return
        pairs = list(zip(points, points[1:]))
        if closed:
            pairs.append((points[-1], points[0]))
        for a, b in pairs:
            ax, ay = tx(a)
            bx, by = tx(b)
            length = math.hypot(bx - ax, by - ay)
            steps = max(1, int(math.ceil(length / max(radius * 0.65, 1.0))))
            for i in range(steps + 1):
                t = i / steps
                draw_disc(ax + (bx - ax) * t, ay + (by - ay) * t, radius, color)

    if draw_contours:
        for loop in result.contours:
            color = (154, 183, 216) if loop.level_index % 2 == 0 else (183, 205, 163)
            draw_polyline(loop.points, color, radius=0.75, closed=True)

    draw_polyline(result.grid.model.outer, (36, 36, 36), radius=1.6, closed=True)
    for hole in result.grid.model.holes:
        draw_polyline(hole, (82, 82, 82), radius=1.3, closed=True)
    draw_polyline(result.path, (209, 47, 31), radius=1.8, closed=False)

    if result.path:
        sx, sy = tx(result.path[0])
        ex, ey = tx(result.path[-1])
        draw_disc(sx, sy, 6.0, (18, 135, 63))
        draw_disc(ex, ey, 6.0, (110, 53, 184))

    write_png_file(path, width, height, pixels)


def write_json(result: LayerResult, path: Path) -> None:
    path.write_text(json.dumps(result.to_json(), indent=2), encoding="utf-8")


def print_summary(result: LayerResult) -> None:
    m = result.metrics
    status = "PASS" if m.ok else "FAIL"
    print(
        f"{result.shape:12s} {status:4s} "
        f"strategy={m.strategy:24s} contours={m.contour_count:3d} levels={m.contour_levels:2d} "
        f"points={m.path_points:5d} coverage={m.coverage_ratio:5.3f} "
        f"cross={m.self_intersections:4d} close={m.spacing_violations:4d} "
        f"minsep={m.min_nonlocal_spacing:5.3f} viol={m.containment_violations:4d} "
        f"gap={dist(result.path[0], result.path[-1]) if result.path else 0.0:5.3f} "
        f"outer=({m.start_outer_distance:4.2f},{m.end_outer_distance:4.2f}) "
        f"clear={m.clearance_warnings:4d} unsafe={m.unsafe_connectors:2d} "
        f"time={m.elapsed_seconds:6.3f}s"
    )
    for diagnostic in result.diagnostics:
        print(f"  - {diagnostic}")


def parse_args(argv: Sequence[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--shape", choices=sorted(available_shapes()), help="single built-in shape to slice")
    parser.add_argument("--all", action="store_true", help="run every built-in shape")
    parser.add_argument("--list-shapes", action="store_true", help="list built-in shape names and exit")
    parser.add_argument("--out", type=Path, default=Path("build") / "continuous_fermat", help="output directory")
    parser.add_argument("--grid", type=int, default=180, help="SDF grid cells along the longest axis")
    parser.add_argument("--coverage-grid", type=int, default=170, help="coverage verifier grid cells along the longest axis")
    parser.add_argument("--line-width", type=float, default=1.2, help="extrusion line width")
    parser.add_argument("--spacing", type=float, default=1.2, help="spacing between contour centerlines")
    parser.add_argument("--max-levels", type=int, default=256, help="maximum offset contour levels")
    parser.add_argument("--coverage-threshold", type=float, default=0.82, help="minimum filled-area estimate")
    parser.add_argument("--start-fraction", type=float, default=0.0, help="outer contour fraction for the layer cycle cut/start")
    parser.add_argument("--exit-fraction", type=float, default=0.5, help="outer contour fraction for the second boundary slot")
    parser.add_argument(
        "--spacing-tolerance",
        type=float,
        default=0.25,
        help="allowed fractional undershoot before non-adjacent red lines count as too close",
    )
    parser.add_argument("--draw-contours", action="store_true", help="include raw offset contours in the SVG")
    parser.add_argument("--no-svg", action="store_true", help="skip SVG writing")
    parser.add_argument("--no-png", action="store_true", help="skip PNG preview writing")
    return parser.parse_args(argv)


def main(argv: Sequence[str]) -> int:
    args = parse_args(argv)
    shapes = available_shapes()

    if args.list_shapes:
        for name in sorted(shapes):
            print(name)
        return 0

    if args.all or not args.shape:
        selected = sorted(shapes)
    else:
        selected = [args.shape]

    args.out.mkdir(parents=True, exist_ok=True)
    any_failed = False

    for name in selected:
        result = plan_one_layer(
            shapes[name],
            grid_cells=args.grid,
            line_width=args.line_width,
            spacing=args.spacing,
            max_levels=args.max_levels,
            coverage_cells=args.coverage_grid,
            coverage_threshold=args.coverage_threshold,
            spacing_tolerance=args.spacing_tolerance,
            start_fraction=args.start_fraction,
            exit_fraction=args.exit_fraction,
        )
        print_summary(result)
        write_json(result, args.out / f"{name}.json")
        if not args.no_svg:
            write_svg(result, args.out / f"{name}.svg", draw_contours=args.draw_contours)
        if not args.no_png:
            write_png_preview(result, args.out / f"{name}.png", draw_contours=args.draw_contours)
        any_failed = any_failed or not result.metrics.ok

    return 2 if any_failed else 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
