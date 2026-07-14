# Continuous Slicing Development Audit

## Current status

Continuous slicing now preserves and exports structurally valid paths even when geometric quality checks raise critical warnings. Those warnings are shown in Orca, stored on the protected extrusion path, and embedded in G-code as `_CONTINUOUS_FERMAT_VALIDATION_WARNING` comments. Machine-envelope, finite-motion, metadata, physical-width, and volumetric-flow checks remain hard failures. Cached or imported Continuous Fermat entities remain blocked because their geometry/metadata pair cannot be recertified against current slices.

This is experimental validation, not a claim that arbitrary hardware is physically certified. The software cannot infer every printer's home position, clips, toolhead envelope, leveling behavior, material response, or mechanical condition.

## User-facing eligibility contract

For every non-empty model layer:

1. there is exactly one connected island; holes do not create additional islands;
2. layer zero starts at the build plate; and
3. each later layer has positive-area geometric overlap with the immediately lower layer.

Layers need not be identical. Translation, rotation, taper, growth, shrinkage, and outline changes are admitted. The emitter may perform one non-extruding in-polygon XY approach after the positive Z transition when adjacent layers choose different seams.

Eligibility is intentionally broader than quality certification. An eligible layer is unsupported only when no structurally bounded path can be emitted. Containment, crossings, turns, coverage, outside-model area, finite material deviation, and redeposition remain visible quality diagnostics but do not discard a previewable path. This distinction is necessary because the eligibility rules alone allow arbitrarily narrow features or arbitrarily small interlayer overlap.

## Algorithm changes from the restrictive baseline

- Removed exact cross-section equality and vertical-prism admission.
- Replaced same-XY interlayer handoff with positive Z followed by at most one non-extruding XY approach.
- Removed the blanket printer-ready export interlock; structural and machine validators remain authoritative.
- Added `continuous_max_line_width`, user-configurable up to the immutable two-nozzle ceiling.
- Derived physical width from rounded-rectangle bead area and slowed widened segments so volumetric deposition remains consistent.
- Replaced pairwise bead-overlap rejection with exact unioned-footprint and global redeposition checks; ordinary adjacent bead overlap is expected.
- Added flow-relative contour phases only as bounded fallbacks.
- Added one geometry-derived width retry for safe exact/material underfill: the rounded-corner bead model converts target area per centerline length into a physical width, and the complete offsets are regenerated at that width and matching spacing before unchanged validation.
- Added footprint-directed width transfer: material is removed from nonlocal overlapping donor segments and assigned to segments adjacent to exact uncovered regions while total modeled material remains constant.
- Added validator-guided 2-opt uncrossing, touching-connector simplification, and local sharp-turn shortcut repair.
- Enabled connected islands containing holes when their final path passes exact hole containment and footprint validation.
- Preserved every contour-tree branch, including branches born after level zero, and constrained branch exits to replace only a local parent arc.
- Added narrow-feature centerlines for cone and pyramid tips.
- Replaced feature-size-relaxed quality thresholds with fixed metrics and a deterministic 2.5-line-width erosion category for thin sections, necks, and hole webs.
- Split validation into structurally emittable and geometric-quality results; advisory results survive preview, JSON cache serialization, and G-code generation.
- Bounded dense-path topology work: oversampled routes above 5,000 points are simplified to 0.01 mm before compact repair, cutting the measured 37.2-second stress layer to about 2.4 seconds while retaining explicit quality diagnostics.

No model index, shape name, or fixture-specific exception exists in the planner.

## Structural and advisory geometry contract

Every emitted layer has a closed path, matching finite metadata, bounded physical width, and a non-empty extrusion footprint. The fixed targets below are advisory:

- one exactly closed polyline and matching finite segment metadata;
- physical width between the continuous-flow floor (5% nozzle diameter or the rounded-corner floor) and the lower of the configured ceiling or 2.0 nozzle diameters;
- zero centerline containment violations;
- zero nonlocal crossings;
- zero physically relevant turns sharper than 135 degrees;
- exact polygon-area coverage at least 0.980; 70-cell sampled coverage remains diagnostic;
- outside swept area at most 0.020 of printable area;
- deposited material from 0.980 through 1.020 of solid-layer material; and
- global redeposition at most 0.040 of printable area.

Pairwise close-spacing and bead-overlap counts are retained for diagnosis but are not hard failures. Exact unioned area, material, and global redeposition measure the actual safety/quality properties those proxies attempted to protect.

Normal entities are cleared only after the final generated path is structurally emittable. A quality miss creates a critical warning and embedded diagnostics. Protected paths are excluded from later generic simplification, clipping, arc fitting, flow boosts, small-area compensation, overhang resampling, and resonance speed changes.

## Serialized-output contract

The strict emitter:

- uses deterministic Marlin 2 or Klipper state normalization;
- establishes absolute XYZ, relative E, unity overrides, disabled pressure advance, homing, and heater-ready state;
- omits arbitrary profile start/end, layer-change, timelapse, and extrusion-role hooks;
- emits exactly one marked closed extrusion section per nondegenerate layer, allowing only one final zero-area apex omission;
- permits only positive Z plus at most one non-extruding XY approach between sections;
- suppresses retraction, wipe, Z-hop, object cancellation, power-loss recovery, pressure equalization, firmware-limit output, and auxiliary extrusion features;
- checks quantized nonzero XY, complete variable-width footprint containment in the transformed bed polygon, bounded positive Z, quantized-E material, and a downward-rounded volumetric feed cap; and
- uses a no-XY shutdown with at most one bounded final Z lift.

The Python auditor independently parses this structure for Marlin 2 and Klipper and has adversarial tests for modal tricks, malformed markers, arcs, non-finite values, unmarked extrusion, retraction, invalid transitions, and unsafe postambles. It does not duplicate C++ geometry/material calculations.

## Measured evidence (2026-07-12)

The deterministic `[corpus-100]` Catch2 regression contains 100 distinct connected, hole-free extruded-part footprints. Vertex count, aspect ratio, lobe count, modulation, radius, and phase vary independently. With 1.2 mm nominal line/nozzle width and a 2.4 mm maximum continuous width:

| Iteration | Accepted | Principal remaining failures |
| --- | ---: | --- |
| Initial strict baseline | 33/100 | coverage, crossings, pair-overlap proxy |
| Exact footprint + global redeposition | 52/100 | residual offset-phase gaps |
| Dense candidate and targeted widening | 82/100 | redeposition/coverage and crossings |
| 2-opt topology repair | 83/100 | 12 redeposition, 3 crossing, 2 coverage |
| Bounded offset phases | 88/100 | 9 redeposition, 2 coverage, 1 crossing |
| Footprint-directed width transfer | 96/100 | 2 redeposition, 1 crossing, 1 turnback |
| Iterated transfer and local topology repair | **100/100** | none |

That table records the earlier all-quality optimization track. The current mandatory corpus assertion is structural emittability (`accepted == 100`); quality misses remain separately counted and visible. A separate isolated touching-connector regression preserves the formerly last failing topology.

Additional automated evidence includes:

- common 0.4 mm nozzle and variable-width rectangle regressions;
- an eight-shape sharp/hole corpus (acute triangles, wedge, diamond, star, annulus, square hole, and two holes), structurally emittable 8/8 with fixed quality diagnostics retained;
- open, underfilled, self-crossing, sub-resolution crossing, hole-excursion, hairpin, metadata, and width-bound failures;
- a 200-layer tapering pyramid exported successfully through the complete printer-ready G-code path, plus pairwise-layer-overlap FFF integration;
- Marlin 2/Klipper printer-ready serialization, quantized material, flow calibration, volumetric cap, and forbidden-feature suppression; and
- 25 Python serialized-G-code auditor tests.

## Remaining limitations

- The complete layer uses external-perimeter flow/role semantics rather than bridge, top, bottom, or infill-specific behavior.
- The job is currently limited to one object, instance, region/material, configured nozzle, and extruder. Multiple regions may require flow/tool changes inside the one-stroke contract.
- The rounded-rectangle bead model does not simulate pressure advance, acceleration corner effects, die swell, cooling, adhesion, or extrusion force.
- The 100-part corpus is broad deterministic procedural evidence, not proof over all possible connected polygons or all nozzle/material combinations. The internet corpus separately records `passed` (emittable) and `quality_ok`.
- Runtime on difficult fallback cases remains higher because contour construction and exact capsule unions are repeated, although dense-path compact repair removed the largest observed superlinear hotspot.
- The configured printable polygon is an extrusion envelope, not a complete collision-free travel envelope. Generic `G28` and the first approach remain machine-dependent.
- Physical qualification still requires a reviewed printer/profile, conservative settings, supervised trials, and independent emergency-stop readiness.

## Research direction

The implementation is inspired by Connected Fermat Spirals, but a complete general solution should still move toward deterministic region decomposition and graph traversal:

1. decompose connected geometry into spiral-fillable subregions;
2. reserve printable connector corridors;
3. generate low-curvature fills and connect them through a graph without crossings;
4. optimize variable width against exact uncovered and multiply-covered area; and
5. preserve structural and serialized machine checks as hard authority while retaining geometric metrics as reviewable advisories.
