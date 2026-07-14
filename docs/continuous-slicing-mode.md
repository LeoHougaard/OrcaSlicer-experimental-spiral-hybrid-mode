# Continuous Slicing Mode (Experimental)

Continuous slicing (`spiral_hybrid_non_crossing`) replaces every eligible non-empty layer with one closed, continuously extruded, variable-width path. It is independent of classic Spiral vase mode. Structurally valid output remains available for preview and export when geometric quality checks raise warnings; there is no blanket research-mode export interlock.

This remains experimental software, not physical-printer certification. The slicer validates the geometry and serialized commands it can model, but cannot prove a machine's homing behavior, collision clearances, leveling state, material response, or hardware condition.

## Eligibility

A model is geometrically eligible when:

- every non-empty layer has exactly one connected island (the island may contain holes);
- layer zero begins at the build plate; and
- every later layer has positive-area geometric overlap with the layer immediately below.

Cross-sections may translate, rotate, taper, grow, shrink, or change outline. Identical layers and vertical prisms are not required. An open path, corrupt segment metadata, non-finite flow, invalid physical width, or empty extrusion footprint remains unsupported. Coverage, model-relative containment, crossings, turnbacks, redeposition, and finite material-ratio deviations are instead reported as critical geometric-quality warnings so the path can still be inspected.

The current whole-job contract also requires one object, one instance, one print region/material/extruder, tool 0, by-layer printing, relative E, and either Marlin 2 or Klipper semantics. Support, raft, skirt, draft shield, brim, prime tower, ironing, per-layer custom G-code, and calibration modes are disabled or rejected because they would introduce additional extrusion strokes or state changes. Cached/imported Continuous Fermat entities remain rejected because they cannot be revalidated against the current slices.

## Enabling and line width

```ini
spiral_mode = 0
spiral_hybrid_non_crossing = 1
```

**Maximum continuous line width** (`continuous_max_line_width`) controls how far the planner may widen individual segments to close residual gaps. It defaults to 150% of nozzle diameter and is hard-clamped to 200%. Width is derived from the rounded-rectangle bead cross-section, not by multiplying the nominal width directly.

Wider segments retain the same volumetric deposition by slowing feed rate by the corresponding extrusion multiplier. Narrower segments are not accelerated beyond nominal speed. The hard maximum is 2.0 nozzle diameters and the configured ceiling. Narrow terminal features may reduce flow down to the larger of 5% nozzle diameter or the rounded-corner cross-section floor.

The planner uses variable width in three stages:

1. when a safe nominal path underfills both exact footprint and material, it derives one physical bead width from target area and centerline length, then regenerates the complete contour geometry at that width and matching effective spacing;
2. it directs available material toward exact uncovered regions; and
3. for difficult residual phases, it transfers material from nonlocal overlapping segments to nearby uncovered regions while keeping total modeled material constant.

All adjusted candidates are rechecked by the same fixed quality metrics. Once an emittable path exists, the planner performs at most one geometry-derived width retry for safe underfill rather than delaying preview with unbounded spacing sweeps.

## Generated path and layer transitions

Each nondegenerate layer contains one exactly closed marked section:

```gcode
;_CONTINUOUS_FERMAT_BEGIN
...
;_CONTINUOUS_FERMAT_END
```

An advisory layer is annotated immediately before or inside its section:

```gcode
;_CONTINUOUS_FERMAT_VALIDATION_WARNING layer=42 ... exact_coverage=0.91 ...
```

Extrusion stops between layers. The emitter first performs a finite positive Z move. If the next layer's selected seam differs, it then performs at most one non-extruding XY approach inside the configured printable polygon before beginning the next closed section. Retraction, wipe, Z-hop, arbitrary layer-change hooks, and unmarked extrusion are not allowed in this transition.

This is one closed extrusion section per layer, not a continuously extruding three-dimensional helix.

## Structural requirements and geometric advisories

Every emitted layer must have:

- exactly closed path with matching per-segment metadata;
- finite extrusion multipliers and physical width within the configured/nozzle limits;
- a non-empty physical swept footprint; and
- finite positive serialized motion and extrusion within the configured machine envelope.

The following fixed targets produce a critical warning rather than deleting an otherwise emittable path: zero model-relative containment violations, zero nonlocal crossings, zero physically relevant sharp turnbacks, exact swept-bead coverage at least 0.980, outside-model swept area at most 0.020, material ratio from 0.980 through 1.020, and redeposition at most 0.040. Coarse sampled coverage remains diagnostic only. A final zero-area apex may be omitted when every preceding layer emitted exactly one protected stroke.

Failed geometry is additionally categorized at a 2.5-nominal-line-width threshold. Half-threshold erosion identifies globally thin sections, local necks, and hole webs. This category explains a warning; it never pre-rejects a thin layer that generated a satisfactory path.

Pairwise bead-overlap and close-spacing counts remain diagnostics. Adjacent or locally neighboring beads normally overlap, so those counts are not independent rejection conditions; the exact unioned footprint, material, and global redeposition checks are authoritative.

The emitter additionally checks transformed/quantized XY, the complete variable-width bead footprint against the configured bed polygon, positive bounded Z, non-collapsed serialized segments, and a conservatively rounded volumetric feed cap. A non-finite serialized material ratio is a hard failure; a finite ratio outside the target is embedded as another advisory comment.

## Printer-ready export

Ordinary file export proceeds when the job passes the structural, configuration, machine-envelope, and G-code checks used by this mode. Geometric-quality warnings are shown in Orca and retained in the G-code so the user can inspect the preview and decide whether to use the result. Profile start/end and layer-change macros are omitted in favor of a deterministic firmware-specific sequence; configured pressure advance is reset to zero, firmware machine-limit output is suppressed, and the shutdown sequence does not park in XY.

Marlin 2 normalization includes tool selection, relative E, unity overrides, disabled volumetric E/linear advance/power-loss recovery, all-axis homing, and a cooling-aware nozzle wait. Klipper uses its native relative-E, unity override, zero-pressure-advance, homing, bed-state, and heater-wait commands.

## Independent G-code audit

The Python auditor is an additional manual check, not the export gate:

```powershell
python tools\continuous_fermat\validate_gcode.py path\to\model.gcode `
    --expected-layers 120 --expected-sections 120 --firmware marlin2
```

It checks marker/layer accounting, modal state, finite coordinates, known positive feed and Z, positive-E XY motion inside sections, exact section closure within the configured serialized tolerance, one positive Z transition plus at most one non-extruding XY approach between sections, and a tightly restricted postamble. It does not recompute the C++ swept-footprint or material calculations.

## Regression evidence and limitations

The deterministic production-planner corpus contains 100 distinct connected, hole-free extruded-part footprints spanning vertex counts, aspect ratios, lobes, phases, and concavity. All 100 produce structurally emittable paths. An additional eight-shape corpus covers acute triangles, a wedge, diamond, star, annulus, square hole, and two holes; all eight are emittable, with any remaining geometric defects recorded as advisories. A 200-layer pyramid regression verifies tapering sharp points through printer-ready G-code. Separate tests cover concave/branched shapes, advisory export, quantized material/flow, malformed paths, machine-envelope containment, and G-code modal attacks.

The two eligibility statements alone do not mathematically guarantee physical printability for arbitrary input: an island may contain a feature narrower than the minimum printable bead, and consecutive layers may overlap by an arbitrarily tiny area. Such models remain geometrically eligible and may be exported with prominent geometric warnings when a structurally bounded path exists. The material model is a rounded-rectangle approximation and does not simulate pressure, cooling, die swell, adhesion, or machine collisions.

Use classic Spiral vase with Continuous slicing disabled when a conventional hollow single-wall vase is required. See [Continuous Slicing Development Audit](continuous-slicing-development-audit.md) for implementation evidence and remaining risks.
