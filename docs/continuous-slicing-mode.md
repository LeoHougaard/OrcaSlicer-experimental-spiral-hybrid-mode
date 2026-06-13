# Experimental Continuous Slicing Mode

This documentation describes a fork-specific feature, not standard upstream OrcaSlicer behavior.

This fork adds a new **continuous extrusion slicing mode** on top of OrcaSlicer's existing **Spiral vase** setting. The mode is still highly experimental. It is intended for parts where classic vase mode is too restrictive, but where you still want slicer behavior that avoids disruptive wall crossings as much as possible.

Do not assume this mode produces true vase-mode output or production-ready continuous extrusion G-code. Always inspect the generated toolpath in preview before printing.

In the UI, enable **Spiral vase** first, then enable **Spiral hybrid non-crossing** in the advanced print settings.

## What It Changes

Classic Spiral vase mode forces the print toward a hollow, single-wall object:

- `wall_loops` is forced to `1`.
- `top_shell_layers` is forced to `0`.
- `sparse_infill_density` is forced to `0`.
- model slicing uses the largest positive contour.
- the classic SpiralVase G-code post-processor handles spiralized extrusion.

With **Spiral hybrid non-crossing** enabled, those classic restrictions are relaxed:

- multiple walls are allowed.
- top shell layers are allowed.
- sparse infill is allowed.
- full model contours are sliced instead of only the largest positive contour.
- OrcaSlicer enables `reduce_crossing_wall` to avoid wall crossings more aggressively.
- layer-change retractions are still disabled while Spiral vase is enabled.
- the classic SpiralVase G-code post-processor is bypassed, so the output should be checked in preview before treating it as equivalent to a true vase-mode spiral.

This makes the mode useful for experimental non-crossing prints that need more structure than a hollow vase.

## Settings

### Spiral hybrid non-crossing

Config key: `spiral_hybrid_non_crossing`

Enables the hybrid behavior. This setting only applies when `spiral_mode` is also enabled.

When disabled, Spiral vase keeps the standard OrcaSlicer behavior and continues to force a single wall, no top shell, and no sparse infill.

### Spiral hybrid interior clearance

Config key: `spiral_hybrid_interior_clearance`

Default: `0.2 mm`

Controls how much interior clearance is reserved before infill is generated. When this value is greater than zero, the normal infill-to-wall overlap settings are forced to zero and `filter_out_gap_fill` is raised to at least this clearance value.

Increase this value if internal paths are too close to the outer wall. Decrease it if the print loses too much internal material.

### Spiral hybrid flow mode

Config key: `spiral_hybrid_flow_mode`

Values:

- `adaptive`: keeps the existing Spiral vase start and finish flow transition controls available.
- `constant`: hides those start and finish flow transition controls for the hybrid path.

Use `adaptive` when you want to keep tuning the Spiral vase start and finish flow ratios. Use `constant` when you want to treat the hybrid path as a uniform-flow experiment.

## Recommended Starting Point

Start with:

- `spiral_mode = 1`
- `spiral_hybrid_non_crossing = 1`
- `spiral_hybrid_interior_clearance = 0.2`
- `spiral_hybrid_flow_mode = adaptive`
- `reduce_crossing_wall = 1`

Then tune walls, top shells, and sparse infill as needed for the part.

## Known Limitations

This mode is highly experimental and fork-specific. Preview the generated toolpath before printing, and expect behavior to change as the implementation evolves.

- It is not the same as classic single-wall vase mode.
- The hybrid path currently bypasses the classic SpiralVase G-code post-processor.
- It permits walls and infill, so it can still contain internal path transitions.
- Multi-material prints and complex support-heavy parts need extra preview checks.
- If you need a true hollow vase print, leave **Spiral hybrid non-crossing** disabled and use classic **Spiral vase**.

## CLI/Profile Example

```ini
spiral_mode = 1
spiral_hybrid_non_crossing = 1
spiral_hybrid_interior_clearance = 0.2
spiral_hybrid_flow_mode = adaptive
wall_loops = 2
sparse_infill_density = 10%
top_shell_layers = 2
```
