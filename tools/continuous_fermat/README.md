# Continuous Fermat Prototype and G-code Auditor

This directory contains a legacy standalone geometry prototype and an independent serialized-G-code auditor. Production planning and mandatory geometry validation are implemented in C++; prototype results never substitute for production acceptance.

Printer-ready Continuous slicing export is available from OrcaSlicer after the job passes its mandatory validators. Passing either tool in this directory is additional development evidence, not physical-printer certification.

## Geometry prototype

The dependency-free prototype samples a signed-distance field, extracts offset contours with marching squares, constructs a cyclic contour weave, and reports sampled containment/crossing/spacing/coverage metrics.

```powershell
python tools/continuous_fermat/fermat_layer.py --all --draw-contours
python tools/continuous_fermat/fermat_layer.py --shape annulus --draw-contours
```

The C++ planner uses Clipper geometry, exact swept-bead area, adaptive per-segment widths, and a stricter acceptance contract. Results are intentionally not interchangeable.

## Serialized G-code audit

For a complete exported artifact, provide exact expected layer and section counts:

```powershell
python tools/continuous_fermat/validate_gcode.py path\to\model.gcode `
    --expected-layers 120 --expected-sections 120 --firmware marlin2
```

Use `--firmware klipper` for Klipper output. `--allow-unknown-counts` is diagnostic-only.

The auditor requires firmware-correct absolute XYZ, relative E, unity speed/flow overrides, disabled pressure advance, all-axis homing, and heater-ready state. It checks:

- balanced standalone section markers and exact layer/section counts;
- finite coordinates, positive feed, and positive constant Z in each section;
- strictly positive E on every in-section XY move;
- no travel, retract, Z, arc, reset, or disruptive modal change inside a section;
- closed section endpoints;
- no unmarked model extrusion;
- between sections, one positive Z move followed by at most one non-extruding linear XY approach;
- after the final section, at most one bounded positive Z lift and only allowlisted shutdown commands.

The optional `--allow-unmarked-before-first-section` switch exists for reviewing historical files with a purge before the first marked section. It never permits unmarked extrusion between or after sections.

The auditor does not recompute variable-width swept geometry, quantized-E material, volumetric limits, collision clearance, leveling behavior, or real bead formation. Those geometry/material checks remain in C++, while machine-specific physical qualification remains outside this software audit.

Run its unit tests with:

```powershell
python -m unittest tools.continuous_fermat.test_validate_gcode
```
