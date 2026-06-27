# Continuous Path Lab Iteration Loop

Goal:

- Strict coverage >= 0.990 on every built-in shape.
- Spacing warnings <= 3 on every built-in shape.
- Crossings and containment violations must be zero.
- Fix root causes, not per-shape symptoms.

Current loop:

1. Run the benchmark:

   ```powershell
   python tools\continuous_path_lab\benchmark.py --algorithm contour_tree_v6 --grid 80 --coverage-grid 70 --json tools\continuous_path_lab\v6_benchmark.json
   ```

2. Sort failures by benchmark score.

3. Classify each failure:

   - `coverage`: strict underfill cells/groups remain.
   - `spacing`: true different-contour or connector corridor overlap remains.
   - `crossing`: route connector crosses existing path.
   - `validator`: same-contour arc was counted as non-local.

4. Implement one simple generic rule.

5. Re-run the same benchmark.

6. Keep the rule only if it improves the worst failures without hiding errors.

Current v6 rules:

- Use all printable iso-contours.
- Reserve printed hole-barrier corridors in the planning SDF.
- Use semantic spacing for same-contour arc pairs.
- Add local underfill detours from strict audit samples.

Current measured blocker:

- Local detours plateau around 95-98% coverage on hard shapes.
- Residual SDF contours still produce `0` residual inserts on these shapes.
- Remaining underfill must be filled from strict-audit underfill components, not from another closed iso-contour threshold.

Next v6 instruction:

- Add component-based residual fill from `strict_coverage_audit`.
- For each large underfill component, create one medial stitch through the component.
- Insert the stitch into the nearest valid path interval.
- Validate candidate insertion with semantic spacing, crossings, containment, and strict coverage.
- Keep only insertions that improve strict coverage without raising spacing/crossing over limits.

Do not:

- Add per-shape branches.
- Suppress spacing warnings without semantic proof.
- Use legacy raster coverage as a pass condition.
- Add separate disconnected islands.
