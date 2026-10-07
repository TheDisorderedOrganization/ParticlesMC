# Single-molecule aligning-field check

This is a minimal validation of the aligning term

\[
U_h=-\mathbf h\cdot\mathbf n
\]

for the three-site ortho-terphenyl model. The system contains one isolated
molecule and uses displacement moves only. Its orientation is the ordered plane
normal from `PlaneNormalOrientation(1, 2, 3)`. The field points along `+z`, while
the initial molecular normal is perpendicular to it.

Run the scan from the repository root:

```bash
julia --project=. examples/ortho-terphenyl/5-single-molecule-field-test/run-field-scan.jl
```

The script runs independent simulations for
`|h| = 0, 0.5, 1, 3, 10` at `T = 1`, samples after burn-in, and writes
`alignment-vs-field.csv`. The measured quantity is

\[
\langle \hat{\mathbf h}\cdot\mathbf n\rangle.
\]

It should be close to zero without a field and should approach one as the field
magnitude increases. Finite sampling means individual estimates need not be
strictly monotonic. This test intentionally uses a one-molecule system; applying
a field to a selected molecule inside the full liquid would require a separate
molecule-selection option.

One run with the default seed produced:

| field magnitude | mean alignment | standard error |
|---:|---:|---:|
| 0.0 | 0.0310 | 0.0129 |
| 0.5 | 0.3134 | 0.0117 |
| 1.0 | 0.3415 | 0.0113 |
| 3.0 | 0.6584 | 0.0075 |
| 10.0 | 0.8995 | 0.0022 |

The detailed values are stored in `alignment-vs-field.csv`.

An XYZ trajectory is also saved for each field magnitude at:

```text
output/h_FIELD/chains/1/trajectory.xyz
```

For example, the strongest-field trajectory is
`output/h_10.0/chains/1/trajectory.xyz`. Frames are written every 100 MC steps,
including the initial configuration, and can be opened with OVITO, VMD, or
another XYZ-compatible molecular viewer.
