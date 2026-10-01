# Lattice mixed frame/truss explorer

This is an exploratory stiffness-method project for simplified upper portions of two steel lattice tower types, **TE5** and **Z5**.

It is intended for comparison, debugging, and intuition-building only. It is not a validated design tool and should not be used for actionable structural decisions without independent verification.

## What is modelled

- Main tower members are modelled as **3D frame elements**.
- Added secondary lattice/bracing members are modelled as **3D truss elements**.
- All nodes use a 6-DOF convention internally: `Ux, Uy, Uz, Rx, Ry, Rz`.
- Frame elements use all six DOFs at each end.
- Truss elements use only translational DOFs.
- Bottom/minimum-Z nodes are fixed supports.
- Truss-only node rotations are automatically restrained because they are unused by axial-only truss elements.

## Important crossing-diagonal option

The solver now includes `CROSSING_DIAGONAL_MODE`, with three options:

| Mode | Meaning | Notes |
|---|---|---|
| `continuous` | Crossing diagonals pass over/under each other. Split crossing nodes are removed by merging the two split truss segments back into a continuous bar. | Recommended default unless the real detail has a physical connection at the crossing. |
| `connected_stabilized` | Crossing diagonals share the crossing node. Tiny stabilizing springs are added only along missing truss-only mechanism directions. | Numerically stable and useful if the crossing is physically connected. |
| `connected_unstabilized` | Crossing diagonals share the crossing node, with no stabilization. | Diagnostic mode. Not guaranteed singular for every possible structure, but currently singular for these TE5/Z5 geometries. |

For the current default vertical load case, `continuous` and `connected_stabilized` give effectively identical global results to the shown precision. `connected_unstabilized` is included to make the matrix-rank issue visible.

## Project structure

```text
lattice_mixed_frame_truss_explorer_v2/
├─ solver.py
├─ helpers.py
├─ README.md
├─ requirements.txt
├─ data/
│  ├─ TE5_nodes_aligned.csv
│  ├─ TE5_elements_with_sections.csv
│  ├─ Z5_nodes_aligned.csv
│  ├─ Z5_elements_with_sections.csv
│  └─ default_high_point_loads_by_label.csv
├─ notebooks/
│  ├─ 01_mixed_frame_truss_workflow.ipynb
│  └─ 02_solver_walkthrough.ipynb
└─ results/
   ├─ crossing_mode_comparison.csv
   ├─ default_continuous.xlsx
   └─ default_connected_stabilized.xlsx
```

## Install in an online-only GitHub Codespace

From the project root:

```bash
python -m pip install --user -r requirements.txt
python -m ipykernel install --user --name lattice-mixed-frame-truss --display-name "Python (lattice-mixed-frame-truss)"
```

Then open:

```text
notebooks/01_mixed_frame_truss_workflow.ipynb
```

and select the kernel:

```text
Python (lattice-mixed-frame-truss)
```

No virtual environment is required.

## Main notebook controls

In `01_mixed_frame_truss_workflow.ipynb`, edit:

```python
SECTION_MODE = "principal"  # or "isotropic"
CROSSING_DIAGONAL_MODE = "continuous"  # or "connected_stabilized", "connected_unstabilized"
```

The load case is controlled by:

```python
LOAD_CASE_BY_LABEL = [
    {"node_label": "high_xneg", "Fx": 0.0, "Fy": 0.0, "Fz": -13400.0, "Mx": 0.0, "My": 0.0, "Mz": 0.0},
    {"node_label": "high_xpos", "Fx": 0.0, "Fy": 0.0, "Fz": -13400.0, "Mx": 0.0, "My": 0.0, "Mz": 0.0},
]
```

DOF/component convention:

- `Fx`: global X force, N
- `Fy`: global Y force, N
- `Fz`: global Z force, N
- `Mx`: moment about global X, N·mm
- `My`: moment about global Y, N·mm
- `Mz`: moment about global Z, N·mm

## Result export

The workflow notebook calls:

```python
export_results_to_excel(results, out_xlsx, top_n=TOP_N_CRITICAL, display_columns=DISPLAY_COLS)
```

The generated XLSX workbook includes:

- `Summary`
- `Critical_Ranked`
- `Critical_Matched`
- `All_Elements_Alternating`
- per-model displacement, reaction, load, validation, stabilization, crossing-mode, and section-property sheets

`Critical_Matched` is designed for side-by-side TE5/Z5 comparison: it takes the union of critical element IDs and displays matching model rows in alternating order where possible.

## Signed stress visualization

The notebook includes both:

- absolute stress visualization: `sigma_max_abs_MPa`
- signed tension/compression visualization: `sigma_extreme_signed_MPa`

Positive signed stress is tensile; negative signed stress is compressive.

For trusses, signed stress is the axial stress. For frames, `sigma_extreme_signed_MPa` is the signed controlling extreme-fiber normal stress from axial force plus bending. Use `sigma_axial_signed_MPa` instead if you want a pure axial tension/compression load-path plot.

## Current default result summary

For the default vertical high-point load case:

```text
Fz = -13.4 kN at high_xneg
Fz = -13.4 kN at high_xpos
```

`continuous` and `connected_stabilized` solve successfully and give the same global metrics to the shown precision. `connected_unstabilized` fails for both models because the reduced stiffness matrix is rank deficient at truss-only crossing mechanisms.
