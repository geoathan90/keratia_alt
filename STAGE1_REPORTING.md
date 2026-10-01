# Stage 1 reporting cleanup

This stage adds grouped reporting without changing the stiffness solver.

## New file

```text
reporting.py
```

It provides reporting-only helpers that operate on the existing `results` dictionary returned by `solve_model(...)`.

## Main concepts

### `region` and `subregion`

Elements are grouped into engineering regions such as:

```text
horn_xpos
horn_xneg
top_square
legs_and_faces
bottom_square
other
```

This replaces the old practice of relying mainly on an alternating all-elements table.

### `physical_member_id`

This is the continuous-for-reporting key.

The FE model can still split a real steel member into sub-elements where trusses connect. For output review, the split pieces can be aggregated by `physical_member_id` so you can inspect the maximum demand along the real physical member.

## Companion notebook

```text
notebooks/03_stage1_grouped_reporting.ipynb
```

This notebook solves the current TE5/Z5 models and exports a grouped XLSX workbook using:

```python
from reporting import export_grouped_results_to_excel
```

## Workbook sheets

The grouped workbook includes:

```text
Summary
Region_Summary
Physical_Members
Critical_Ranked
Compression_Sorted
Frame_Bending_Sorted
All_Elements_Grouped
R_horn_xneg
R_horn_xpos
R_top_square
R_legs_and_faces
R_bottom_square
```

The old matched/alternating style can still be generated from the solver, but it is no longer the preferred human-facing format.
