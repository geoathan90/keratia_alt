"""
reporting.py
============

Grouped reporting helpers for the keratia tower model.

The functions in this module do not alter the analysis model. They enrich result
tables with reporting metadata and export grouped XLSX workbooks that are easier
to inspect than one large all-elements table.
"""

from __future__ import annotations

from pathlib import Path
from typing import Mapping, Optional, Sequence

import numpy as np
import pandas as pd

REPORTING_REGIONS = [
    "horn_xneg",
    "horn_xpos",
    "top_square",
    "legs_and_faces",
    "bottom_square",
    "other",
]

DEFAULT_GROUPED_DISPLAY_COLS = [
    "model",
    "region",
    "subregion",
    "physical_member_id",
    "element_id",
    "element_type",
    "start_label",
    "end_label",
    "member_role",
    "section_name",
    "axial_state",
    "sigma_extreme_signed_MPa",
    "sigma_max_abs_MPa",
    "sigma_axial_signed_MPa",
    "axial_force_signed_N",
    "max_abs_N_N",
    "max_M_resultant_Nmm",
    "stress_basis",
]


def _x_side_from_labels(start_label: str, end_label: str, fallback: str = "") -> str:
    """Infer x-positive/x-negative side from canonical labels."""
    txt = f"{start_label} {end_label} {fallback}".lower()
    if "xpos" in txt:
        return "xpos"
    if "xneg" in txt:
        return "xneg"
    return "unknown"


def infer_region_subregion(row: pd.Series) -> tuple[str, str]:
    """Infer a human-facing output group for one element row."""
    role = str(row.get("member_role", "")).lower()
    etype = str(row.get("element_type", "frame")).lower()
    start = str(row.get("start_label", ""))
    end = str(row.get("end_label", ""))
    labels = f"{start} {end}".lower()
    side = _x_side_from_labels(start, end, role)

    if "bottom_square" in role:
        return "bottom_square", "bottom_frame" if etype == "frame" else "bottom_truss"

    if "top_square" in role:
        if "diagonal" in role:
            return "top_square", "planar_diagonal_truss"
        return "top_square", "square_frame" if etype == "frame" else "top_square_truss"

    if (
        "high_xpos" in labels
        or "high_xneg" in labels
        or "load_point" in role
        or "upper_truss" in role
        or "inserted_upper_truss" in role
        or "inserted_diagonal_truss" in role
        or "inserted_center_beam_truss" in role
    ):
        region = "horn_xpos" if side == "xpos" else "horn_xneg" if side == "xneg" else "other"
        if etype == "frame":
            if "center" in role:
                return region, "horn_frame_center_to_load"
            return region, "horn_frame_load_to_top"
        if "center" in role:
            return region, "horn_center_truss"
        if "diagonal" in role:
            return region, "horn_diagonal_truss"
        return region, "horn_truss"

    if "leg" in role or "face" in role or "inscribed" in role or "node_to" in role:
        if etype == "frame":
            return "legs_and_faces", "leg_frame"
        if "xpos" in role:
            return "legs_and_faces", "face_truss_xpos"
        if "xneg" in role:
            return "legs_and_faces", "face_truss_xneg"
        if "ypos" in role:
            return "legs_and_faces", "face_truss_ypos"
        if "yneg" in role:
            return "legs_and_faces", "face_truss_yneg"
        return "legs_and_faces", "face_truss"

    if "center_to_top" in role:
        return "legs_and_faces", "center_to_top_frame"

    return "other", f"{etype}_unclassified"


def infer_physical_member_id(row: pd.Series) -> str:
    """Infer the reporting key that groups split FE elements into one steel piece."""
    etype = str(row.get("element_type", "frame")).lower()
    parent = row.get("parent_element_id", np.nan)
    prefix = "F" if etype == "frame" else "T"
    if pd.notna(parent):
        try:
            return f"{prefix}_parent_{int(float(parent)):03d}"
        except Exception:
            pass
    try:
        return f"{prefix}_element_{int(row.get('element_id')):03d}"
    except Exception:
        return f"{prefix}_unknown"


def add_reporting_metadata(elements: pd.DataFrame) -> pd.DataFrame:
    """Ensure an element table contains region, subregion, and physical_member_id."""
    out = elements.copy()
    for col in ["region", "subregion", "physical_member_id"]:
        if col not in out.columns:
            out[col] = ""

    inferred = out.apply(infer_region_subregion, axis=1, result_type="expand")
    inferred.columns = ["_region", "_subregion"]
    out = pd.concat([out, inferred], axis=1)

    for col, inferred_col in [("region", "_region"), ("subregion", "_subregion")]:
        blank = out[col].isna() | out[col].astype(str).str.strip().eq("")
        out.loc[blank, col] = out.loc[blank, inferred_col]

    blank = out["physical_member_id"].isna() | out["physical_member_id"].astype(str).str.strip().eq("")
    out.loc[blank, "physical_member_id"] = out.loc[blank].apply(infer_physical_member_id, axis=1)
    return out.drop(columns=["_region", "_subregion"])


def enrich_result_for_reporting(result: Mapping[str, object]) -> dict:
    """Merge reporting metadata from the element table into member_forces."""
    out = dict(result)
    elements = add_reporting_metadata(pd.DataFrame(out["elements"]).copy())
    mf = pd.DataFrame(out["member_forces"]).copy()
    meta_cols = ["element_id", "region", "subregion", "physical_member_id"]
    mf = mf.drop(columns=[c for c in meta_cols if c != "element_id" and c in mf.columns], errors="ignore")
    mf = mf.merge(elements[meta_cols], on="element_id", how="left")
    out["elements"] = elements
    out["member_forces"] = mf
    return out


def enrich_results_for_reporting(results: Mapping[str, Mapping[str, object]]) -> dict[str, dict]:
    """Apply enrich_result_for_reporting to every result in a dict."""
    return {name: enrich_result_for_reporting(res) for name, res in results.items()}


def _region_sort(df: pd.DataFrame) -> pd.DataFrame:
    """
    Sort a reporting table by engineering region when possible.

    Some summary tables, such as Region_Summary, do not have element_id or
    physical_member_id columns. Earlier versions assumed element_id always
    existed, which caused KeyError in grouped reporting exports. This sorter only
    uses columns present in the specific table being sorted.
    """
    if df.empty:
        return df.copy()

    out = df.copy()
    order = {name: i for i, name in enumerate(REPORTING_REGIONS)}
    if "region" in out.columns:
        out["_region_order"] = out["region"].map(order).fillna(999).astype(int)
    else:
        out["_region_order"] = 999

    sort_candidates = [
        "model",
        "_region_order",
        "region",
        "subregion",
        "physical_member_id",
        "element_id",
        "critical_element_id",
    ]
    sort_cols = [c for c in sort_candidates if c in out.columns]
    return out.sort_values(sort_cols).drop(columns=["_region_order"]).reset_index(drop=True)


def all_elements_grouped_table(results: Mapping[str, Mapping[str, object]]) -> pd.DataFrame:
    """All element result rows sorted by region/subregion."""
    results = enrich_results_for_reporting(results)
    tables = []
    for name, res in results.items():
        table = pd.DataFrame(res["member_forces"]).copy()
        table.insert(0, "model", name)
        tables.append(table)
    if not tables:
        return pd.DataFrame()
    return _region_sort(pd.concat(tables, ignore_index=True))


def region_summary_table(results: Mapping[str, Mapping[str, object]]) -> pd.DataFrame:
    """Peak demand summary by model/region/subregion."""
    results = enrich_results_for_reporting(results)
    rows = []
    for name, res in results.items():
        mf = pd.DataFrame(res["member_forces"]).copy()
        for (region, subregion), g in mf.groupby(["region", "subregion"], dropna=False):
            if g.empty:
                continue
            crit = g.loc[g["sigma_max_abs_MPa"].idxmax()]
            comp = g.loc[g["axial_force_signed_N"].idxmin()]
            tens = g.loc[g["axial_force_signed_N"].idxmax()]
            rows.append({
                "model": name,
                "region": region,
                "subregion": subregion,
                "element_count": int(len(g)),
                "frame_count": int((g["element_type"] == "frame").sum()),
                "truss_count": int((g["element_type"] == "truss").sum()),
                "max_sigma_MPa": float(g["sigma_max_abs_MPa"].max()),
                "critical_element_id": int(crit["element_id"]),
                "critical_physical_member_id": crit.get("physical_member_id", ""),
                "max_abs_axial_force_kN": float(g["max_abs_N_N"].max() / 1000.0),
                "most_compressive_force_kN": float(comp["axial_force_signed_N"] / 1000.0),
                "most_compressive_element_id": int(comp["element_id"]),
                "largest_tensile_force_kN": float(tens["axial_force_signed_N"] / 1000.0),
                "largest_tensile_element_id": int(tens["element_id"]),
                "max_bending_resultant_kNmm": float(g["max_M_resultant_Nmm"].max() / 1000.0),
            })
    return _region_sort(pd.DataFrame(rows)) if rows else pd.DataFrame()


def physical_member_summary_table(results: Mapping[str, Mapping[str, object]]) -> pd.DataFrame:
    """Aggregate split FE elements into continuous physical members for reporting."""
    results = enrich_results_for_reporting(results)
    rows = []
    for name, res in results.items():
        mf = pd.DataFrame(res["member_forces"]).copy()
        for physical_member_id, g in mf.groupby("physical_member_id", dropna=False):
            stress_row = g.loc[g["sigma_max_abs_MPa"].idxmax()]
            moment_row = g.loc[g["max_M_resultant_Nmm"].idxmax()]
            comp_row = g.loc[g["axial_force_signed_N"].idxmin()]
            tens_row = g.loc[g["axial_force_signed_N"].idxmax()]
            rows.append({
                "model": name,
                "region": g["region"].mode().iloc[0] if not g["region"].mode().empty else "other",
                "subregion": g["subregion"].mode().iloc[0] if not g["subregion"].mode().empty else "other",
                "physical_member_id": physical_member_id,
                "element_type": ",".join(sorted(set(g["element_type"].astype(str)))),
                "member_role": g["member_role"].mode().iloc[0] if not g["member_role"].mode().empty else "",
                "section_name": g["section_name"].mode().iloc[0] if not g["section_name"].mode().empty else "",
                "sub_element_count": int(len(g)),
                "sub_element_ids": ",".join(str(int(x)) for x in sorted(g["element_id"].astype(int))),
                "total_sub_element_length_mm": float(g["length_mm"].sum()),
                "max_sigma_MPa": float(g["sigma_max_abs_MPa"].max()),
                "critical_element_id": int(stress_row["element_id"]),
                "max_bending_resultant_kNmm": float(g["max_M_resultant_Nmm"].max() / 1000.0),
                "max_bending_element_id": int(moment_row["element_id"]),
                "most_compressive_force_kN": float(comp_row["axial_force_signed_N"] / 1000.0),
                "most_compressive_element_id": int(comp_row["element_id"]),
                "largest_tensile_force_kN": float(tens_row["axial_force_signed_N"] / 1000.0),
                "largest_tensile_element_id": int(tens_row["element_id"]),
            })
    return _region_sort(pd.DataFrame(rows)) if rows else pd.DataFrame()


def compression_member_table(results: Mapping[str, Mapping[str, object]], n: Optional[int] = None) -> pd.DataFrame:
    """Elements sorted from most compressive axial force upward."""
    table = all_elements_grouped_table(results)
    if table.empty:
        return table
    out = table.sort_values(["axial_force_signed_N", "model", "element_id"]).reset_index(drop=True)
    return out if n is None else out.head(n).copy()


def frame_bending_table(results: Mapping[str, Mapping[str, object]], n: Optional[int] = None) -> pd.DataFrame:
    """Frame elements sorted by bending resultant demand."""
    table = all_elements_grouped_table(results)
    if table.empty:
        return table
    out = table.loc[table["element_type"].eq("frame")].copy()
    out = out.sort_values("max_M_resultant_Nmm", ascending=False).reset_index(drop=True)
    return out if n is None else out.head(n).copy()


def export_grouped_results_to_excel(
    results: Mapping[str, Mapping[str, object]],
    out_xlsx: str | Path,
    top_n: int = 20,
    display_columns: Optional[Sequence[str]] = None,
) -> None:
    """Export a grouped XLSX report."""
    results = enrich_results_for_reporting(results)
    out_xlsx = Path(out_xlsx)
    out_xlsx.parent.mkdir(parents=True, exist_ok=True)
    display_columns = list(display_columns or DEFAULT_GROUPED_DISPLAY_COLS)

    def cols_existing(df: pd.DataFrame, cols: Sequence[str]) -> list[str]:
        return [c for c in cols if c in df.columns]

    summary_rows = []
    for name, res in results.items():
        disp = pd.DataFrame(res["displacements"])
        mf = pd.DataFrame(res["member_forces"])
        elems = pd.DataFrame(res["elements"])
        summary_rows.append({
            "model": name,
            "nodes": len(res["nodes"]),
            "elements": len(elems),
            "frame_elements": int((elems["element_type"] == "frame").sum()) if "element_type" in elems else 0,
            "truss_elements": int((elems["element_type"] == "truss").sum()) if "element_type" in elems else 0,
            "max_translation_mm": float(disp["translation_mag_mm"].max()),
            "max_abs_Uz_mm": float(disp["Uz_mm"].abs().max()),
            "max_abs_axial_force_kN": float(mf["max_abs_N_N"].max() / 1000.0),
            "max_bending_resultant_kNmm": float(mf["max_M_resultant_Nmm"].max() / 1000.0),
            "max_sigma_MPa": float(mf["sigma_max_abs_MPa"].max()),
        })

    summary = pd.DataFrame(summary_rows)
    all_grouped = all_elements_grouped_table(results)
    region_summary = region_summary_table(results)
    physical_members = physical_member_summary_table(results)
    compression = compression_member_table(results)
    frame_bending = frame_bending_table(results)

    critical_ranked = []
    for name, res in results.items():
        mf = pd.DataFrame(res["member_forces"]).copy()
        mf.insert(0, "model", name)
        critical_ranked.append(mf.sort_values("sigma_max_abs_MPa", ascending=False).head(top_n))
    critical_ranked = pd.concat(critical_ranked, ignore_index=True) if critical_ranked else pd.DataFrame()

    with pd.ExcelWriter(out_xlsx, engine="openpyxl") as writer:
        summary.to_excel(writer, sheet_name="Summary", index=False)
        region_summary.to_excel(writer, sheet_name="Region_Summary", index=False)
        physical_members.to_excel(writer, sheet_name="Physical_Members", index=False)
        critical_ranked[cols_existing(critical_ranked, display_columns)].to_excel(writer, sheet_name="Critical_Ranked", index=False)
        compression[cols_existing(compression, display_columns)].to_excel(writer, sheet_name="Compression_Sorted", index=False)
        frame_bending[cols_existing(frame_bending, display_columns)].to_excel(writer, sheet_name="Frame_Bending_Sorted", index=False)
        all_grouped[cols_existing(all_grouped, display_columns)].to_excel(writer, sheet_name="All_Elements_Grouped", index=False)

        for region in REPORTING_REGIONS:
            part = all_grouped.loc[all_grouped["region"].eq(region)].copy() if "region" in all_grouped else pd.DataFrame()
            if not part.empty:
                part[cols_existing(part, display_columns)].to_excel(writer, sheet_name=f"R_{region}"[:31], index=False)

        for name, res in results.items():
            safe = str(name)[:20]
            for key, sheet_suffix in [
                ("displacements", "Displacements"),
                ("reactions", "Reactions"),
                ("loads_used", "Loads"),
                ("validation", "Validation"),
                ("stabilization", "Stabilization"),
                ("crossing_mode_diagnostics", "Crossing_Mode"),
                ("section_properties", "Sections"),
            ]:
                value = res.get(key)
                if isinstance(value, pd.DataFrame):
                    value.to_excel(writer, sheet_name=f"{safe}_{sheet_suffix}"[:31], index=False)

        wb = writer.book
        for ws in wb.worksheets:
            ws.freeze_panes = "A2"
            if ws.max_row >= 1 and ws.max_column >= 1:
                ws.auto_filter.ref = ws.dimensions
            for col_cells in ws.columns:
                letter = col_cells[0].column_letter
                max_len = 0
                for cell in col_cells[: min(len(col_cells), 80)]:
                    if cell.value is not None:
                        max_len = max(max_len, len(str(cell.value)))
                ws.column_dimensions[letter].width = min(max(max_len + 2, 10), 38)
