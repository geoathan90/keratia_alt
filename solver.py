"""
solver.py
=========

A deliberately verbose, inspectable, educational 3D stiffness-method solver for
small exploratory steel lattice-tower models containing a mixture of:

    * 3D frame/beam elements
    * 3D truss/bar elements

The code is written for the modelling workflow used in the companion notebook:
CSV node and element tables exported/generated from simplified tower geometry.
It is NOT intended to replace a validated structural-analysis package.

Units
-----
This solver assumes a consistent N-mm unit system:

    length      : mm
    force       : N
    moment      : N*mm
    stress      : N/mm^2 = MPa
    E, G        : N/mm^2
    A           : mm^2
    I, J        : mm^4

Node DOFs
---------
Every node is assigned 6 global degrees of freedom, in this fixed order:

    0: Ux  translation along global X [mm]
    1: Uy  translation along global Y [mm]
    2: Uz  translation along global Z [mm]
    3: Rx  rotation about global X [rad]
    4: Ry  rotation about global Y [rad]
    5: Rz  rotation about global Z [rad]

Frame elements use all 6 DOFs at both ends.
Truss elements use only the translational DOFs Ux, Uy, Uz at both ends.

Important point for mixed frame/truss models
--------------------------------------------
A node connected only to truss elements has no bending/torsion stiffness
connected to its rotational DOFs. Those rotations are physically irrelevant
for truss-only joints, but numerically they are free zero-stiffness DOFs. The
solver therefore automatically restrains rotations at truss-only nodes. This
is not adding translational stiffness; it merely removes unused rotational
unknowns.

Support convention
------------------
By default, all nodes with z == z_min are treated as fully fixed supports.
This is important for the current tower-top models because new midpoint nodes
can be created on bottom-square members, and those should also be fixed if the
entire bottom square is intended to represent a rigid connection to the omitted
lower tower.

Element sign convention
-----------------------
For an element's local axial force:

    positive axial_force_signed_N = tension
    negative axial_force_signed_N = compression

For a frame element, the local end-force vector produced by k_local @ u_local
has N_i and N_j values such that:

    tension     -> N_i < 0 and N_j > 0
    compression -> N_i > 0 and N_j < 0

Therefore the reported signed axial force is:

    axial_force_signed_N = 0.5 * (N_j - N_i)

For a truss element, the same convention is used directly from the axial
extension.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

import math
import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# Basic constants and DOF bookkeeping
# ---------------------------------------------------------------------------

DOF_NAMES = ("Ux", "Uy", "Uz", "Rx", "Ry", "Rz")
TRANSLATION_DOF_INDICES = (0, 1, 2)
ROTATION_DOF_INDICES = (3, 4, 5)

DEFAULT_E = 210_000.0  # N/mm^2, structural steel
DEFAULT_G = 81_000.0   # N/mm^2, approximately E/[2(1+nu)] with nu≈0.30


@dataclass(frozen=True)
class ModelData:
    """Container for one tower model read from CSV files."""

    name: str
    nodes: pd.DataFrame
    elements: pd.DataFrame


@dataclass(frozen=True)
class SolveOptions:
    """
    Options controlling the modelling assumptions.

    crossing_diagonal_mode:
        "continuous" -> crossing diagonals are assumed to pass over/under each
                        other without force transfer. Split crossing nodes are
                        removed from the analysis by merging their truss
                        segments back into continuous bars. This is often the
                        safest default unless the real tower detail shows a
                        physical connection at the crossing.

        "connected_stabilized" -> crossing diagonals share the inserted
                                  intersection node. Truss-only mechanism
                                  directions are given tiny numerical springs
                                  so the matrix remains solvable.

        "connected_unstabilized" -> crossing diagonals share the inserted
                                    intersection node, but no tiny springs are
                                    added. This is useful as a diagnostic. It
                                    is not guaranteed to be singular in every
                                    possible structure, but it often is for
                                    planar truss-only crossing joints.

    section_mode:
        "principal" -> use calculated principal inertias of the L section.
                       By default Iy receives Imax and Iz receives Imin.
        "isotropic" -> use Iy = Iz = (Imax + Imin)/2. This is less realistic,
                       but it removes local-axis sensitivity.

    principal_orientation:
        "strong_inward" -> local y axis uses Imax; local z axis uses Imin.
        "weak_inward"   -> local y axis uses Imin; local z axis uses Imax.
        Only meaningful when section_mode="principal".

    fix_min_z_nodes:
        If True, all nodes at the minimum Z coordinate are fully fixed.

    fix_nodes_marked_as_support:
        If True, nodes whose is_fixed_support column is truthy are fully fixed.

    auto_fix_truss_only_rotations:
        If True, rotational DOFs at truss-only nodes are restrained. This is
        normally required for a 6-DOF mixed model.

    stabilize_truss_only_translations:
        Truss-only crossing nodes can have translational mechanism directions.
        Example: two crossing diagonals lying in one plane provide no first-order
        axial stiffness normal to that plane. If True, the solver adds very tiny
        numerical springs only along the missing translational directions of
        truss-only nodes. This is a numerical regularization, not a real design
        assumption. It should have negligible effect unless loads are applied
        directly in those mechanism directions.
    """

    crossing_diagonal_mode: str = "continuous"
    section_mode: str = "principal"
    principal_orientation: str = "strong_inward"
    fix_min_z_nodes: bool = True
    fix_nodes_marked_as_support: bool = True
    auto_fix_truss_only_rotations: bool = True
    stabilize_truss_only_translations: bool = True
    truss_only_stabilization_k_N_per_mm: float = 1.0e-3
    min_z_tolerance_mm: float = 1.0e-6
    length_tolerance_mm: float = 1.0e-9
    axial_state_tolerance_N: float = 1.0


# ---------------------------------------------------------------------------
# CSV loading and normalization
# ---------------------------------------------------------------------------


def _as_bool_series(s: pd.Series) -> pd.Series:
    """Convert common CSV truthy/falsy values to bool."""
    if s.dtype == bool:
        return s.fillna(False)
    return (
        s.fillna(False)
        .astype(str)
        .str.strip()
        .str.lower()
        .isin(["true", "1", "yes", "y", "fixed", "support"])
    )


def read_model(name: str, nodes_csv: str | Path, elements_csv: str | Path) -> ModelData:
    """
    Read one tower model from node and element CSV files.

    The function is intentionally forgiving about extra columns. It keeps
    additional metadata columns intact because they are useful in the notebook
    for filtering, debugging, and visualization.
    """
    nodes = pd.read_csv(nodes_csv)
    elements = pd.read_csv(elements_csv)

    required_node_cols = {"node_id", "x", "y", "z", "node_label"}
    required_element_cols = {
        "element_id", "start_node", "end_node", "start_label", "end_label",
        "element_type", "member_role", "section_name", "b1_mm", "b2_mm", "t_mm",
    }

    missing_nodes = required_node_cols - set(nodes.columns)
    missing_elems = required_element_cols - set(elements.columns)
    if missing_nodes:
        raise ValueError(f"{name}: node CSV missing required columns: {sorted(missing_nodes)}")
    if missing_elems:
        raise ValueError(f"{name}: element CSV missing required columns: {sorted(missing_elems)}")

    # Normalize IDs and coordinates.
    nodes = nodes.copy()
    elements = elements.copy()
    nodes["node_id"] = nodes["node_id"].astype(int)
    elements["element_id"] = elements["element_id"].astype(int)
    elements["start_node"] = elements["start_node"].astype(int)
    elements["end_node"] = elements["end_node"].astype(int)
    for c in ["x", "y", "z"]:
        nodes[c] = nodes[c].astype(float)

    # Normalize element_type. Missing/blank means frame, because older input
    # files initially contained only frame elements.
    elements["element_type"] = (
        elements["element_type"].fillna("frame").astype(str).str.strip().str.lower()
    )

    # Normalize material columns. If a row does not specify E/G, use structural
    # steel defaults.
    if "E_N_per_mm2" not in elements.columns:
        elements["E_N_per_mm2"] = DEFAULT_E
    if "G_N_per_mm2" not in elements.columns:
        elements["G_N_per_mm2"] = DEFAULT_G
    elements["E_N_per_mm2"] = pd.to_numeric(elements["E_N_per_mm2"], errors="coerce").fillna(DEFAULT_E)
    elements["G_N_per_mm2"] = pd.to_numeric(elements["G_N_per_mm2"], errors="coerce").fillna(DEFAULT_G)

    # Make sure the optional support/load columns exist.
    if "is_fixed_support" not in nodes.columns:
        nodes["is_fixed_support"] = False
    else:
        nodes["is_fixed_support"] = _as_bool_series(nodes["is_fixed_support"])

    if "is_load_point" not in nodes.columns:
        nodes["is_load_point"] = False
    else:
        nodes["is_load_point"] = _as_bool_series(nodes["is_load_point"])

    # Sort for reproducibility.
    nodes = nodes.sort_values("node_id").reset_index(drop=True)
    elements = elements.sort_values("element_id").reset_index(drop=True)
    return ModelData(name=name, nodes=nodes, elements=elements)




# ---------------------------------------------------------------------------
# Crossing-diagonal modelling modes
# ---------------------------------------------------------------------------

CROSSING_DIAGONAL_MODES = {
    "continuous",
    "connected_stabilized",
    "connected_unstabilized",
}


def normalize_crossing_diagonal_mode(mode: str) -> str:
    """
    Normalize user-friendly crossing-diagonal mode names.

    Accepted canonical values:
        continuous
        connected_stabilized
        connected_unstabilized

    A few readable aliases are accepted so notebook users can type the option
    naturally without remembering the exact underscore spelling.
    """
    raw = str(mode).strip().lower().replace("-", "_").replace(" ", "_")
    aliases = {
        "connected_with_stabilization": "connected_stabilized",
        "connected_with_stabilisation": "connected_stabilized",
        "connected_stabilised": "connected_stabilized",
        "connected": "connected_stabilized",
        "split_stabilized": "connected_stabilized",
        "split_stabilised": "connected_stabilized",
        "connected_without_stabilization": "connected_unstabilized",
        "connected_without_stabilisation": "connected_unstabilized",
        "connected_unstabilised": "connected_unstabilized",
        "split_unstabilized": "connected_unstabilized",
        "split_unstabilised": "connected_unstabilized",
        "pass_through": "continuous",
        "pass_through_continuous": "continuous",
    }
    out = aliases.get(raw, raw)
    if out not in CROSSING_DIAGONAL_MODES:
        raise ValueError(
            f"Unknown crossing_diagonal_mode={mode!r}. "
            f"Expected one of {sorted(CROSSING_DIAGONAL_MODES)}."
        )
    return out


def _is_truss_intersection_node_row(row: pd.Series) -> bool:
    """Return True when a node row represents a generated truss crossing node."""
    label = str(row.get("node_label", "")).lower()
    group = str(row.get("node_group", "")).lower()
    return ("truss_intersection" in group) or ("truss_xing" in label)


def _is_true_crossing_split_element_row(row: pd.Series) -> bool:
    """
    Return True for truss rows created by splitting a real 3D truss-truss crossing.

    We intentionally use both source_note/member_role and parent_element_id.
    This avoids accidentally merging ordinary frame/truss splits that were made
    to create real connection nodes along a member.
    """
    if str(row.get("element_type", "")).lower() != "truss":
        return False
    parent = row.get("parent_element_id", np.nan)
    if pd.isna(parent):
        return False
    note = str(row.get("source_note", "")).lower()
    role = str(row.get("member_role", "")).lower()
    return ("split at true 3d truss-truss intersection" in note) or ("split_at_crossing" in role)


def merge_crossing_truss_splits_to_continuous(model: ModelData) -> tuple[ModelData, pd.DataFrame]:
    """
    Return a copy of the model where crossing diagonals are continuous bars.

    The input CSVs are usually stored in the safer explicit/split form, where
    crossing diagonals are divided at generated `truss_intersection` nodes. This
    function reverses only those split-at-crossing trusses:

        A -- X and X -- B  ->  A -------- B

    It does NOT merge ordinary splits made to connect real brace nodes or frame
    midpoint nodes. It also removes now-unused crossing nodes from the model so
    the global stiffness matrix does not contain isolated zero-stiffness DOFs.
    """
    nodes = model.nodes.copy()
    elements = model.elements.copy()

    crossing_node_ids = set(
        int(r["node_id"])
        for _, r in nodes.iterrows()
        if _is_truss_intersection_node_row(r)
    )
    if not crossing_node_ids:
        return model, pd.DataFrame()

    split_mask = elements.apply(_is_true_crossing_split_element_row, axis=1)
    split_rows = elements.loc[split_mask].copy()
    keep_rows = elements.loc[~split_mask].copy()

    node_labels = nodes.set_index("node_id")["node_label"].to_dict()
    merged_records = []
    diagnostics = []

    if split_rows.empty:
        # There are crossing-looking nodes but no recoverable split rows. Keep
        # model unchanged; report this explicitly in diagnostics.
        diagnostics.append({
            "model": model.name,
            "mode_action": "continuous_merge",
            "status": "no_split_rows_found",
            "details": f"crossing_node_count={len(crossing_node_ids)}",
        })
        return model, pd.DataFrame(diagnostics)

    # parent_element_id names the original truss element before it was split.
    # Groups normally contain exactly two segments, but the algorithm is robust
    # to more than two collinear split segments.
    for parent_id, group in split_rows.groupby("parent_element_id", sort=True):
        parent_id_int = int(float(parent_id))
        endpoints = []
        for _, row in group.iterrows():
            endpoints.extend([int(row["start_node"]), int(row["end_node"])])
        non_crossing = [n for n in endpoints if n not in crossing_node_ids]
        # Unique while preserving order.
        unique_non_crossing = []
        for n in non_crossing:
            if n not in unique_non_crossing:
                unique_non_crossing.append(n)

        if len(unique_non_crossing) != 2:
            diagnostics.append({
                "model": model.name,
                "mode_action": "continuous_merge",
                "status": "skipped_unexpected_endpoint_count",
                "parent_element_id": parent_id_int,
                "details": f"non_crossing_endpoints={unique_non_crossing}",
            })
            # If something is ambiguous, keep the original split segments rather
            # than silently deleting structural connectivity.
            keep_rows = pd.concat([keep_rows, group], ignore_index=True)
            continue

        a, b = unique_non_crossing
        template = group.sort_values(["split_order", "element_id"], na_position="last").iloc[0].copy()
        template["element_id"] = parent_id_int
        template["start_node"] = a
        template["end_node"] = b
        template["start_label"] = node_labels.get(a, str(a))
        template["end_label"] = node_labels.get(b, str(b))
        template["member_role"] = str(template.get("member_role", "")).replace("_split_at_crossing", "")
        template["source_note"] = "continuous crossing-diagonal mode: merged split truss segments"
        template["parent_element_id"] = np.nan
        template["split_order"] = np.nan
        merged_records.append(template)
        diagnostics.append({
            "model": model.name,
            "mode_action": "continuous_merge",
            "status": "merged",
            "parent_element_id": parent_id_int,
            "start_node": a,
            "end_node": b,
            "start_label": node_labels.get(a, str(a)),
            "end_label": node_labels.get(b, str(b)),
            "split_segment_count": int(len(group)),
        })

    if merged_records:
        merged_df = pd.DataFrame(merged_records)
        new_elements = pd.concat([keep_rows, merged_df], ignore_index=True, sort=False)
    else:
        new_elements = keep_rows.copy()

    # Remove nodes that no longer have any element incident to them. This removes
    # crossing nodes in continuous mode, but keeps any inserted/support nodes that
    # still belong to frame/truss elements.
    referenced = set(new_elements["start_node"].astype(int)) | set(new_elements["end_node"].astype(int))
    removed_nodes = set(nodes["node_id"].astype(int)) - referenced
    new_nodes = nodes.loc[nodes["node_id"].astype(int).isin(referenced)].copy()

    if removed_nodes:
        diagnostics.append({
            "model": model.name,
            "mode_action": "remove_unreferenced_nodes",
            "status": "removed",
            "details": ", ".join(str(n) for n in sorted(removed_nodes)),
        })

    new_elements["element_id"] = new_elements["element_id"].astype(int)
    new_elements["start_node"] = new_elements["start_node"].astype(int)
    new_elements["end_node"] = new_elements["end_node"].astype(int)
    new_nodes = new_nodes.sort_values("node_id").reset_index(drop=True)
    new_elements = new_elements.sort_values("element_id").reset_index(drop=True)

    return ModelData(model.name, new_nodes, new_elements), pd.DataFrame(diagnostics)


def prepare_model_for_crossing_diagonal_mode(model: ModelData, options: SolveOptions) -> tuple[ModelData, SolveOptions, pd.DataFrame]:
    """
    Apply the selected crossing-diagonal modelling mode.

    The returned options may differ from the input options because the crossing
    mode deliberately controls stabilization:

        continuous              -> stabilization off
        connected_stabilized    -> stabilization on
        connected_unstabilized  -> stabilization off
    """
    mode = normalize_crossing_diagonal_mode(options.crossing_diagonal_mode)

    if mode == "continuous":
        working_model, diagnostics = merge_crossing_truss_splits_to_continuous(model)
        working_options = replace(options, crossing_diagonal_mode=mode, stabilize_truss_only_translations=False)
    elif mode == "connected_stabilized":
        working_model = model
        diagnostics = pd.DataFrame([{
            "model": model.name,
            "mode_action": "crossing_mode",
            "status": "connected_stabilized",
            "details": "crossing truss nodes retained; tiny truss-only stabilization enabled",
        }])
        working_options = replace(options, crossing_diagonal_mode=mode, stabilize_truss_only_translations=True)
    elif mode == "connected_unstabilized":
        working_model = model
        diagnostics = pd.DataFrame([{
            "model": model.name,
            "mode_action": "crossing_mode",
            "status": "connected_unstabilized",
            "details": "crossing truss nodes retained; no truss-only stabilization applied",
        }])
        working_options = replace(options, crossing_diagonal_mode=mode, stabilize_truss_only_translations=False)
    else:  # defensive; normalize already validates
        raise AssertionError(mode)

    return working_model, working_options, diagnostics


# ---------------------------------------------------------------------------
# L-section properties
# ---------------------------------------------------------------------------


def _composite_rect_section_properties(b1: float, b2: float, t: float) -> Dict[str, float]:
    """
    Compute exact sharp-corner L-section area, centroid, and non-principal
    centroidal inertias using a composite-rectangle model.

    Geometry convention used only for property calculation:

        * horizontal leg: width b1 in local y direction, thickness t in z
        * vertical leg:   height b2 in local z direction, thickness t in y
        * overlap square: t x t is subtracted once

    The outside heel/vertex of the angle is at (y,z) = (0,0), with both legs
    extending into positive y and positive z. Real rolled angles have root
    radii; this idealization intentionally ignores them.
    """
    if b1 <= 0 or b2 <= 0 or t <= 0:
        raise ValueError(f"Invalid L-section dimensions: b1={b1}, b2={b2}, t={t}")
    if t >= min(b1, b2):
        raise ValueError(f"Thickness must be smaller than both legs: b1={b1}, b2={b2}, t={t}")

    # Each tuple is (sign, width_y, height_z, centroid_y, centroid_z).
    # We add the two rectangles and subtract the overlap square.
    parts = [
        (+1.0, b1, t,  b1 / 2.0, t / 2.0),
        (+1.0, t,  b2, t / 2.0,  b2 / 2.0),
        (-1.0, t,  t,  t / 2.0,  t / 2.0),
    ]

    A = 0.0
    Sy = 0.0
    Sz = 0.0
    for sign, w, h, yc, zc in parts:
        a = sign * w * h
        A += a
        Sy += a * yc
        Sz += a * zc
    ybar = Sy / A
    zbar = Sz / A

    Iy = 0.0   # about centroidal y-axis: ∫ z^2 dA
    Iz = 0.0   # about centroidal z-axis: ∫ y^2 dA
    Iyz = 0.0  # product: ∫ y z dA about centroidal axes

    for sign, w, h, yc, zc in parts:
        a_abs = w * h
        a = sign * a_abs
        # Centroidal rectangle inertias about axes parallel to y/z.
        Iy_c = w * h**3 / 12.0
        Iz_c = h * w**3 / 12.0
        dy = yc - ybar
        dz = zc - zbar
        Iy += sign * Iy_c + a * dz**2
        Iz += sign * Iz_c + a * dy**2
        Iyz += a * dy * dz

    # Principal moments. The sign of Iyz affects the angle, but the eigenvalues
    # always come from the same invariant expression.
    Iavg = 0.5 * (Iy + Iz)
    radius = math.sqrt((0.5 * (Iy - Iz)) ** 2 + Iyz**2)
    Imax = Iavg + radius
    Imin = Iavg - radius

    # Thin-walled open-section torsional constant approximation. This is the
    # St. Venant J for a thin open strip; adequate for exploratory modelling.
    J = (t**3 / 3.0) * (b1 + b2 - t)

    # Estimate extreme fiber distances in principal axes by projecting the
    # idealized L-section polygon vertices onto the principal directions.
    # The inertia-axis matrix is set up so that e.T @ M @ e gives moment of
    # inertia about an axis with unit vector e in the original y-z plane.
    M = np.array([[Iy, -Iyz], [-Iyz, Iz]], dtype=float)
    eigvals, eigvecs = np.linalg.eigh(M)  # ascending eigenvalues
    idx_min = int(np.argmin(eigvals))
    idx_max = int(np.argmax(eigvals))
    e_min = eigvecs[:, idx_min]
    e_max = eigvecs[:, idx_max]

    polygon = np.array(
        [
            [0.0, 0.0],
            [b1, 0.0],
            [b1, t],
            [t, t],
            [t, b2],
            [0.0, b2],
        ],
        dtype=float,
    )
    centered = polygon - np.array([ybar, zbar])
    # If local y is the strong principal axis, c_y corresponds to projections
    # on e_max and c_z on e_min. We also store the alternative.
    proj_max = centered @ e_max
    proj_min = centered @ e_min
    c_max = float(np.max(np.abs(proj_max)))
    c_min = float(np.max(np.abs(proj_min)))

    return {
        "A_mm2": A,
        "centroid_y_mm": ybar,
        "centroid_z_mm": zbar,
        "Iy_leg_axes_mm4": Iy,
        "Iz_leg_axes_mm4": Iz,
        "Iyz_leg_axes_mm4": Iyz,
        "Imax_mm4": Imax,
        "Imin_mm4": Imin,
        "Iavg_mm4": Iavg,
        "J_mm4": J,
        "c_strong_mm": c_max,
        "c_weak_mm": c_min,
        "c_avg_mm": max(c_max, c_min),
    }


def section_properties_for_element(row: pd.Series, options: SolveOptions) -> Dict[str, float]:
    """
    Return section properties used by the solver for one element row.

    The same section-property engine is used for frame and truss elements, but
    the truss stiffness only consumes E and A. The bending/torsion properties
    are still reported because they are helpful for tables and debugging.
    """
    b1 = float(row["b1_mm"])
    b2 = float(row["b2_mm"])
    t = float(row["t_mm"])
    props = _composite_rect_section_properties(b1, b2, t)

    if options.section_mode not in {"principal", "isotropic"}:
        raise ValueError("section_mode must be 'principal' or 'isotropic'")

    if options.section_mode == "isotropic":
        Iy = props["Iavg_mm4"]
        Iz = props["Iavg_mm4"]
        cy = props["c_avg_mm"]
        cz = props["c_avg_mm"]
        section_basis = "isotropic: Iy=Iz=Iavg"
    else:
        if options.principal_orientation == "strong_inward":
            Iy = props["Imax_mm4"]
            Iz = props["Imin_mm4"]
            cy = props["c_strong_mm"]
            cz = props["c_weak_mm"]
            section_basis = "principal: local y=Imax, local z=Imin"
        elif options.principal_orientation == "weak_inward":
            Iy = props["Imin_mm4"]
            Iz = props["Imax_mm4"]
            cy = props["c_weak_mm"]
            cz = props["c_strong_mm"]
            section_basis = "principal: local y=Imin, local z=Imax"
        else:
            raise ValueError("principal_orientation must be 'strong_inward' or 'weak_inward'")

    return {
        **props,
        "E_N_per_mm2": float(row.get("E_N_per_mm2", DEFAULT_E) or DEFAULT_E),
        "G_N_per_mm2": float(row.get("G_N_per_mm2", DEFAULT_G) or DEFAULT_G),
        "A_mm2": props["A_mm2"],
        "Iy_mm4": Iy,
        "Iz_mm4": Iz,
        "J_mm4": props["J_mm4"],
        "c_y_mm": cy,
        "c_z_mm": cz,
        "section_basis": section_basis,
    }


def build_section_table(elements: pd.DataFrame, options: SolveOptions) -> pd.DataFrame:
    """Return one row of calculated section properties per element."""
    records = []
    for _, row in elements.iterrows():
        props = section_properties_for_element(row, options)
        records.append(
            {
                "element_id": int(row["element_id"]),
                "element_type": row["element_type"],
                "section_name": row["section_name"],
                "b1_mm": float(row["b1_mm"]),
                "b2_mm": float(row["b2_mm"]),
                "t_mm": float(row["t_mm"]),
                **props,
            }
        )
    return pd.DataFrame(records)


# ---------------------------------------------------------------------------
# Coordinate systems and element stiffness matrices
# ---------------------------------------------------------------------------


def _unit(v: np.ndarray, tol: float = 1e-12) -> np.ndarray:
    """Normalize a vector, raising a clear error if it is too short."""
    n = float(np.linalg.norm(v))
    if n < tol:
        raise ValueError("Cannot normalize a near-zero vector")
    return v / n


def frame_local_axes(p_i: np.ndarray, p_j: np.ndarray) -> np.ndarray:
    """
    Build a 3x3 local-to-global axis matrix for a frame element.

    Returned rows are local unit axes expressed in global coordinates:

        R[0] = local x axis, along element i -> j
        R[1] = local y axis, transverse reference axis
        R[2] = local z axis, completing right-handed system

    The local y axis is chosen to point as closely as possible toward the tower
    centreline (global x=y=0 at the same z elevation). This approximates the
    "angles face inward" modelling convention. If the inward vector is nearly
    parallel to the member axis, a global fallback vector is used.
    """
    x_local = _unit(p_j - p_i)

    midpoint = 0.5 * (p_i + p_j)
    inward = np.array([-midpoint[0], -midpoint[1], 0.0], dtype=float)

    # Fallback if a member sits on/near the centreline, or if the inward vector
    # happens to be almost parallel to the member axis.
    if np.linalg.norm(inward) < 1e-9:
        inward = np.array([0.0, 1.0, 0.0], dtype=float)

    # Project inward onto the plane perpendicular to the member axis.
    y_raw = inward - np.dot(inward, x_local) * x_local

    if np.linalg.norm(y_raw) < 1e-9:
        # Try a global Z fallback. If the member is vertical-ish, use global Y.
        trial = np.array([0.0, 0.0, 1.0], dtype=float)
        y_raw = trial - np.dot(trial, x_local) * x_local
        if np.linalg.norm(y_raw) < 1e-9:
            trial = np.array([0.0, 1.0, 0.0], dtype=float)
            y_raw = trial - np.dot(trial, x_local) * x_local

    y_local = _unit(y_raw)
    z_local = _unit(np.cross(x_local, y_local))
    # Re-orthogonalize y to remove numerical drift.
    y_local = _unit(np.cross(z_local, x_local))
    return np.vstack([x_local, y_local, z_local])


def frame_local_stiffness(E: float, G: float, A: float, Iy: float, Iz: float, J: float, L: float) -> np.ndarray:
    """
    Standard 12x12 Euler-Bernoulli 3D frame-element local stiffness matrix.

    Local DOF order:

        [u_i, v_i, w_i, rx_i, ry_i, rz_i,
         u_j, v_j, w_j, rx_j, ry_j, rz_j]

    where:

        u  = displacement along local x
        v  = displacement along local y
        w  = displacement along local z
        rx = twist about local x
        ry = rotation about local y
        rz = rotation about local z

    Bending associations:

        E*Iz controls bending in the local x-y plane: v and rz terms.
        E*Iy controls bending in the local x-z plane: w and ry terms.
    """
    k = np.zeros((12, 12), dtype=float)

    EA_L = E * A / L
    GJ_L = G * J / L

    EIy = E * Iy
    EIz = E * Iz

    # Axial terms.
    k[0, 0] = k[6, 6] = EA_L
    k[0, 6] = k[6, 0] = -EA_L

    # Torsion terms.
    k[3, 3] = k[9, 9] = GJ_L
    k[3, 9] = k[9, 3] = -GJ_L

    # Bending about local z: v/rz DOFs -> EIz.
    c1 = 12.0 * EIz / L**3
    c2 = 6.0 * EIz / L**2
    c3 = 4.0 * EIz / L
    c4 = 2.0 * EIz / L
    dof = [1, 5, 7, 11]  # v_i, rz_i, v_j, rz_j
    sub = np.array(
        [
            [ c1,  c2, -c1,  c2],
            [ c2,  c3, -c2,  c4],
            [-c1, -c2,  c1, -c2],
            [ c2,  c4, -c2,  c3],
        ]
    )
    for a in range(4):
        for b in range(4):
            k[dof[a], dof[b]] = sub[a, b]

    # Bending about local y: w/ry DOFs -> EIy.
    c1 = 12.0 * EIy / L**3
    c2 = 6.0 * EIy / L**2
    c3 = 4.0 * EIy / L
    c4 = 2.0 * EIy / L
    dof = [2, 4, 8, 10]  # w_i, ry_i, w_j, ry_j
    sub = np.array(
        [
            [ c1, -c2, -c1, -c2],
            [-c2,  c3,  c2,  c4],
            [-c1,  c2,  c1,  c2],
            [-c2,  c4,  c2,  c3],
        ]
    )
    for a in range(4):
        for b in range(4):
            k[dof[a], dof[b]] = sub[a, b]

    return k


def transformation_12x12(R: np.ndarray) -> np.ndarray:
    """
    Build the 12x12 displacement transformation matrix for a frame element.

    R maps a 3-vector from global components to local components:

        v_local = R @ v_global

    T is block diagonal with R repeated for translations and rotations at the
    two end nodes:

        u_local_12 = T @ u_global_12
    """
    T = np.zeros((12, 12), dtype=float)
    for block_start in (0, 3, 6, 9):
        T[block_start:block_start + 3, block_start:block_start + 3] = R
    return T


def truss_global_stiffness(E: float, A: float, p_i: np.ndarray, p_j: np.ndarray) -> Tuple[np.ndarray, np.ndarray, float]:
    """
    Return a 12x12 global stiffness matrix for a 3D truss embedded in 6-DOF nodes.

    Only the translational DOFs receive nonzero stiffness. Rotational rows and
    columns stay zero. This lets frame and truss elements share the same global
    DOF numbering without special global assembly logic.
    """
    x = p_j - p_i
    L = float(np.linalg.norm(x))
    n = x / L
    k3 = (E * A / L) * np.outer(n, n)
    k = np.zeros((12, 12), dtype=float)
    # translation DOFs at node i: 0,1,2; at node j: 6,7,8
    k[0:3, 0:3] = k3
    k[0:3, 6:9] = -k3
    k[6:9, 0:3] = -k3
    k[6:9, 6:9] = k3
    return k, n, L


# ---------------------------------------------------------------------------
# Global assembly, boundary conditions, and solving
# ---------------------------------------------------------------------------


def node_id_to_row_index(nodes: pd.DataFrame) -> Dict[int, int]:
    """Map node_id values to zero-based row indices used in global DOF numbering."""
    return {int(nid): i for i, nid in enumerate(nodes["node_id"].astype(int).tolist())}


def dof_index(node_row_index: int, local_dof: int) -> int:
    """Return the global DOF index for a node row index and local DOF number 0..5."""
    return 6 * node_row_index + local_dof


def element_dof_indices(start_idx: int, end_idx: int) -> List[int]:
    """Return the 12 global DOF indices for an element connecting two node row indices."""
    return [
        dof_index(start_idx, 0), dof_index(start_idx, 1), dof_index(start_idx, 2),
        dof_index(start_idx, 3), dof_index(start_idx, 4), dof_index(start_idx, 5),
        dof_index(end_idx, 0),   dof_index(end_idx, 1),   dof_index(end_idx, 2),
        dof_index(end_idx, 3),   dof_index(end_idx, 4),   dof_index(end_idx, 5),
    ]


def validate_model(model: ModelData, options: Optional[SolveOptions] = None) -> pd.DataFrame:
    """
    Perform basic model checks and return a table of findings.

    This is intentionally not silent: when you are constructing geometry by hand,
    it is better to see more diagnostics than fewer.
    """
    nodes = model.nodes
    elements = model.elements
    records = []

    def add(check: str, status: str, details: str):
        records.append({"model": model.name, "check": check, "status": status, "details": details})

    duplicate_node_ids = nodes["node_id"].duplicated().sum()
    duplicate_labels = nodes["node_label"].duplicated().sum()
    add("duplicate node IDs", "OK" if duplicate_node_ids == 0 else "ERROR", str(int(duplicate_node_ids)))
    add("duplicate node labels", "OK" if duplicate_labels == 0 else "WARNING", str(int(duplicate_labels)))

    known_nodes = set(nodes["node_id"].astype(int))
    missing_refs = elements.loc[
        ~elements["start_node"].astype(int).isin(known_nodes) |
        ~elements["end_node"].astype(int).isin(known_nodes)
    ]
    add("missing element node references", "OK" if missing_refs.empty else "ERROR", str(len(missing_refs)))

    bad_types = sorted(set(elements["element_type"]) - {"frame", "truss"})
    add("element_type values", "OK" if not bad_types else "ERROR", ", ".join(bad_types) if bad_types else "frame/truss")

    missing_dims = elements[["b1_mm", "b2_mm", "t_mm"]].isna().any(axis=1).sum()
    add("missing section dimensions", "OK" if missing_dims == 0 else "ERROR", str(int(missing_dims)))

    coords = nodes.set_index("node_id")[["x", "y", "z"]]
    lengths = []
    for _, e in elements.iterrows():
        pi = coords.loc[int(e["start_node"])].to_numpy(float)
        pj = coords.loc[int(e["end_node"])].to_numpy(float)
        lengths.append(np.linalg.norm(pj - pi))
    zero_len = sum(L < 1e-9 for L in lengths)
    add("zero-length elements", "OK" if zero_len == 0 else "ERROR", str(int(zero_len)))

    pairs = elements.apply(lambda r: tuple(sorted((int(r["start_node"]), int(r["end_node"])))), axis=1)
    duplicates = pairs.duplicated().sum()
    add("duplicate node-pair elements", "OK" if duplicates == 0 else "WARNING", str(int(duplicates)))

    # Connected components using all elements as graph edges.
    # Kept simple to avoid adding networkx as a dependency.
    adjacency: Dict[int, set] = {int(n): set() for n in nodes["node_id"]}
    for _, e in elements.iterrows():
        a, b = int(e["start_node"]), int(e["end_node"])
        adjacency[a].add(b)
        adjacency[b].add(a)
    seen = set()
    n_components = 0
    for n in adjacency:
        if n in seen:
            continue
        n_components += 1
        stack = [n]
        seen.add(n)
        while stack:
            u = stack.pop()
            for v in adjacency[u]:
                if v not in seen:
                    seen.add(v)
                    stack.append(v)
    add("connected components", "OK" if n_components == 1 else "WARNING", str(n_components))

    if options is not None:
        restrained = determine_restrained_dofs(model, options)
        fixed_nodes = sorted({d // 6 for d in restrained if d % 6 in TRANSLATION_DOF_INDICES})
        add("restrained DOF count", "INFO", str(len(restrained)))
        add("nodes with translational restraints", "INFO", str(len(fixed_nodes)))

    return pd.DataFrame(records)


def determine_restrained_dofs(model: ModelData, options: SolveOptions) -> List[int]:
    """
    Decide which global DOFs are restrained.

    Full supports:
        * nodes at minimum Z, if options.fix_min_z_nodes is True
        * nodes marked is_fixed_support, if options.fix_nodes_marked_as_support is True

    Rotational-only restraints:
        * nodes connected only to truss elements, if requested
    """
    nodes = model.nodes
    elements = model.elements
    node_index = node_id_to_row_index(nodes)
    restrained: set[int] = set()

    support_node_ids: set[int] = set()

    if options.fix_min_z_nodes:
        zmin = float(nodes["z"].min())
        at_min_z = nodes.loc[(nodes["z"] - zmin).abs() <= options.min_z_tolerance_mm, "node_id"]
        support_node_ids.update(int(n) for n in at_min_z)

    if options.fix_nodes_marked_as_support and "is_fixed_support" in nodes.columns:
        marked = nodes.loc[_as_bool_series(nodes["is_fixed_support"]), "node_id"]
        support_node_ids.update(int(n) for n in marked)

    # Fully fix support nodes: translations + rotations.
    for node_id in support_node_ids:
        idx = node_index[node_id]
        for d in range(6):
            restrained.add(dof_index(idx, d))

    if options.auto_fix_truss_only_rotations:
        connected_types: Dict[int, set[str]] = {int(n): set() for n in nodes["node_id"]}
        for _, e in elements.iterrows():
            etype = str(e["element_type"]).lower()
            connected_types[int(e["start_node"])].add(etype)
            connected_types[int(e["end_node"])].add(etype)

        for node_id, types in connected_types.items():
            # Nodes connected only to trusses have meaningless rotations in a
            # 6-DOF model. Restraining those rotations prevents singular zero
            # rows/columns, without affecting translational stiffness.
            if types and types.issubset({"truss"}) and node_id not in support_node_ids:
                idx = node_index[node_id]
                for d in ROTATION_DOF_INDICES:
                    restrained.add(dof_index(idx, d))

    return sorted(restrained)


def loads_by_label_to_vector(nodes: pd.DataFrame, load_rows: Sequence[Mapping[str, float | str]]) -> Tuple[np.ndarray, pd.DataFrame]:
    """
    Convert a list/dict-style load case into the global force vector.

    Each load row should contain:

        node_label, Fx, Fy, Fz, Mx, My, Mz

    Missing force/moment components default to zero.
    """
    n_dof = 6 * len(nodes)
    F = np.zeros(n_dof, dtype=float)
    label_to_index = {str(label): i for i, label in enumerate(nodes["node_label"].astype(str))}
    records = []

    for row in load_rows:
        label = str(row["node_label"])
        if label not in label_to_index:
            raise ValueError(f"Load references unknown node_label: {label!r}")
        idx = label_to_index[label]
        comps = [
            float(row.get("Fx", 0.0) or 0.0),
            float(row.get("Fy", 0.0) or 0.0),
            float(row.get("Fz", 0.0) or 0.0),
            float(row.get("Mx", 0.0) or 0.0),
            float(row.get("My", 0.0) or 0.0),
            float(row.get("Mz", 0.0) or 0.0),
        ]
        for d, value in enumerate(comps):
            F[dof_index(idx, d)] += value
        rec = {"node_label": label, "node_id": int(nodes.iloc[idx]["node_id"])}
        rec.update(dict(zip(["Fx", "Fy", "Fz", "Mx", "My", "Mz"], comps)))
        records.append(rec)

    return F, pd.DataFrame(records)


def read_loads_csv(loads_csv: str | Path) -> List[Mapping[str, float | str]]:
    """Read a node-label load CSV into a list of records."""
    return pd.read_csv(loads_csv).to_dict("records")


def assemble_global_stiffness(model: ModelData, options: SolveOptions) -> Tuple[np.ndarray, pd.DataFrame, Dict[int, np.ndarray], Dict[int, np.ndarray]]:
    """
    Assemble the global stiffness matrix K.

    Returns:
        K                  full global stiffness matrix
        section_table      calculated per-element section properties
        local_stiffnesses  dict element_id -> local/global local k used for force recovery
        rotations          dict element_id -> R matrix for frame force recovery

    For truss elements, local_stiffnesses stores a 2x2 axial stiffness in a
    small convention-independent way? To keep recovery simple, truss end forces
    are computed directly from displacement, so only frame local matrices are
    stored here. The dict still records frame matrices keyed by element_id.
    """
    nodes = model.nodes
    elements = model.elements
    node_index = node_id_to_row_index(nodes)
    coords_by_id = nodes.set_index("node_id")[["x", "y", "z"]]

    n_dof = 6 * len(nodes)
    K = np.zeros((n_dof, n_dof), dtype=float)
    section_table = build_section_table(elements, options)
    sec_by_eid = section_table.set_index("element_id")

    frame_k_local: Dict[int, np.ndarray] = {}
    frame_R: Dict[int, np.ndarray] = {}

    for _, e in elements.iterrows():
        eid = int(e["element_id"])
        ni = int(e["start_node"])
        nj = int(e["end_node"])
        i_idx = node_index[ni]
        j_idx = node_index[nj]
        edofs = element_dof_indices(i_idx, j_idx)

        p_i = coords_by_id.loc[ni].to_numpy(float)
        p_j = coords_by_id.loc[nj].to_numpy(float)
        L = float(np.linalg.norm(p_j - p_i))
        if L <= options.length_tolerance_mm:
            raise ValueError(f"Element {eid} has near-zero length: {L}")

        sec = sec_by_eid.loc[eid]
        E = float(sec["E_N_per_mm2"])
        G = float(sec["G_N_per_mm2"])
        A = float(sec["A_mm2"])
        Iy = float(sec["Iy_mm4"])
        Iz = float(sec["Iz_mm4"])
        J = float(sec["J_mm4"])

        etype = str(e["element_type"]).lower()
        if etype == "frame":
            R = frame_local_axes(p_i, p_j)
            T = transformation_12x12(R)
            k_local = frame_local_stiffness(E, G, A, Iy, Iz, J, L)
            k_global = T.T @ k_local @ T
            frame_k_local[eid] = k_local
            frame_R[eid] = R
        elif etype == "truss":
            k_global, _, _ = truss_global_stiffness(E, A, p_i, p_j)
        else:
            raise ValueError(f"Element {eid} has unknown element_type: {etype!r}")

        # Add element matrix into global matrix.
        ix = np.ix_(edofs, edofs)
        K[ix] += k_global

    return K, section_table, frame_k_local, frame_R



def apply_truss_only_translation_stabilization(
    K: np.ndarray,
    model: ModelData,
    options: SolveOptions,
) -> pd.DataFrame:
    """
    Add tiny springs in missing translational directions at truss-only nodes.

    Why this exists
    ---------------
    A truss element has axial stiffness only. At a node connected only to truss
    elements, the translational stiffness rank equals the number of independent
    element directions meeting there. If the connected trusses lie in a plane,
    the out-of-plane translation has zero first-order axial stiffness. This is a
    legitimate truss mechanism, but it produces a singular linear stiffness
    matrix even when no external load is applied to that mode.

    Current tower relevance
    -----------------------
    We deliberately created nodes at crossing diagonals. Some of those nodes are
    connected only to truss members in one face. In a real tower, local plates,
    eccentricities, bolt details, and member flexural stiffness prevent a truly
    free out-of-plane infinitesimal mechanism. In this exploratory mixed model,
    the cleanest numerical treatment is to add extremely small stabilizing
    springs only in the mathematically missing directions.

    What this function does
    -----------------------
    For each truss-only node:

        1. Collect unit vectors of all incident truss elements.
        2. Build S = Σ n nᵀ.
        3. Eigenvectors with near-zero eigenvalues are mechanism directions.
        4. Add k * v vᵀ to that node's translational 3x3 block.

    The default k is 1e-3 N/mm. Compare this with ordinary member axial
    stiffnesses, which are typically many thousands of N/mm. This should be
    negligible for normal load paths but enough to remove zero pivots.
    """
    records = []
    if not options.stabilize_truss_only_translations:
        return pd.DataFrame(records)

    nodes = model.nodes
    elements = model.elements
    node_index = node_id_to_row_index(nodes)
    coords = nodes.set_index("node_id")[["x", "y", "z"]]

    connected_types: Dict[int, set[str]] = {int(n): set() for n in nodes["node_id"]}
    incident_dirs: Dict[int, List[np.ndarray]] = {int(n): [] for n in nodes["node_id"]}

    for _, e in elements.iterrows():
        etype = str(e["element_type"]).lower()
        ni = int(e["start_node"])
        nj = int(e["end_node"])
        connected_types[ni].add(etype)
        connected_types[nj].add(etype)
        if etype != "truss":
            continue
        pi = coords.loc[ni].to_numpy(float)
        pj = coords.loc[nj].to_numpy(float)
        n = _unit(pj - pi)
        incident_dirs[ni].append(n)
        incident_dirs[nj].append(-n)

    k_stab = float(options.truss_only_stabilization_k_N_per_mm)
    for node_id, types in connected_types.items():
        if not types or not types.issubset({"truss"}):
            continue
        dirs = incident_dirs[node_id]
        if not dirs:
            continue
        S = np.zeros((3, 3), dtype=float)
        for n in dirs:
            S += np.outer(n, n)
        eigvals, eigvecs = np.linalg.eigh(S)
        max_eval = max(float(np.max(eigvals)), 1.0)
        missing = []
        for i, lam in enumerate(eigvals):
            if lam < 1.0e-9 * max_eval:
                v = eigvecs[:, i]
                missing.append(v)
                idx = node_index[node_id]
                transl_dofs = [dof_index(idx, d) for d in TRANSLATION_DOF_INDICES]
                K[np.ix_(transl_dofs, transl_dofs)] += k_stab * np.outer(v, v)
        if missing:
            records.append({
                "node_id": node_id,
                "node_label": nodes.loc[nodes["node_id"] == node_id, "node_label"].iloc[0],
                "incident_truss_count": len(dirs),
                "translational_rank_from_trusses": int(np.linalg.matrix_rank(S, tol=1e-9)),
                "stabilized_missing_directions": len(missing),
                "stabilization_k_N_per_mm": k_stab,
            })

    return pd.DataFrame(records)


def solve_model(
    model: ModelData,
    load_rows: Sequence[Mapping[str, float | str]],
    options: Optional[SolveOptions] = None,
) -> Dict[str, object]:
    """
    Solve one model for one load case.

    Returns a dictionary containing the original data, displacement/reaction
    tables, member-force table, section-property table, and diagnostic metadata.
    """
    if options is None:
        options = SolveOptions()

    # Apply the selected crossing-diagonal modelling convention before any
    # stiffness assembly occurs. This allows the same explicit/split CSV files
    # to support either physically connected crossings or continuous pass-through
    # diagonals without maintaining two parallel geometry datasets.
    working_model, working_options, crossing_mode_diagnostics = prepare_model_for_crossing_diagonal_mode(model, options)

    K, section_table, frame_k_local, frame_R = assemble_global_stiffness(working_model, working_options)
    stabilization = apply_truss_only_translation_stabilization(K, working_model, working_options)
    F, loads_used = loads_by_label_to_vector(working_model.nodes, load_rows)
    restrained = determine_restrained_dofs(working_model, working_options)

    n_dof = K.shape[0]
    all_dofs = np.arange(n_dof)
    free = np.setdiff1d(all_dofs, np.array(restrained, dtype=int))

    Kff = K[np.ix_(free, free)]
    Ff = F[free]

    rank = int(np.linalg.matrix_rank(Kff, tol=1e-8))
    if rank < len(free):
        raise np.linalg.LinAlgError(
            f"Reduced stiffness matrix is singular or rank deficient: rank {rank}/{len(free)}. "
            "Check supports, disconnected geometry, truss-only rotations, and missing frame connectivity."
        )

    U = np.zeros(n_dof, dtype=float)
    U[free] = np.linalg.solve(Kff, Ff)

    # Reactions are K*u - F. At free DOFs this should be near zero; at restrained
    # DOFs it is the support reaction needed to enforce the restraint.
    R = K @ U - F

    displacements = displacement_dataframe(working_model.nodes, U)
    reactions = reaction_dataframe(working_model.nodes, R, restrained)
    member_forces = member_forces_dataframe(working_model, U, section_table, frame_k_local, frame_R, working_options)

    return {
        "model": working_model.name,
        "input_nodes": model.nodes,
        "input_elements": model.elements,
        "nodes": working_model.nodes,
        "elements": working_model.elements,
        "crossing_mode_diagnostics": crossing_mode_diagnostics,
        "K": K,
        "F": F,
        "U": U,
        "R": R,
        "restrained_dofs": restrained,
        "free_dofs": free,
        "reduced_stiffness_rank": rank,
        "reduced_stiffness_size": len(free),
        "loads_used": loads_used,
        "section_properties": section_table,
        "displacements": displacements,
        "reactions": reactions,
        "member_forces": member_forces,
        "validation": validate_model(working_model, working_options),
        "stabilization": stabilization,
        "options": working_options,
    }


# ---------------------------------------------------------------------------
# Results tables
# ---------------------------------------------------------------------------


def displacement_dataframe(nodes: pd.DataFrame, U: np.ndarray) -> pd.DataFrame:
    """Return nodal displacements/rotations as a tidy DataFrame."""
    records = []
    for i, row in nodes.reset_index(drop=True).iterrows():
        comps = U[6 * i:6 * i + 6]
        records.append(
            {
                "node_id": int(row["node_id"]),
                "node_label": row["node_label"],
                "x": row["x"], "y": row["y"], "z": row["z"],
                "Ux_mm": comps[0], "Uy_mm": comps[1], "Uz_mm": comps[2],
                "Rx_rad": comps[3], "Ry_rad": comps[4], "Rz_rad": comps[5],
                "translation_mag_mm": float(np.linalg.norm(comps[:3])),
                "rotation_mag_rad": float(np.linalg.norm(comps[3:])),
            }
        )
    return pd.DataFrame(records)


def reaction_dataframe(nodes: pd.DataFrame, R: np.ndarray, restrained_dofs: Sequence[int]) -> pd.DataFrame:
    """Return reactions only at restrained DOFs, plus one row per restrained node."""
    restrained = set(int(d) for d in restrained_dofs)
    records = []
    for i, row in nodes.reset_index(drop=True).iterrows():
        dofs = [6 * i + d for d in range(6)]
        if not any(d in restrained for d in dofs):
            continue
        comps = [R[d] if d in restrained else 0.0 for d in dofs]
        records.append(
            {
                "node_id": int(row["node_id"]),
                "node_label": row["node_label"],
                "Rxn_Fx_N": comps[0], "Rxn_Fy_N": comps[1], "Rxn_Fz_N": comps[2],
                "Rxn_Mx_Nmm": comps[3], "Rxn_My_Nmm": comps[4], "Rxn_Mz_Nmm": comps[5],
                "force_reaction_mag_N": float(np.linalg.norm(comps[:3])),
                "moment_reaction_mag_Nmm": float(np.linalg.norm(comps[3:])),
            }
        )
    return pd.DataFrame(records)


def _frame_stress_from_end_forces(N: float, My: float, Mz: float, A: float, Iy: float, Iz: float, cy: float, cz: float) -> float:
    """
    Conservative-ish frame normal stress estimate at one end of an element.

    The exact signs of bending stress at individual fibers are less important
    for the current comparison than the absolute maximum. We therefore evaluate
    the four extreme combinations of y=±cy and z=±cz and return the largest
    absolute normal stress:

        sigma = N/A + Mz*y/Iz - My*z/Iy

    where local y and z are the bending axes used by the element stiffness.
    """
    candidates = []
    for y in (-cy, cy):
        for z in (-cz, cz):
            sigma = N / A + Mz * y / Iz - My * z / Iy
            candidates.append(abs(sigma))
    return float(max(candidates))





def _frame_stress_extremes_from_end_forces(N: float, My: float, Mz: float, A: float, Iy: float, Iz: float, cy: float, cz: float) -> tuple[float, float, float]:
    """
    Return signed extreme normal stresses at one frame-element end.

    Output:
        max_tension_MPa      largest positive sigma among the checked fibers
        max_compression_MPa  most negative sigma among the checked fibers
        controlling_signed   whichever of the two has larger absolute value

    This is useful for diverging tension/compression plots. A single frame
    member under bending can have tension on one side of the angle and
    compression on the other. Therefore this signed quantity should be read as
    "the sign of the controlling extreme fiber stress", not as a statement that
    the whole frame member is globally in tension or compression.
    """
    sigmas = []
    for y in (-cy, cy):
        for z in (-cz, cz):
            sigmas.append(float(N / A + Mz * y / Iz - My * z / Iy))
    max_tension = max(sigmas)
    max_compression = min(sigmas)
    controlling = max_tension if abs(max_tension) >= abs(max_compression) else max_compression
    return float(max_tension), float(max_compression), float(controlling)


def member_forces_dataframe(
    model: ModelData,
    U: np.ndarray,
    section_table: pd.DataFrame,
    frame_k_local: Mapping[int, np.ndarray],
    frame_R: Mapping[int, np.ndarray],
    options: SolveOptions,
) -> pd.DataFrame:
    """Recover local element forces and stress measures for all elements."""
    nodes = model.nodes
    elements = model.elements
    node_index = node_id_to_row_index(nodes)
    coords = nodes.set_index("node_id")[["x", "y", "z"]]
    sec_by_eid = section_table.set_index("element_id")

    records = []
    for _, e in elements.iterrows():
        eid = int(e["element_id"])
        ni = int(e["start_node"])
        nj = int(e["end_node"])
        i_idx = node_index[ni]
        j_idx = node_index[nj]
        edofs = element_dof_indices(i_idx, j_idx)
        u_global = U[edofs]
        p_i = coords.loc[ni].to_numpy(float)
        p_j = coords.loc[nj].to_numpy(float)
        L = float(np.linalg.norm(p_j - p_i))

        sec = sec_by_eid.loc[eid]
        A = float(sec["A_mm2"])
        E = float(sec["E_N_per_mm2"])
        Iy = float(sec["Iy_mm4"])
        Iz = float(sec["Iz_mm4"])
        cy = float(sec["c_y_mm"])
        cz = float(sec["c_z_mm"])

        etype = str(e["element_type"]).lower()
        if etype == "frame":
            Rmat = frame_R[eid]
            T = transformation_12x12(Rmat)
            u_local = T @ u_global
            f_local = frame_k_local[eid] @ u_local

            N_i, Vy_i, Vz_i, T_i, My_i, Mz_i = f_local[0:6]
            N_j, Vy_j, Vz_j, T_j, My_j, Mz_j = f_local[6:12]
            axial_signed = 0.5 * (N_j - N_i)
            sigma_i = _frame_stress_from_end_forces(axial_signed, My_i, Mz_i, A, Iy, Iz, cy, cz)
            sigma_j = _frame_stress_from_end_forces(axial_signed, My_j, Mz_j, A, Iy, Iz, cy, cz)
            sigma_i_tens, sigma_i_comp, sigma_i_ctrl = _frame_stress_extremes_from_end_forces(axial_signed, My_i, Mz_i, A, Iy, Iz, cy, cz)
            sigma_j_tens, sigma_j_comp, sigma_j_ctrl = _frame_stress_extremes_from_end_forces(axial_signed, My_j, Mz_j, A, Iy, Iz, cy, cz)
            sigma_max = max(sigma_i, sigma_j)
            sigma_tension_max = max(sigma_i_tens, sigma_j_tens)
            sigma_compression_min = min(sigma_i_comp, sigma_j_comp)
            sigma_extreme_signed = sigma_tension_max if abs(sigma_tension_max) >= abs(sigma_compression_min) else sigma_compression_min
            max_M_resultant = max(math.hypot(My_i, Mz_i), math.hypot(My_j, Mz_j))
            stress_basis = "frame axial+bending"
        elif etype == "truss":
            xhat = (p_j - p_i) / L
            ui = U[6 * i_idx:6 * i_idx + 3]
            uj = U[6 * j_idx:6 * j_idx + 3]
            axial_extension = float(np.dot(uj - ui, xhat))
            axial_signed = E * A / L * axial_extension
            N_i = -axial_signed
            N_j = axial_signed
            Vy_i = Vz_i = T_i = My_i = Mz_i = 0.0
            Vy_j = Vz_j = T_j = My_j = Mz_j = 0.0
            sigma_i = sigma_j = abs(axial_signed) / A
            sigma_tension_max = max(axial_signed / A, 0.0)
            sigma_compression_min = min(axial_signed / A, 0.0)
            sigma_extreme_signed = axial_signed / A
            sigma_max = sigma_i
            max_M_resultant = 0.0
            stress_basis = "truss axial only"
        else:
            raise ValueError(f"Unknown element_type {etype!r} for element {eid}")

        if axial_signed > options.axial_state_tolerance_N:
            axial_state = "tension"
        elif axial_signed < -options.axial_state_tolerance_N:
            axial_state = "compression"
        else:
            axial_state = "near zero"

        records.append(
            {
                "element_id": eid,
                "start_node": ni,
                "end_node": nj,
                "start_label": e["start_label"],
                "end_label": e["end_label"],
                "member_role": e["member_role"],
                "element_type": etype,
                "section_name": e["section_name"],
                "length_mm": L,
                "N_i": N_i,
                "N_j": N_j,
                "Vy_i": Vy_i, "Vy_j": Vy_j,
                "Vz_i": Vz_i, "Vz_j": Vz_j,
                "T_i": T_i, "T_j": T_j,
                "My_i": My_i, "My_j": My_j,
                "Mz_i": Mz_i, "Mz_j": Mz_j,
                "axial_force_signed_N": axial_signed,
                "axial_state": axial_state,
                "max_abs_N_N": abs(axial_signed),
                "max_M_resultant_Nmm": max_M_resultant,
                "sigma_axial_signed_MPa": axial_signed / A,
                "sigma_i_abs_MPa": sigma_i,
                "sigma_j_abs_MPa": sigma_j,
                "sigma_tension_max_MPa": sigma_tension_max,
                "sigma_compression_min_MPa": sigma_compression_min,
                "sigma_extreme_signed_MPa": sigma_extreme_signed,
                "sigma_max_abs_MPa": sigma_max,
                "stress_basis": stress_basis,
                "A_mm2": A,
                "Iy_mm4": Iy,
                "Iz_mm4": Iz,
                "J_mm4": float(sec["J_mm4"]),
            }
        )

    return pd.DataFrame(records)


# ---------------------------------------------------------------------------
# Convenience summary/export functions
# ---------------------------------------------------------------------------


def summarize_result(result: Mapping[str, object]) -> Dict[str, float | str | int]:
    """Return one compact summary row for a solved model."""
    disp: pd.DataFrame = result["displacements"]  # type: ignore[assignment]
    mf: pd.DataFrame = result["member_forces"]    # type: ignore[assignment]
    return {
        "model": str(result["model"]),
        "crossing_diagonal_mode": str(result["options"].crossing_diagonal_mode),
        "section_mode": str(result["options"].section_mode),
        "nodes": int(len(result["nodes"])),
        "elements": int(len(result["elements"])),
        "frame_elements": int((result["elements"]["element_type"] == "frame").sum()),
        "truss_elements": int((result["elements"]["element_type"] == "truss").sum()),
        "restrained_dofs": int(len(result["restrained_dofs"])),
        "free_dofs": int(len(result["free_dofs"])),
        "reduced_stiffness_rank": int(result["reduced_stiffness_rank"]),
        "max_translation_mm": float(disp["translation_mag_mm"].max()),
        "max_abs_Uz_mm": float(disp["Uz_mm"].abs().max()),
        "max_abs_axial_force_kN": float(mf["max_abs_N_N"].max() / 1000.0),
        "max_bending_resultant_kNmm": float(mf["max_M_resultant_Nmm"].max() / 1000.0),
        "max_sigma_MPa": float(mf["sigma_max_abs_MPa"].max()),
    }


def critical_member_table(results: Mapping[str, Mapping[str, object]], n: int = 10) -> pd.DataFrame:
    """Top-N critical members by sigma_max_abs_MPa for each model."""
    tables = []
    for name, res in results.items():
        table = res["member_forces"].sort_values("sigma_max_abs_MPa", ascending=False).head(n).copy()
        table.insert(0, "model", name)
        tables.append(table)
    return pd.concat(tables, ignore_index=True)





def critical_member_table_matched_alternating(results: Mapping[str, Mapping[str, object]], n: int = 10) -> pd.DataFrame:
    """
    Return a matched, alternating critical-member table.

    Procedure:
        1. Find the top `n` elements by sigma_max_abs_MPa within each model.
        2. Take the union of those element_ids.
        3. Pull those element_ids from every model where they exist.
        4. Sort by element_id first, model second, so analogous IDs appear as
           adjacent rows: TE5 e0, Z5 e0, TE5 e1, Z5 e1, etc.

    The `in_top_n_for_model` flag tells you whether that row was actually in
    that model's top-n list, or whether it was included only as the matching
    counterpart of a critical element from another model.
    """
    model_order = list(results.keys())
    top_ids_by_model: Dict[str, set[int]] = {}
    critical_ids: set[int] = set()
    for name, res in results.items():
        mf = res["member_forces"]  # type: ignore[index]
        top = mf.sort_values("sigma_max_abs_MPa", ascending=False).head(n).copy()
        ids = set(top["element_id"].astype(int))
        top_ids_by_model[name] = ids
        critical_ids |= ids

    rows = []
    for name, res in results.items():
        mf = res["member_forces"].copy()  # type: ignore[index]
        part = mf.loc[mf["element_id"].astype(int).isin(critical_ids)].copy()
        part.insert(0, "model", name)
        part["in_top_n_for_model"] = part["element_id"].astype(int).isin(top_ids_by_model[name])
        rows.append(part)

    if not rows:
        return pd.DataFrame()
    out = pd.concat(rows, ignore_index=True)
    out["model"] = pd.Categorical(out["model"], categories=model_order, ordered=True)
    return out.sort_values(["element_id", "model"]).reset_index(drop=True)


def alternating_element_table(results: Mapping[str, Mapping[str, object]]) -> pd.DataFrame:
    """
    Return all element rows with models alternating by element_id.

    Sort order:
        element_id first, model second using the insertion order of results.

    If one model has extra element IDs, they naturally appear at the end of the
    relevant element_id groups.
    """
    model_order = list(results.keys())
    tables = []
    for name, res in results.items():
        table = res["member_forces"].copy()
        table.insert(0, "model", name)
        tables.append(table)
    out = pd.concat(tables, ignore_index=True)
    out["model"] = pd.Categorical(out["model"], categories=model_order, ordered=True)
    return out.sort_values(["element_id", "model"]).reset_index(drop=True)


def export_result(result: Mapping[str, object], out_dir: str | Path, prefix: Optional[str] = None) -> None:
    """Write standard result CSV files for one solved model."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    if prefix is None:
        prefix = str(result["model"])
    for key, filename in [
        ("displacements", "displacements.csv"),
        ("reactions", "reactions.csv"),
        ("member_forces", "member_forces.csv"),
        ("section_properties", "section_properties.csv"),
        ("loads_used", "loads.csv"),
        ("validation", "validation.csv"),
        ("stabilization", "stabilization.csv"),
        ("crossing_mode_diagnostics", "crossing_mode_diagnostics.csv"),
    ]:
        df = result[key]
        if isinstance(df, pd.DataFrame):
            df.to_csv(out_dir / f"{prefix}_{filename}", index=False)




def export_results_to_excel(
    results: Mapping[str, Mapping[str, object]],
    out_xlsx: str | Path,
    top_n: int = 10,
    display_columns: Optional[Sequence[str]] = None,
) -> None:
    """
    Export a multi-model result package to one XLSX workbook.

    This is intentionally placed in solver.py so the notebook can call one
    function and get a complete report artifact. The function uses pandas'
    ExcelWriter because it is portable in ordinary Codespaces/Jupyter setups.

    Sheets created:
        Summary
        Critical_Ranked
        Critical_Matched
        All_Elements_Alternating
        <Model>_Displacements
        <Model>_Reactions
        <Model>_Loads
        <Model>_Validation
        <Model>_Stabilization
        <Model>_Crossing_Mode
        <Model>_Sections

    Notes:
        * Critical_Matched uses matched/alternating element_id ordering.
        * The workbook is an engineering exploration output, not a formal
          design report.
    """
    out_xlsx = Path(out_xlsx)
    out_xlsx.parent.mkdir(parents=True, exist_ok=True)

    if display_columns is None:
        display_columns = [
            "model", "element_id", "element_type", "start_label", "end_label",
            "member_role", "section_name", "axial_state", "sigma_extreme_signed_MPa",
            "sigma_max_abs_MPa", "sigma_axial_signed_MPa", "axial_force_signed_N",
            "max_abs_N_N", "max_M_resultant_Nmm", "stress_basis",
        ]

    summary = pd.DataFrame([summarize_result(res) for res in results.values()])
    critical_ranked = critical_member_table(results, n=top_n)
    critical_matched = critical_member_table_matched_alternating(results, n=top_n)
    all_alternating = alternating_element_table(results)

    def cols_existing(df: pd.DataFrame, cols: Sequence[str]) -> list[str]:
        return [c for c in cols if c in df.columns]

    with pd.ExcelWriter(out_xlsx, engine="openpyxl") as writer:
        summary.to_excel(writer, sheet_name="Summary", index=False)
        critical_ranked[cols_existing(critical_ranked, display_columns)].to_excel(writer, sheet_name="Critical_Ranked", index=False)
        matched_cols = list(display_columns)
        if "in_top_n_for_model" not in matched_cols:
            matched_cols.append("in_top_n_for_model")
        critical_matched[cols_existing(critical_matched, matched_cols)].to_excel(writer, sheet_name="Critical_Matched", index=False)
        all_alternating[cols_existing(all_alternating, display_columns)].to_excel(writer, sheet_name="All_Elements_Alternating", index=False)

        for name, res in results.items():
            safe = str(name)[:20]
            res["displacements"].to_excel(writer, sheet_name=f"{safe}_Displacements", index=False)  # type: ignore[index]
            res["reactions"].to_excel(writer, sheet_name=f"{safe}_Reactions", index=False)          # type: ignore[index]
            res["loads_used"].to_excel(writer, sheet_name=f"{safe}_Loads", index=False)             # type: ignore[index]
            res["validation"].to_excel(writer, sheet_name=f"{safe}_Validation", index=False)        # type: ignore[index]
            res["stabilization"].to_excel(writer, sheet_name=f"{safe}_Stabilization", index=False)  # type: ignore[index]
            res["crossing_mode_diagnostics"].to_excel(writer, sheet_name=f"{safe}_Crossing_Mode", index=False)  # type: ignore[index]
            res["section_properties"].to_excel(writer, sheet_name=f"{safe}_Sections", index=False)  # type: ignore[index]

        # Simple formatting: freeze panes, autofilter, and reasonable widths.
        wb = writer.book
        for ws in wb.worksheets:
            ws.freeze_panes = "A2"
            if ws.max_row >= 1 and ws.max_column >= 1:
                ws.auto_filter.ref = ws.dimensions
            for col_cells in ws.columns:
                letter = col_cells[0].column_letter
                max_len = 0
                for cell in col_cells[: min(len(col_cells), 80)]:
                    value = cell.value
                    if value is not None:
                        max_len = max(max_len, len(str(value)))
                ws.column_dimensions[letter].width = min(max(max_len + 2, 10), 38)


# ---------------------------------------------------------------------------
# Simple command-line entrypoint for quick checks
# ---------------------------------------------------------------------------


def _demo_cli() -> None:
    """Minimal command-line interface for quick local tests."""
    import argparse

    parser = argparse.ArgumentParser(description="Solve one mixed frame/truss tower model.")
    parser.add_argument("--name", required=True)
    parser.add_argument("--nodes", required=True)
    parser.add_argument("--elements", required=True)
    parser.add_argument("--loads", required=True)
    parser.add_argument("--out", default="results_cli")
    parser.add_argument("--section-mode", choices=["principal", "isotropic"], default="principal")
    parser.add_argument("--crossing-diagonal-mode", choices=sorted(CROSSING_DIAGONAL_MODES), default="continuous")
    parser.add_argument("--xlsx", default=None, help="Optional XLSX report path for this single model.")
    args = parser.parse_args()

    model = read_model(args.name, args.nodes, args.elements)
    options = SolveOptions(section_mode=args.section_mode, crossing_diagonal_mode=args.crossing_diagonal_mode)
    loads = read_loads_csv(args.loads)
    result = solve_model(model, loads, options)
    export_result(result, args.out)
    if args.xlsx:
        export_results_to_excel({args.name: result}, args.xlsx)
    print(pd.DataFrame([summarize_result(result)]).to_string(index=False))


if __name__ == "__main__":
    _demo_cli()
