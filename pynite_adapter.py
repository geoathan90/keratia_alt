"""
pynite_adapter.py
=================

Experimental PyNite backend for the keratia lattice tower models.

This module translates the existing CSV/ModelData workflow into a PyNite
`FEModel3D`, runs linear-elastic nodal-load cases, and converts the results
back into tables shaped similarly to our educational stiffness solver outputs.

Important modelling note
------------------------
PyNite members are beam-column/physical-member objects. PyNite does not expose
exactly the same custom mixed frame/truss element formulation as our local
solver. In this adapter:

* CSV `frame` elements are added as ordinary PyNite members.
* CSV `truss` elements are represented as pin-ended PyNite members by releasing
  local Ry/Rz rotations at both ends.

That is a useful engineering comparison, but it is not mathematically identical
in every detail to the in-house truss element. Treat this backend as a second
analysis route/check, not as a bit-for-bit replacement.

Recommended comparison setting
------------------------------
For direct solver-to-solver comparisons, start with:

    section_mode = "isotropic"
    crossing_diagonal_mode = "continuous"

The isotropic section mode removes most local-axis-orientation sensitivity.
PyNite uses its own local-axis conventions and a member rotation angle, while
our educational solver builds local axes from a tower-inward direction.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping, Optional, Sequence

import math
import numpy as np
import pandas as pd

from solver import (
    ModelData,
    SolveOptions,
    build_section_table,
    merge_crossing_truss_splits_to_continuous,
    normalize_crossing_diagonal_mode,
    section_properties_for_element,
)

try:
    from Pynite import FEModel3D
except Exception as exc:  # pragma: no cover - only triggered when dependency missing
    FEModel3D = None
    _PYNITE_IMPORT_ERROR = exc
else:
    _PYNITE_IMPORT_ERROR = None


@dataclass(frozen=True)
class PyNiteRunOptions:
    """Options for the PyNite translation layer."""

    section_mode: str = "isotropic"
    principal_orientation: str = "strong_inward"
    crossing_diagonal_mode: str = "continuous"
    fix_min_z_nodes: bool = True
    fix_nodes_marked_as_support: bool = True
    release_truss_end_bending: bool = True
    auto_fix_truss_only_rotations: bool = True
    min_z_tolerance_mm: float = 1.0e-6
    axial_state_tolerance_N: float = 1.0
    stress_sample_points: int = 7


def _require_pynite() -> None:
    """Raise a clear installation error if PyNite is unavailable."""
    if FEModel3D is None:
        raise ImportError(
            "PyNite is not installed. Install it with:\n\n"
            "    python -m pip install --user PyniteFEA\n\n"
            "or install all repo requirements with:\n\n"
            "    python -m pip install --user -r requirements.txt"
        ) from _PYNITE_IMPORT_ERROR


def _truthy(value: object) -> bool:
    """Interpret common CSV boolean-like values."""
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"true", "1", "yes", "y", "fixed", "support"}


def _safe_combo_get(d: Mapping[str, float], combo_name: str) -> float:
    """Read a PyNite result dictionary; return 0.0 if the key is missing."""
    try:
        return float(d.get(combo_name, 0.0))
    except AttributeError:
        try:
            return float(d[combo_name])
        except Exception:
            return 0.0


def _pynite_node_name(node_id: int) -> str:
    return f"N{int(node_id)}"


def _pynite_member_name(element_id: int) -> str:
    return f"E{int(element_id)}"


def _working_model_for_pynite(model: ModelData, options: PyNiteRunOptions) -> tuple[ModelData, pd.DataFrame]:
    """Apply only the crossing-diagonal preprocessing needed before PyNite build."""
    mode = normalize_crossing_diagonal_mode(options.crossing_diagonal_mode)
    diagnostics = pd.DataFrame([{"mode": mode, "note": "no preprocessing"}])
    working_model = model

    if mode == "continuous":
        working_model, diagnostics = merge_crossing_truss_splits_to_continuous(model)
    elif mode == "connected_unstabilized":
        diagnostics = pd.DataFrame([{
            "mode": mode,
            "note": (
                "connected_unstabilized is a diagnostic mode. PyNite may report "
                "an instability for truss-only crossing nodes."
            ),
        }])
    elif mode == "connected_stabilized":
        diagnostics = pd.DataFrame([{
            "mode": mode,
            "note": (
                "The PyNite adapter does not currently reproduce the in-house "
                "tiny directional stabilization springs. Use continuous for the "
                "initial PyNite comparison, or use the in-house solver for the "
                "connected_stabilized assumption."
            ),
        }])

    return working_model, diagnostics


def _support_flags_for_nodes(model: ModelData, options: PyNiteRunOptions) -> dict[int, tuple[bool, bool, bool, bool, bool, bool]]:
    """Return support flags DX,DY,DZ,RX,RY,RZ by node_id."""
    nodes = model.nodes
    elements = model.elements
    z_min = float(nodes["z"].min())

    fixed_nodes: set[int] = set()
    if options.fix_min_z_nodes:
        fixed_nodes.update(
            int(r.node_id) for r in nodes.itertuples(index=False)
            if abs(float(r.z) - z_min) <= options.min_z_tolerance_mm
        )
    if options.fix_nodes_marked_as_support and "is_fixed_support" in nodes.columns:
        fixed_nodes.update(
            int(row["node_id"]) for _, row in nodes.iterrows()
            if _truthy(row.get("is_fixed_support", False))
        )

    truss_only_nodes: set[int] = set()
    if options.auto_fix_truss_only_rotations:
        incident: dict[int, list[str]] = {int(n): [] for n in nodes["node_id"]}
        for _, e in elements.iterrows():
            etype = str(e.get("element_type", "frame")).strip().lower()
            incident[int(e["start_node"])].append(etype)
            incident[int(e["end_node"])].append(etype)
        for node_id, types in incident.items():
            if types and all(t == "truss" for t in types):
                truss_only_nodes.add(node_id)

    flags: dict[int, tuple[bool, bool, bool, bool, bool, bool]] = {}
    for node_id in nodes["node_id"].astype(int):
        if int(node_id) in fixed_nodes:
            flags[int(node_id)] = (True, True, True, True, True, True)
        elif int(node_id) in truss_only_nodes:
            flags[int(node_id)] = (False, False, False, True, True, True)
        else:
            flags[int(node_id)] = (False, False, False, False, False, False)
    return flags


def build_pynite_model(
    model: ModelData,
    load_rows: Sequence[Mapping[str, float | str]],
    options: Optional[PyNiteRunOptions] = None,
    case_name: str = "Case_1",
    combo_name: Optional[str] = None,
):
    """Build and return a PyNite FEModel3D plus translated metadata."""
    _require_pynite()
    options = options or PyNiteRunOptions()
    combo_name = combo_name or case_name

    working_model, crossing_diag = _working_model_for_pynite(model, options)
    solver_options = SolveOptions(
        section_mode=options.section_mode,
        principal_orientation=options.principal_orientation,
        crossing_diagonal_mode=options.crossing_diagonal_mode,
        fix_min_z_nodes=options.fix_min_z_nodes,
        fix_nodes_marked_as_support=options.fix_nodes_marked_as_support,
        auto_fix_truss_only_rotations=options.auto_fix_truss_only_rotations,
        axial_state_tolerance_N=options.axial_state_tolerance_N,
    )

    fe = FEModel3D()

    for _, n in working_model.nodes.iterrows():
        fe.add_node(_pynite_node_name(int(n["node_id"])), float(n["x"]), float(n["y"]), float(n["z"]))

    material_names: dict[tuple[float, float], str] = {}
    for _, e in working_model.elements.iterrows():
        E = float(e.get("E_N_per_mm2", 210_000.0) or 210_000.0)
        G = float(e.get("G_N_per_mm2", 81_000.0) or 81_000.0)
        key = (E, G)
        if key not in material_names:
            nu = max(0.0, min(0.49, E / (2.0 * G) - 1.0)) if G > 0 else 0.30
            mat_name = f"Steel_{len(material_names) + 1}"
            fe.add_material(mat_name, E=E, G=G, nu=nu, rho=0.0)
            material_names[key] = mat_name

    section_records = []
    for _, e in working_model.elements.iterrows():
        props = section_properties_for_element(e, solver_options)
        eid = int(e["element_id"])
        sec_name = f"S{eid}"
        fe.add_section(
            sec_name,
            A=float(props["A_mm2"]),
            Iy=float(props["Iy_mm4"]),
            Iz=float(props["Iz_mm4"]),
            J=float(props["J_mm4"]),
        )
        E = float(e.get("E_N_per_mm2", 210_000.0) or 210_000.0)
        G = float(e.get("G_N_per_mm2", 81_000.0) or 81_000.0)
        member_name = _pynite_member_name(eid)
        fe.add_member(
            member_name,
            _pynite_node_name(int(e["start_node"])),
            _pynite_node_name(int(e["end_node"])),
            material_names[(E, G)],
            sec_name,
        )
        if str(e.get("element_type", "frame")).lower() == "truss" and options.release_truss_end_bending:
            fe.def_releases(
                member_name,
                Dxi=False, Dyi=False, Dzi=False, Rxi=False, Ryi=True, Rzi=True,
                Dxj=False, Dyj=False, Dzj=False, Rxj=False, Ryj=True, Rzj=True,
            )
        section_records.append({
            "element_id": eid,
            "pynite_member": member_name,
            "pynite_section": sec_name,
            **props,
        })

    for node_id, flags in _support_flags_for_nodes(working_model, options).items():
        fe.def_support(
            _pynite_node_name(node_id),
            support_DX=flags[0], support_DY=flags[1], support_DZ=flags[2],
            support_RX=flags[3], support_RY=flags[4], support_RZ=flags[5],
        )

    label_to_node = dict(zip(working_model.nodes["node_label"], working_model.nodes["node_id"]))
    loads_used = []
    for lr in load_rows:
        label = str(lr["node_label"])
        if label not in label_to_node:
            raise ValueError(f"Unknown node_label in load case: {label!r}")
        node_name = _pynite_node_name(int(label_to_node[label]))
        for user_key, pynite_dir in [
            ("Fx", "FX"), ("Fy", "FY"), ("Fz", "FZ"),
            ("Mx", "MX"), ("My", "MY"), ("Mz", "MZ"),
        ]:
            value = float(lr.get(user_key, 0.0) or 0.0)
            if abs(value) > 0.0:
                fe.add_node_load(node_name, pynite_dir, value, case=case_name)
            loads_used.append({
                "node_label": label,
                "node_id": int(label_to_node[label]),
                "pynite_node": node_name,
                "dof": user_key,
                "pynite_direction": pynite_dir,
                "value": value,
                "case": case_name,
            })
    fe.add_load_combo(combo_name, {case_name: 1.0})

    section_table = pd.DataFrame(section_records)
    return fe, working_model, section_table, crossing_diag, pd.DataFrame(loads_used)


def analyze_pynite_model(
    model: ModelData,
    load_rows: Sequence[Mapping[str, float | str]],
    options: Optional[PyNiteRunOptions] = None,
    case_name: str = "Case_1",
    combo_name: Optional[str] = None,
) -> dict:
    """Build, solve, and convert a PyNite model into our standard result shape."""
    options = options or PyNiteRunOptions()
    combo_name = combo_name or case_name
    fe, working_model, section_table, crossing_diag, loads_used = build_pynite_model(
        model, load_rows, options=options, case_name=case_name, combo_name=combo_name
    )
    fe.analyze_linear(log=False, check_stability=True, check_statics=False, sparse=True)

    displacements = pynite_displacements_dataframe(fe, working_model, combo_name)
    reactions = pynite_reactions_dataframe(fe, working_model, combo_name)
    member_forces = pynite_member_forces_dataframe(fe, working_model, section_table, options, combo_name)

    validation = pd.DataFrame([{
        "model": working_model.name,
        "backend": "PyNite",
        "case": case_name,
        "combo": combo_name,
        "nodes": len(working_model.nodes),
        "elements": len(working_model.elements),
        "frame_elements": int((working_model.elements["element_type"] == "frame").sum()),
        "truss_like_members": int((working_model.elements["element_type"] == "truss").sum()),
        "section_mode": options.section_mode,
        "crossing_diagonal_mode": options.crossing_diagonal_mode,
        "truss_representation": "pin-ended PyNite member with local Ry/Rz releases",
        "axial_sign_convention": "converted to keratia convention: positive=tension, negative=compression",
    }])

    return {
        "name": working_model.name,
        "backend": "PyNite",
        "pynite_model": fe,
        "nodes": working_model.nodes.copy(),
        "elements": working_model.elements.copy(),
        "displacements": displacements,
        "reactions": reactions,
        "member_forces": member_forces,
        "section_properties": section_table,
        "loads_used": loads_used,
        "validation": validation,
        "crossing_mode_diagnostics": crossing_diag,
        "options": options,
    }


def pynite_displacements_dataframe(fe, model: ModelData, combo_name: str) -> pd.DataFrame:
    rows = []
    for _, n in model.nodes.iterrows():
        node = fe.nodes[_pynite_node_name(int(n["node_id"]))]
        dx = _safe_combo_get(node.DX, combo_name)
        dy = _safe_combo_get(node.DY, combo_name)
        dz = _safe_combo_get(node.DZ, combo_name)
        rx = _safe_combo_get(node.RX, combo_name)
        ry = _safe_combo_get(node.RY, combo_name)
        rz = _safe_combo_get(node.RZ, combo_name)
        rows.append({
            "node_id": int(n["node_id"]),
            "node_label": n["node_label"],
            "x": float(n["x"]),
            "y": float(n["y"]),
            "z": float(n["z"]),
            "Ux_mm": dx,
            "Uy_mm": dy,
            "Uz_mm": dz,
            "Rx_rad": rx,
            "Ry_rad": ry,
            "Rz_rad": rz,
            "translation_mag_mm": float(math.sqrt(dx * dx + dy * dy + dz * dz)),
        })
    return pd.DataFrame(rows)


def pynite_reactions_dataframe(fe, model: ModelData, combo_name: str) -> pd.DataFrame:
    rows = []
    for _, n in model.nodes.iterrows():
        node = fe.nodes[_pynite_node_name(int(n["node_id"]))]
        rx = _safe_combo_get(node.RxnFX, combo_name)
        ry = _safe_combo_get(node.RxnFY, combo_name)
        rz = _safe_combo_get(node.RxnFZ, combo_name)
        mx = _safe_combo_get(node.RxnMX, combo_name)
        my = _safe_combo_get(node.RxnMY, combo_name)
        mz = _safe_combo_get(node.RxnMZ, combo_name)
        if max(abs(rx), abs(ry), abs(rz), abs(mx), abs(my), abs(mz)) > 1e-9:
            rows.append({
                "node_id": int(n["node_id"]),
                "node_label": n["node_label"],
                "Rxn_Fx_N": rx,
                "Rxn_Fy_N": ry,
                "Rxn_Fz_N": rz,
                "Rxn_Mx_Nmm": mx,
                "Rxn_My_Nmm": my,
                "Rxn_Mz_Nmm": mz,
                "force_reaction_resultant_N": float(math.sqrt(rx * rx + ry * ry + rz * rz)),
            })
    return pd.DataFrame(rows)


def _axial_state(N: float, tol: float) -> str:
    if N > tol:
        return "tension"
    if N < -tol:
        return "compression"
    return "near zero"


def pynite_member_forces_dataframe(
    fe,
    model: ModelData,
    section_table: pd.DataFrame,
    options: PyNiteRunOptions,
    combo_name: str,
) -> pd.DataFrame:
    """Recover approximate member result rows in the same style as solver.py."""
    sec = section_table.set_index("element_id")
    rows = []
    n_points = max(2, int(options.stress_sample_points))

    for _, e in model.elements.iterrows():
        eid = int(e["element_id"])
        member_name = _pynite_member_name(eid)
        m = fe.members[member_name]
        L = float(m.L())
        props = sec.loc[eid]
        A = float(props["A_mm2"])
        Iy = float(props["Iy_mm4"])
        Iz = float(props["Iz_mm4"])
        cy = float(props["c_y_mm"])
        cz = float(props["c_z_mm"])

        xs = np.linspace(0.0, L, n_points)
        N_vals = []
        My_vals = []
        Mz_vals = []
        sigma_candidates = []
        for x in xs:
            try:
                # PyNite's member.axial(...) sign is opposite to the convention
                # used throughout this project. Convert immediately so every
                # downstream table/plot keeps the keratia convention:
                #     positive = tension, negative = compression.
                N = -float(m.axial(x, combo_name))
            except Exception:
                N = 0.0
            try:
                My = float(m.moment("My", x, combo_name))
            except Exception:
                My = 0.0
            try:
                Mz = float(m.moment("Mz", x, combo_name))
            except Exception:
                Mz = 0.0
            N_vals.append(N)
            My_vals.append(My)
            Mz_vals.append(Mz)

            if str(e.get("element_type", "frame")).lower() == "truss":
                sigma_candidates.append(N / A if A else 0.0)
            else:
                for y in (-cy, cy):
                    for z in (-cz, cz):
                        sigma_candidates.append(N / A + Mz * y / Iz - My * z / Iy)

        N_arr = np.asarray(N_vals, dtype=float)
        My_arr = np.asarray(My_vals, dtype=float)
        Mz_arr = np.asarray(Mz_vals, dtype=float)
        sig_arr = np.asarray(sigma_candidates, dtype=float)

        N_signed = float(np.nanmean(N_arr)) if len(N_arr) else 0.0
        max_abs_N = float(np.nanmax(np.abs(N_arr))) if len(N_arr) else abs(N_signed)
        max_abs_My = float(np.nanmax(np.abs(My_arr))) if len(My_arr) else 0.0
        max_abs_Mz = float(np.nanmax(np.abs(Mz_arr))) if len(Mz_arr) else 0.0
        max_M_resultant = float(math.sqrt(max_abs_My * max_abs_My + max_abs_Mz * max_abs_Mz))
        if len(sig_arr):
            i_sig = int(np.nanargmax(np.abs(sig_arr)))
            sigma_signed = float(sig_arr[i_sig])
            sigma_abs = float(abs(sigma_signed))
        else:
            sigma_signed = 0.0
            sigma_abs = 0.0

        rows.append({
            "element_id": eid,
            "pynite_member": member_name,
            "element_type": e.get("element_type", "frame"),
            "start_node": int(e["start_node"]),
            "end_node": int(e["end_node"]),
            "start_label": e.get("start_label", ""),
            "end_label": e.get("end_label", ""),
            "member_role": e.get("member_role", ""),
            "section_name": e.get("section_name", ""),
            "length_mm": L,
            "axial_force_signed_N": N_signed,
            "max_abs_N_N": max_abs_N,
            "axial_state": _axial_state(N_signed, options.axial_state_tolerance_N),
            "sigma_axial_signed_MPa": N_signed / A if A else 0.0,
            "max_abs_My_Nmm": max_abs_My,
            "max_abs_Mz_Nmm": max_abs_Mz,
            "max_M_resultant_Nmm": max_M_resultant,
            "sigma_extreme_signed_MPa": sigma_signed,
            "sigma_max_abs_MPa": sigma_abs,
            "stress_basis": "PyNite member sampling; axial+bending for frames, axial for truss-like members",
        })
    return pd.DataFrame(rows)


def summarize_pynite_result(result: Mapping[str, object]) -> dict:
    """Compact summary similar to solver.summarize_result."""
    disp = pd.DataFrame(result["displacements"])
    mf = pd.DataFrame(result["member_forces"])
    elems = pd.DataFrame(result["elements"])
    return {
        "model": result.get("name", ""),
        "backend": "PyNite",
        "nodes": int(len(result["nodes"])),
        "elements": int(len(elems)),
        "frame_elements": int((elems["element_type"] == "frame").sum()),
        "truss_like_members": int((elems["element_type"] == "truss").sum()),
        "max_translation_mm": float(disp["translation_mag_mm"].max()),
        "max_abs_Uz_mm": float(disp["Uz_mm"].abs().max()),
        "max_abs_axial_force_kN": float(mf["max_abs_N_N"].max() / 1000.0),
        "max_bending_resultant_kNmm": float(mf["max_M_resultant_Nmm"].max() / 1000.0),
        "max_sigma_MPa": float(mf["sigma_max_abs_MPa"].max()),
    }
