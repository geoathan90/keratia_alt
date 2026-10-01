"""
helpers.py
==========

Small plotting and notebook convenience functions for the mixed frame/truss
lattice-tower exploration project.

The structural solver itself lives in solver.py. This file intentionally keeps
plotting and display helpers separate so solver.py remains easier to inspect.
"""

from __future__ import annotations

from pathlib import Path
from typing import Iterable, Mapping, Optional

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from mpl_toolkits.mplot3d.art3d import Line3DCollection


DEFAULT_DISPLAY_COLS = [
    "model",
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


def select_elements(elements: pd.DataFrame, element_type: str = "all") -> pd.DataFrame:
    """
    Filter elements by type for visualization.

    element_type:
        "all"   -> all elements
        "frame" -> only frame elements
        "truss" -> only truss elements
    """
    elements = elements.copy()
    if "element_type" not in elements.columns:
        elements["element_type"] = "frame"
    if element_type == "all":
        return elements
    return elements.loc[elements["element_type"].fillna("frame").eq(element_type)].copy()


def select_member_forces_for_elements(member_forces: pd.DataFrame, elements_subset: pd.DataFrame) -> pd.DataFrame:
    """Filter member_forces to the element IDs present in elements_subset."""
    keep_ids = set(elements_subset["element_id"].astype(int))
    return member_forces.loc[member_forces["element_id"].astype(int).isin(keep_ids)].copy()


def _set_axes_equal(ax):
    """Make 3D axes use roughly equal scale."""
    x_limits = ax.get_xlim3d()
    y_limits = ax.get_ylim3d()
    z_limits = ax.get_zlim3d()
    x_range = abs(x_limits[1] - x_limits[0])
    y_range = abs(y_limits[1] - y_limits[0])
    z_range = abs(z_limits[1] - z_limits[0])
    max_range = max(x_range, y_range, z_range)
    x_mid = np.mean(x_limits)
    y_mid = np.mean(y_limits)
    z_mid = np.mean(z_limits)
    ax.set_xlim3d([x_mid - max_range/2, x_mid + max_range/2])
    ax.set_ylim3d([y_mid - max_range/2, y_mid + max_range/2])
    ax.set_zlim3d([z_mid - max_range/2, z_mid + max_range/2])


def plot_geometry(
    nodes: pd.DataFrame,
    elements: pd.DataFrame,
    title: str = "Geometry",
    annotate: bool = True,
    annotate_elements: bool = True,
    show_frame: bool = True,
    show_truss: bool = True,
    figsize=(9, 8),
):
    """
    Plot 3D geometry with frame/truss members visually separated.

    Uses matplotlib only, so it works reliably in browser notebooks/Codespaces.
    """
    nodes_by_id = nodes.set_index("node_id")
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111, projection="3d")

    for etype, lw, alpha, label in [
        ("frame", 2.2, 0.95, "frame"),
        ("truss", 1.0, 0.70, "truss"),
    ]:
        if etype == "frame" and not show_frame:
            continue
        if etype == "truss" and not show_truss:
            continue
        subset = elements.loc[elements["element_type"].fillna("frame").eq(etype)]
        segments = []
        for _, e in subset.iterrows():
            a = nodes_by_id.loc[int(e["start_node"]), ["x", "y", "z"]].to_numpy(float)
            b = nodes_by_id.loc[int(e["end_node"]), ["x", "y", "z"]].to_numpy(float)
            segments.append([a, b])
            if annotate_elements:
                mid = 0.5 * (a + b)
                ax.text(mid[0], mid[1], mid[2], str(int(e["element_id"])), fontsize=7)
        if segments:
            coll = Line3DCollection(segments, linewidths=lw, alpha=alpha, label=label)
            ax.add_collection3d(coll)

    ax.scatter(nodes["x"], nodes["y"], nodes["z"], s=14)
    if annotate:
        for _, n in nodes.iterrows():
            ax.text(n["x"], n["y"], n["z"], str(n["node_label"]), fontsize=7)

    # Highlight supports/load points if metadata columns are present.
    if "is_fixed_support" in nodes.columns:
        mask = nodes["is_fixed_support"].fillna(False).astype(str).str.lower().isin(["true", "1", "yes"])
        # Also highlight any min-z nodes, because supports are normally auto-detected that way.
        mask = mask | ((nodes["z"] - nodes["z"].min()).abs() < 1e-6)
        if mask.any():
            ax.scatter(nodes.loc[mask, "x"], nodes.loc[mask, "y"], nodes.loc[mask, "z"], marker="s", s=40, label="fixed/min-z")
    if "is_load_point" in nodes.columns:
        mask = nodes["is_load_point"].fillna(False).astype(str).str.lower().isin(["true", "1", "yes"])
        if mask.any():
            ax.scatter(nodes.loc[mask, "x"], nodes.loc[mask, "y"], nodes.loc[mask, "z"], marker="^", s=50, label="load point")

    ax.set_title(title)
    ax.set_xlabel("X [mm]")
    ax.set_ylabel("Y [mm]")
    ax.set_zlabel("Z [mm]")
    ax.legend(loc="best")
    _set_axes_equal(ax)
    return fig, ax


def plot_stress_3d(
    nodes: pd.DataFrame,
    elements: pd.DataFrame,
    member_forces: pd.DataFrame,
    title: str = "Stress visualization",
    annotate_load_points: bool = True,
    annotate_elements: bool = False,
    show_frame: bool = True,
    show_truss: bool = True,
    figsize=(9, 8),
):
    """Plot members colored by sigma_max_abs_MPa."""
    nodes_by_id = nodes.set_index("node_id")
    forces = member_forces.set_index("element_id")
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111, projection="3d")

    plot_elems = elements.copy()
    if "element_type" not in plot_elems.columns:
        plot_elems["element_type"] = "frame"
    mask = pd.Series(False, index=plot_elems.index)
    if show_frame:
        mask |= plot_elems["element_type"].eq("frame")
    if show_truss:
        mask |= plot_elems["element_type"].eq("truss")
    plot_elems = plot_elems.loc[mask]

    segments = []
    values = []
    linewidths = []
    for _, e in plot_elems.iterrows():
        eid = int(e["element_id"])
        if eid not in forces.index:
            continue
        a = nodes_by_id.loc[int(e["start_node"]), ["x", "y", "z"]].to_numpy(float)
        b = nodes_by_id.loc[int(e["end_node"]), ["x", "y", "z"]].to_numpy(float)
        segments.append([a, b])
        values.append(float(forces.loc[eid, "sigma_max_abs_MPa"]))
        linewidths.append(2.4 if e["element_type"] == "frame" else 1.2)
        if annotate_elements:
            mid = 0.5 * (a + b)
            ax.text(mid[0], mid[1], mid[2], str(eid), fontsize=7)

    if segments:
        coll = Line3DCollection(segments, array=np.array(values), cmap="viridis", linewidths=linewidths)
        ax.add_collection3d(coll)
        cbar = fig.colorbar(coll, ax=ax, shrink=0.65, pad=0.10)
        cbar.set_label("sigma_max_abs_MPa")

    ax.scatter(nodes["x"], nodes["y"], nodes["z"], s=8)
    if annotate_load_points and "is_load_point" in nodes.columns:
        mask = nodes["is_load_point"].fillna(False).astype(str).str.lower().isin(["true", "1", "yes"])
        for _, n in nodes.loc[mask].iterrows():
            ax.text(n["x"], n["y"], n["z"], str(n["node_label"]), fontsize=8)
            ax.scatter([n["x"]], [n["y"]], [n["z"]], marker="^", s=50)

    ax.set_title(title)
    ax.set_xlabel("X [mm]")
    ax.set_ylabel("Y [mm]")
    ax.set_zlabel("Z [mm]")
    _set_axes_equal(ax)
    return fig, ax



def plot_signed_stress_3d(
    nodes: pd.DataFrame,
    elements: pd.DataFrame,
    member_forces: pd.DataFrame,
    title: str = "Signed tension/compression stress visualization",
    stress_column: str = "sigma_extreme_signed_MPa",
    annotate_load_points: bool = True,
    annotate_elements: bool = False,
    show_frame: bool = True,
    show_truss: bool = True,
    figsize=(9, 8),
):
    """
    Plot members using a diverging color map for signed stress.

    Positive values are tensile. Negative values are compressive.

    Recommended stress columns:
        sigma_extreme_signed_MPa
            For frame elements, this is the signed controlling extreme-fiber
            stress from axial + bending. For truss elements, it is the axial
            stress. Good for seeing which members are governed by tension or
            compression at the maximum-stress fiber.

        sigma_axial_signed_MPa
            Axial stress only for both frame and truss elements. Good for seeing
            the pure axial tension/compression load path.
    """
    nodes_by_id = nodes.set_index("node_id")
    forces = member_forces.set_index("element_id")
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111, projection="3d")

    plot_elems = elements.copy()
    if "element_type" not in plot_elems.columns:
        plot_elems["element_type"] = "frame"
    mask = pd.Series(False, index=plot_elems.index)
    if show_frame:
        mask |= plot_elems["element_type"].eq("frame")
    if show_truss:
        mask |= plot_elems["element_type"].eq("truss")
    plot_elems = plot_elems.loc[mask]

    segments = []
    values = []
    linewidths = []
    for _, e in plot_elems.iterrows():
        eid = int(e["element_id"])
        if eid not in forces.index or stress_column not in forces.columns:
            continue
        a = nodes_by_id.loc[int(e["start_node"]), ["x", "y", "z"]].to_numpy(float)
        b = nodes_by_id.loc[int(e["end_node"]), ["x", "y", "z"]].to_numpy(float)
        segments.append([a, b])
        values.append(float(forces.loc[eid, stress_column]))
        linewidths.append(2.4 if e["element_type"] == "frame" else 1.2)
        if annotate_elements:
            mid = 0.5 * (a + b)
            ax.text(mid[0], mid[1], mid[2], str(eid), fontsize=7)

    if segments:
        values_arr = np.array(values, dtype=float)
        vmax = max(float(np.nanmax(np.abs(values_arr))), 1e-12)
        coll = Line3DCollection(
            segments,
            array=values_arr,
            cmap="coolwarm_r",
            linewidths=linewidths,
            clim=(-vmax, vmax),
        )
        ax.add_collection3d(coll)
        cbar = fig.colorbar(coll, ax=ax, shrink=0.65, pad=0.10)
        cbar.set_label(f"{stress_column} (+ tension, - compression) [MPa]")

    ax.scatter(nodes["x"], nodes["y"], nodes["z"], s=8)
    if annotate_load_points and "is_load_point" in nodes.columns:
        mask = nodes["is_load_point"].fillna(False).astype(str).str.lower().isin(["true", "1", "yes"])
        for _, n in nodes.loc[mask].iterrows():
            ax.text(n["x"], n["y"], n["z"], str(n["node_label"]), fontsize=8)
            ax.scatter([n["x"]], [n["y"]], [n["z"]], marker="^", s=50)

    ax.set_title(title)
    ax.set_xlabel("X [mm]")
    ax.set_ylabel("Y [mm]")
    ax.set_zlabel("Z [mm]")
    _set_axes_equal(ax)
    return fig, ax


def plot_deformed_shape(
    nodes: pd.DataFrame,
    elements: pd.DataFrame,
    displacements: pd.DataFrame,
    scale: float = 100.0,
    title: str = "Deformed shape",
    show_truss: bool = True,
    figsize=(9, 8),
):
    """Plot undeformed and exaggerated deformed shape."""
    d = displacements.set_index("node_id")
    def node_xyz(nid):
        r = nodes.loc[nodes["node_id"] == nid].iloc[0]
        return np.array([r["x"], r["y"], r["z"]], dtype=float)
    def node_xyz_def(nid):
        base = node_xyz(nid)
        u = d.loc[nid, ["Ux_mm", "Uy_mm", "Uz_mm"]].to_numpy(float)
        return base + scale * u

    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111, projection="3d")
    for _, e in elements.iterrows():
        if e["element_type"] == "truss" and not show_truss:
            continue
        a = node_xyz(int(e["start_node"])); b = node_xyz(int(e["end_node"]))
        ad = node_xyz_def(int(e["start_node"])); bd = node_xyz_def(int(e["end_node"]))
        ax.plot([a[0], b[0]], [a[1], b[1]], [a[2], b[2]], alpha=0.25, linewidth=0.8)
        ax.plot([ad[0], bd[0]], [ad[1], bd[1]], [ad[2], bd[2]], linewidth=1.4)
    ax.set_title(f"{title} — scale ×{scale:g}")
    ax.set_xlabel("X [mm]"); ax.set_ylabel("Y [mm]"); ax.set_zlabel("Z [mm]")
    _set_axes_equal(ax)
    return fig, ax


def display_columns() -> list[str]:
    return DEFAULT_DISPLAY_COLS.copy()
