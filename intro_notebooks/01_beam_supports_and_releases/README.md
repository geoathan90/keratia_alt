# 01 — Beam supports and moment releases

## Model

- beam length: 6000 mm
- midpoint node: 3000 mm
- midpoint vertical load: -10,000 N in global Y
- material: E = 210,000 N/mm²
- bending inertia for the in-plane problem: Iz = 80,000,000 mm⁴

The beam is represented by two ordinary PyNite members, `M1` and `M2`, joined rigidly at node `MID`.

## Input files

- `nodes.csv`: coordinates only.
- `materials.csv`: PyNite material properties.
- `sections.csv`: section properties.
- `members.csv`: member connectivity and material/section assignment.
- `configurations.csv`: names and descriptions of the five structural configurations.
- `supports.csv`: node restraint flags, by configuration.
- `member_releases.csv`: only the active member-end releases. Missing rows mean no releases.
- `nodal_loads.csv`: the physical load case, common to every structural configuration.

### Boolean convention

The support and release CSVs use `0` = false and `1` = true, but the meaning depends on the field:

- `support_RZ = 1` means node rotation RZ is **restrained**.
- `Rzi = 1` in the release CSV means the member-end RZ transfer is **released**.

## Why the out-of-plane DOFs are restrained

PyNite is a 3D frame solver. This lesson is intentionally 2D in the global X-Y plane. Therefore every node is restrained in:

- DZ
- RX
- RY

These are plane-enforcement restraints, not intended as physical supports along the beam.

The in-plane DOFs used in the lesson are:

- DX — horizontal translation
- DY — vertical translation
- RZ — in-plane rotation

## Configuration map

| configuration | node/support idea | member-end release idea |
|---|---|---|
| `FF` | both ends fixed | none |
| `PR` | A pin, B roller | none |
| `FP_SUPPORT` | A fixed, B rotationally free | none |
| `FF_REL_RIGHT` | A and B nodes fixed | release M2 RZj at B |
| `FF_REL_BOTH` | A and B nodes fixed | release M1 RZi at A and M2 RZj at B |

The two especially useful comparisons are:

- `FP_SUPPORT` vs `FF_REL_RIGHT`
- `PR` vs `FF_REL_BOTH`

Under this purely vertical load, each pair should have essentially the same bending response even though the modelling statements are different.

## Hand-calculation reference

For a midpoint load P on a beam of length L:

- fixed-fixed: each vertical reaction is P/2 and each end-moment magnitude is PL/8;
- pin-roller: each vertical reaction is P/2 and the end moments are zero;
- fixed-pinned: reactions are 11P/16 and 5P/16, with fixed-end moment magnitude 3PL/16.

For P = 10 kN and L = 6 m:

- fixed-fixed: 5 / 5 kN reactions, 7.5 kN·m end-moment magnitude;
- pin-roller: 5 / 5 kN reactions, zero end moments;
- fixed-pinned: 6.875 / 3.125 kN reactions, 11.25 kN·m fixed-end moment magnitude.
