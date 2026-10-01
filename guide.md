# Guide to `Z5_07_frame_end_releases.csv`

This file controls **frame-member end releases** for:

```text
notebooks/07_z5_pynite_release_sensitivity.ipynb
```

It is intended for sensitivity studies of joint behaviour in the Z5 PyNite model.

The committed CSV reproduces the notebook-05 frame assumptions:

- every frame end is `rigid`;
- no user-defined frame-end releases are active;
- truss members are **not** listed here because their pin-like PyNite treatment is applied automatically by the adapter.

The central modelling idea is:

> A release belongs to a **specific end of a specific member**, not to a node as a whole.

That distinction lets a continuous physical member remain continuous through an intermediate node while another physical member meeting at the same node can be pin-like.

---

## 1. CSV columns

The columns are:

```text
element_id
start_node
start_label
end_node
end_label
parent_element_id
i_end_type
j_end_type
custom_i_Dx
custom_i_Dy
custom_i_Dz
custom_i_Rx
custom_i_Ry
custom_i_Rz
custom_j_Dx
custom_j_Dy
custom_j_Dz
custom_j_Rx
custom_j_Ry
custom_j_Rz
notes
```

### `element_id`

The finite-element ID used by the Z5 model.

Example:

```text
138
```

means PyNite member:

```text
E138
```

Do not change this ID unless the structural input model itself changes.

---

### `start_node` and `end_node`

These identify the two nodes of the element.

The notebook treats:

- `start_node` as the member's **i-end**;
- `end_node` as the member's **j-end**.

Therefore:

```text
element_id = 138
start_node = 3
end_node = 29
```

means:

```text
i-end = node 3
j-end = node 29
```

This direction matters because releases are specified independently for the i and j ends.

---

### `start_label` and `end_label`

Human-readable versions of the two node IDs.

They are included to make the CSV easier to inspect.

For example:

```text
3, top_xneg_ypos
29, z5_e10_mid
```

The notebook uses the numerical element/node definitions for the model. These labels should normally be left unchanged.

---

### `parent_element_id`

The ID of the **physical parent member** from which a finite-element segment came.

This is especially important when one physical angle has been subdivided so braces can connect at intermediate nodes.

For example:

```text
element 138: node 3 -> node 29, parent_element_id = 10
element 139: node 29 -> node 4, parent_element_id = 10
```

Elements 138 and 139 are therefore two FE segments of one physical parent member:

```text
node 3 -------- node 29 -------- node 4
          138             139
             parent 10
```

If parent 10 is physically continuous through node 29, the normal assumption is:

```text
138 j-end = rigid
139 i-end = rigid
```

A brace merely attaching at node 29 is **not** by itself a reason to introduce a hinge into parent 10.

---

## 2. End-type fields

The main editable fields are:

```text
i_end_type
j_end_type
```

Each accepts one of:

```text
rigid
bending_pin
custom
```

### `rigid`

No member-end DOF is released.

Equivalent release flags:

```text
Dx = False
Dy = False
Dz = False
Rx = False
Ry = False
Rz = False
```

This is the current notebook-05 baseline for every frame member.

A rigid end can transmit:

- axial force;
- local shear forces;
- torsion;
- bending moments.

Within the FE idealization, the member end shares the relevant nodal translations and rotations without an end release.

---

### `bending_pin`

This is the convenient preset for a pin-like bending connection.

It is interpreted as:

```text
Dx = False
Dy = False
Dz = False
Rx = False
Ry = True
Rz = True
```

Therefore:

- translations remain connected;
- local torsional rotation `Rx` remains connected;
- local bending rotations `Ry` and `Rz` are released.

Conceptually:

```text
My = 0
Mz = 0
```

at that idealized member end, within numerical tolerance.

This is the preset intended for the first gusset/pin sensitivity studies.

It deliberately **does not release torsion**. That keeps the first sensitivity experiment narrower and reduces the chance of introducing unnecessary mechanisms.

---

### `custom`

When an end is set to `custom`, the notebook reads the six corresponding Boolean fields.

For the i-end:

```text
custom_i_Dx
custom_i_Dy
custom_i_Dz
custom_i_Rx
custom_i_Ry
custom_i_Rz
```

For the j-end:

```text
custom_j_Dx
custom_j_Dy
custom_j_Dz
custom_j_Rx
custom_j_Ry
custom_j_Rz
```

For these fields:

```text
False = DOF remains connected
True  = DOF is released
```

So, for example:

```text
custom_i_Ry = True
custom_i_Rz = False
```

releases only local i-end rotation about local y.

### Important

The custom Boolean fields are only used when the corresponding end type is:

```text
custom
```

If:

```text
i_end_type = bending_pin
```

then the notebook applies the `bending_pin` preset regardless of the values currently written in `custom_i_*`.

Likewise, if:

```text
i_end_type = rigid
```

the custom i-end fields are ignored.

---

## 3. Meaning of the six custom DOFs

These are **member-local** PyNite end DOFs.

They are not the same thing as the global X/Y/Z directions of the tower.

### Translational DOFs

```text
Dx
Dy
Dz
```

- `Dx`: local translation along the member x-axis;
- `Dy`: local translation along the member y-axis;
- `Dz`: local translation along the member z-axis.

For an ordinary frame connection, these should normally remain:

```text
False
```

because releasing translations changes force transfer much more radically than a simple moment release.

For example:

```text
Dx = True
```

allows that member end to slide axially relative to the node.

That is **not** an ordinary pin joint.

---

### Rotational DOFs

```text
Rx
Ry
Rz
```

- `Rx`: rotation about the member's local longitudinal x-axis; associated primarily with torsion;
- `Ry`: rotation about local y; associated with one bending axis;
- `Rz`: rotation about local z; associated with the other bending axis.

The `bending_pin` preset is therefore:

```text
Rx = False
Ry = True
Rz = True
```

---

## 4. Local axes are important

The release fields refer to **local element axes**.

For a member running in an arbitrary 3D direction:

```text
local x = along the member from i to j
local y = PyNite member local y
local z = PyNite member local z
```

Therefore:

```text
custom_i_Ry = True
```

does **not** mean:

> release rotation about global Y.

It means:

> release rotation about member E###'s local y-axis at its i-end.

For the initial sensitivity studies, using `bending_pin` avoids having to choose only one local bending axis.

---

## 5. `notes`

Free-form human-readable comments.

The notebook does not use this field mechanically.

Examples:

```text
Notebook-05 baseline: no frame-end release
Physical gusset candidate
Sensitivity test only
Keep continuous with parent 14
```

Use this field to document why a row was changed.

---

# 6. Worked examples with arbitrary geometry

The following examples are schematic and are **not instructions to modify those exact Z5 elements**. They illustrate how the CSV behaves.

---

## Example A — two independent frame members meeting at a pin-like gusset

Suppose:

```text
A ----- E200 ----- B ----- E201 ----- C
```

and E200 and E201 are physically separate steel members meeting at a gusset at node B.

Assume:

```text
E200: start_node=A, end_node=B
E201: start_node=B, end_node=C
```

Then B is:

- the **j-end** of E200;
- the **i-end** of E201.

If both physical members are intended to be bending-pin connected to the gusset at B:

```text
element_id,i_end_type,j_end_type
200,rigid,bending_pin
201,bending_pin,rigid
```

Conceptually:

```text
A ======= E200 ======= o B o ======= E201 ======= C
                       ^   ^
                    released bending
                    at each member end
```

This removes bending-moment transfer from E200 into the joint and from the joint into E201.

---

## Example B — one continuous physical member split for FE connectivity

Suppose one physical angle runs continuously from A to C but has been split at B because a brace connects there:

```text
A ===== E300 ===== B ===== E301 ===== C
                       |
                       |
                    brace
```

and:

```text
E300 parent_element_id = 50
E301 parent_element_id = 50
```

If the angle is physically continuous through B, do **not** create a hinge merely because B exists as an FE node.

Use:

```text
E300 j_end_type = rigid
E301 i_end_type = rigid
```

The brace can still be axial/pin-like independently.

This situation is directly analogous to many of the Z5 split parent members.

---

## Example C — continuous chord plus a separate frame member at the same node

Suppose:

```text
A ===== E400 ===== B ===== E401 ===== C
                       \
                        \
                         E402
                          \
                           D
```

E400 and E401 are two segments of the same physical chord.

E402 is a separate physical frame member attached at B.

Then a plausible idealization is:

```text
E400 j-end = rigid
E401 i-end = rigid
E402 i-end = bending_pin
```

assuming E402 is defined B -> D.

This is exactly why releases are stored **per member end** instead of by node.

Node B can simultaneously represent:

- continuity of one physical member;
- a pin-like end of another physical member.

---

## Example D — one-end bending release

Suppose E500 is:

```text
node 70 ---- E500 ---- node 71
```

and only the node-71 connection is intended to be pin-like.

If:

```text
start_node = 70
end_node   = 71
```

then:

```text
i_end_type = rigid
j_end_type = bending_pin
```

No custom fields need to be changed.

---

## Example E — custom one-axis moment release

Suppose an experimental connection should release bending about local z but retain bending about local y and retain torsion.

Set:

```text
i_end_type = custom
```

with:

```text
custom_i_Dx = False
custom_i_Dy = False
custom_i_Dz = False
custom_i_Rx = False
custom_i_Ry = False
custom_i_Rz = True
```

That is more specialized than the standard `bending_pin` preset.

Use it only when the connection idealization provides a reason to distinguish the two local bending axes.

---

## Example F — full rotational release

A mathematically more complete rotational pin at one end could be represented as:

```text
i_end_type = custom

custom_i_Dx = False
custom_i_Dy = False
custom_i_Dz = False

custom_i_Rx = True
custom_i_Ry = True
custom_i_Rz = True
```

This releases torsion as well as both bending rotations.

This is **not** the current `bending_pin` preset.

Releasing `Rx` broadly can introduce additional zero-stiffness rotational modes or otherwise make the model less robust, so it should be tested separately.

---

# 7. Example using the existing Z5 split parent members

Consider:

```text
element 138: node 3  -> node 29
element 139: node 29 -> node 4
parent_element_id = 10 for both
```

The geometry is:

```text
node 3 ===== E138 ===== node 29 ===== E139 ===== node 4
                         parent 10
```

In the committed 07 baseline:

```text
E138: rigid, rigid
E139: rigid, rigid
```

This preserves bending continuity through node 29.

Changing:

```text
E138 j-end = bending_pin
E139 i-end = bending_pin
```

would intentionally insert a hinge into physical parent 10 at node 29.

That is usually **not** what is intended merely because braces also attach to node 29.

By contrast, if a separate physical frame member terminated at node 29, releasing that separate member end could be reasonable while keeping E138/E139 continuous.

---

# 8. Example using a Z5 physical joint

Consider element 6:

```text
element 6
start_node = 0  (center)
end_node   = 3  (top_xneg_ypos)
```

If a sensitivity study assumes element 6 is bending-pin connected at node 3 but remains rigid at node 0:

```text
element_id = 6
i_end_type = rigid
j_end_type = bending_pin
```

If instead both physical ends are assumed pin-like:

```text
i_end_type = bending_pin
j_end_type = bending_pin
```

This change affects only element 6.

It does **not** automatically release every other frame member meeting nodes 0 or 3.

Each of those member ends must be specified independently.

---

# 9. Recommended editing workflow

For a controlled sensitivity study:

1. Start from the committed baseline CSV.
2. Change only the specific frame end or ends being investigated.
3. Add a short explanation in `notes`.
4. Run notebook 07.
5. Check:
   - whether PyNite reports instability;
   - deformed shape;
   - maximum translations;
   - axial-force redistribution;
   - bending moments at released ends;
   - reactions;
   - whether unexpected local mechanisms appear.
6. Compare against the unchanged notebook-05-equivalent baseline.
7. Commit the configuration separately if it represents a meaningful modelling case.

Avoid changing many unrelated joints at once. A staged sensitivity study makes it much easier to identify which modelling assumption caused a change.

---

# 10. Practical presets

## Current 05-equivalent baseline

```text
i_end_type = rigid
j_end_type = rigid
```

for every frame element.

---

## Bending pin at i-end

```text
i_end_type = bending_pin
j_end_type = rigid
```

Equivalent i-end flags:

```text
Dx=False
Dy=False
Dz=False
Rx=False
Ry=True
Rz=True
```

---

## Bending pin at j-end

```text
i_end_type = rigid
j_end_type = bending_pin
```

---

## Bending pin at both ends

```text
i_end_type = bending_pin
j_end_type = bending_pin
```

---

## Custom release

Set the relevant end to:

```text
custom
```

and edit only that end's six `custom_*` fields.

---

# 11. Important modelling cautions

### A release is not the same as a support

`Z5_07_frame_end_releases.csv` controls **member-to-node behaviour**.

`Z5_07_supports.csv` controls **node-to-ground restraints**.

For example:

- changing node 7 from `fixed` to `pinned` is a support change;
- setting E142's i-end to `bending_pin` is a member-end connection change.

Those are different operations.

### Do not release both sides of a computational split without a physical reason

If E138 and E139 are one continuous physical angle, releasing:

```text
E138 j-end
E139 i-end
```

creates an actual hinge at the FE subdivision point.

That can remove real bending stiffness and can create mechanisms, as the notebook-06 experiment demonstrated more generally.

### Translational releases are high-impact

`Dx/Dy/Dz=True` permits relative translation between the member end and node in that local direction.

That is generally much more drastic than a moment release and should not be used as shorthand for a pin.

### Release directions are local

`Ry` and `Rz` refer to PyNite's **member-local** bending axes.

They are not global Y/Z tower rotations.

### A stable numerical result is not by itself proof of a correct connection model

The purpose of notebook 07 is sensitivity analysis. Release choices should ultimately be tied to the actual steel detail, gusset arrangement, continuity of the angle, bolts/welds, and the intended structural idealization.
