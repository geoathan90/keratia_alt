# PyNite comparison workflow

This project now includes an experimental PyNite backend for comparing the in-house educational stiffness solver against an established finite-element library.

## Install

From the repository root in your Codespace:

```bash
python -m pip install --user -r requirements.txt
```

The PyPI package is named:

```text
PyniteFEA
```

The Python import is:

```python
from Pynite import FEModel3D
```

Do **not** install the unrelated package named `Pynite`.

## New files

```text
pynite_adapter.py
notebooks/04_pynite_comparison_workflow.ipynb
```

## Recommended first comparison settings

```python
SECTION_MODE = "isotropic"
CROSSING_DIAGONAL_MODE = "continuous"
```

Why isotropic first?

The in-house solver and PyNite do not construct asymmetric L-angle local axes in exactly the same way. Using `Iy = Iz = Iavg` removes most local-axis-orientation sensitivity, so the first comparison focuses more on solver/model behaviour than angle-axis convention.

Why continuous first?

The PyNite adapter does not currently reproduce the in-house solver's tiny directional stabilization springs for connected truss-only crossing nodes. The `continuous` crossing assumption is therefore the cleaner starting point.

## Truss representation caveat

The in-house solver has true axial-only truss elements.

In the PyNite adapter, CSV `truss` elements are represented as pin-ended PyNite members by releasing local `Ry/Rz` rotations at both ends. This is a useful comparison model, but it is not mathematically identical to the in-house axial-only truss element.

Treat the PyNite notebook as a second analysis route/check, not as a strict replacement.
