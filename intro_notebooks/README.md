# PyNite intro notebooks

This directory is a small, CSV-driven learning track for PyNite and basic finite-element modelling.

The intent is deliberately different from the main `notebooks/` directory:

- keep each structural model very small;
- expose the inputs in CSV files;
- separate **geometry**, **supports**, **member releases**, and **loads**;
- compare PyNite results with hand-calculation benchmarks wherever practical;
- use the same vocabulary and data-flow style that later appears in the Z5 tower work.

## Lessons

### 01 — Beam supports and moment releases

`01_beam_supports_and_releases/`

A 6 m straight beam is split at its midpoint and loaded there by 10 kN downward. Five model configurations compare:

1. fixed-fixed supports;
2. pin-roller supports;
3. fixed-left / pinned-right support behaviour;
4. fixed-fixed support nodes with a right-end member moment release;
5. fixed-fixed support nodes with moment releases at both beam ends.

The important distinction is:

> A **support restraint** acts on a node degree of freedom. A **member release** acts on the force/moment transfer between one member end and that node.

That distinction becomes essential in lattice towers because several members can meet at one node and do not necessarily have identical connection behaviour.

## Units

These notebooks use the same explicit unit convention as the current keratia PyNite workflow:

- length: mm
- force: N
- moment: N·mm
- elastic modulus: N/mm² (MPa)

PyNite does not impose a unit system; consistency is the user's responsibility.
