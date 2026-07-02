# Regression tests

Golden-file characterization tests for the CardioMechanics binaries. Each test
runs a binary on a fixed input and compares its output to a committed reference
with numerical tolerances (never byte-exact).

Coverage:
- **CellModelTest** — all runnable single-cell ionic models (auto-discovered),
  plus coupled ionic+tension runs for Land17. Serial.
- **CardioMechanics** — benchmark2015 Problem1 (Static solver) at `mpirun -np 4`:
  the `Pressure.dat` time series and the final deformed geometry (last VTU point
  coordinates, via meshio). Marked `mpi`.
- **BidomainMatrixGenerator** — assembles the EM01 mono-domain matrices (serial)
  and compares structural/numeric invariants (dims, nnz, Frobenius norm, sums)
  of the stiffness/mass matrices and material vector, read directly from the
  PETSc binary format (no petsc4py needed).
- **acCELLerate** — standalone EP solve on the EM01 cube at `mpirun -np 4`;
  compares the P8 sensor traces (Vm, Cai). Marked `mpi slow`.
- **CardioMechanics EM01** — full electromechanics (`NewmarkBeta` solver +
  acCELLerate plugin + Land17) at `mpirun -np 4`; compares final deformation and
  the coupled P8 sensor traces. Marked `mpi slow`.

The three EM01 tests share one staged tree and a single matrix-assembly step
(the `em01_root` fixture). The EM01 electromechanics runs are shortened to
`EM01_SIM_LENGTH` (0.05 s) — the full 1 s beat is ~30+ min; 0.05 s still captures
the wavefront reaching the far sensor P8 (~42 ms) and the first clear mechanical
deformation (~30 ms).

## Setup

```sh
python3 -m venv .venv && source .venv/bin/activate
pip install -r tests/requirements.txt
```

## Run

```sh
pytest tests/                 # all
pytest tests/ -k cellmodel    # one binary
pytest -m "not slow"          # skip long integration cases
pytest -m "not mpi"           # skip parallel cases
```

## Regenerate goldens

After an *intended* behaviour change, rebaseline and review the diff before
committing:

```sh
pytest tests/ --update-golden
git diff tests/**/golden
```

## Notes

- **Binary discovery:** set `$CM_BIN_DIR` to pin a build directory, otherwise the
  newest match of `_build/*/bin/**/<name>` is used (release/debug/installed
  agnostic).
- **kaRootDir:** the run fixture sets `kaRootDir` to the repo root in the
  subprocess environment. Without it the binaries segfault at startup.
- **Goldens** were generated on this machine (PETSc 3.24, serial). Cross-environment
  drift is absorbed by tolerances, not per-platform goldens.

## Known non-functional models

The following cell models do not run standalone in the current repo and are
excluded from the CellModelTest suite (`KNOWN_BROKEN` in
`cellmodel/test_cellmodel.py`). These are pre-existing legacy defects, not
regressions — to be investigated separately:

- `HimenoEtAl_Endo`, `HimenoEtAl_Epi`, `HimenoEtAl_Mid` — each tries to open a
  base `HimenoEtAl.ev` that is not shipped in `electrophysiology/data/`.
- `Kurata` — fails parameter initialization ("Init elphy parameters failed").
