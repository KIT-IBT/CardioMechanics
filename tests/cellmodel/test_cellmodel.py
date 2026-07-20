from pathlib import Path

import pytest

from helpers.compare import compare_columns, read_golden, read_trace, write_golden
from helpers.run import run_binary

REPO_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = REPO_ROOT / "electrophysiology" / "data"
GOLDEN_DIR = Path(__file__).parent / "golden"

TEND = "1.0"
TINC = "1e-5"
TOUTINC = "500"  # output every 500 steps -> ~200 sampled rows
RTOL = 1e-6
ATOL = 1e-9

# Models that do not run standalone in the current repo (legacy defects, not
# regressions): HimenoEtAl_* expect a base HimenoEtAl.ev that is not shipped;
# Kurata fails parameter initialization.
KNOWN_BROKEN = {"HimenoEtAl_Endo", "HimenoEtAl_Epi", "HimenoEtAl_Mid", "Kurata"}


def _discover_ev_models():
    # Non-symlink .ev files only: the symlinks (e.g. OHaraRudy -> OHaraRudy_endo)
    # would just duplicate their targets' traces.
    stems = sorted(p.stem for p in DATA_DIR.glob("*.ev") if not p.is_symlink())
    return [s for s in stems if s not in KNOWN_BROKEN]


# Coupled ionic+tension runs (-evfile X.ev -fvfile Y.fv in one process): the
# cell's calcium transient forward-drives the tension model.
COUPLED = [
    ("TenTusscher2_original_Endo", "Land17"),
    ("OHaraRudy_endo", "Land17"),
    ("CourtemancheEtAl", "Land17Atrial"),
]

CASES = [(m, m, None) for m in _discover_ev_models()]
CASES += [(f"{ev}+{fv}", ev, fv) for ev, fv in COUPLED]


@pytest.mark.parametrize("case_id,ev,fv", CASES, ids=[c[0] for c in CASES])
def test_cellmodel_trace(case_id, ev, fv, binary, cm_env, tmp_path, update_golden):
    evfile = DATA_DIR / f"{ev}.ev"
    assert evfile.is_file(), f"missing ev file {evfile}"
    args = ["-evfile", evfile, "-tend", TEND, "-tinc", TINC, "-toutinc", TOUTINC]
    if fv is not None:
        fvfile = DATA_DIR / f"{fv}.fv"
        assert fvfile.is_file(), f"missing fv file {fvfile}"
        args += ["-fvfile", fvfile]
    outfile = tmp_path / f"{case_id}.txt"
    args += ["-outfile", outfile]

    run_binary(binary("CellModelTest"), args, cwd=tmp_path, env=cm_env)
    actual = read_trace(outfile)

    golden_path = GOLDEN_DIR / f"{case_id}.ref.csv"
    if update_golden:
        GOLDEN_DIR.mkdir(exist_ok=True)
        write_golden(golden_path, *actual)
        pytest.skip(f"updated golden {golden_path.name}")

    assert golden_path.is_file(), f"missing golden {golden_path}; run with --update-golden"
    compare_columns(actual, read_golden(golden_path), rtol=RTOL, atol=ATOL)
