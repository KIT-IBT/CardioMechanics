import shutil
from pathlib import Path

import numpy as np
import pytest

from helpers.compare import compare_columns, read_golden, read_table, read_vtu_points, write_golden
from helpers.run import run_binary

REPO_ROOT = Path(__file__).resolve().parents[2]
SRC = REPO_ROOT / "examples" / "benchmark2015"
GOLDEN_DIR = Path(__file__).parent / "golden"

NP = 4  # CardioMechanics is tested only in parallel (see REGRESSION_PLAN §15.3)
PRESSURE_RTOL, PRESSURE_ATOL = 1e-4, 1e-6
DEFORM_RTOL, DEFORM_ATOL = 1e-4, 1e-8

pytestmark = pytest.mark.mpi


@pytest.fixture(scope="module")
def problem1_out(binary, cm_env, tmp_path_factory):
    """Run Problem1 (Static solver) once at np=NP; return its output directory.

    Problem1.xml uses paths relative to the working dir, so inputs are staged
    into an isolated tmp dir that also holds the ./Problem1 output folder.
    """
    if shutil.which("mpirun") is None:
        pytest.skip("mpirun not found")
    wd = tmp_path_factory.mktemp("problem1")
    shutil.copy(SRC / "Problem1.xml", wd)
    (wd / "tetgen").mkdir()
    for f in SRC.glob("tetgen/bar.*"):
        shutil.copy(f, wd / "tetgen")
    (wd / "Problem1").mkdir()
    run_binary(
        binary("CardioMechanics"),
        ["-settings", "Problem1.xml"],
        cwd=wd,
        env=cm_env,
        np=NP,
        timeout=600,
    )
    return wd / "Problem1"


def test_problem1_pressure(problem1_out, update_golden):
    actual = read_table(problem1_out / "Pressure.dat")
    golden_path = GOLDEN_DIR / "problem1_pressure.csv"
    if update_golden:
        GOLDEN_DIR.mkdir(exist_ok=True)
        write_golden(golden_path, *actual)
        pytest.skip(f"updated golden {golden_path.name}")
    assert golden_path.is_file(), f"missing golden {golden_path}; run with --update-golden"
    compare_columns(actual, read_golden(golden_path), rtol=PRESSURE_RTOL, atol=PRESSURE_ATOL)


def _last_vtu(vtu_dir):
    return max(vtu_dir.glob("Bar.*.vtu"), key=lambda p: int(p.stem.split(".")[1]))


def test_problem1_deformation(problem1_out, update_golden):
    pytest.importorskip("meshio")
    pid, pts = read_vtu_points(_last_vtu(problem1_out / "Bar_vtu"))
    golden_path = GOLDEN_DIR / "problem1_deformation.npz"
    if update_golden:
        GOLDEN_DIR.mkdir(exist_ok=True)
        np.savez_compressed(golden_path, pointid=pid, points=pts)
        pytest.skip(f"updated golden {golden_path.name}")
    assert golden_path.is_file(), f"missing golden {golden_path}; run with --update-golden"
    g = np.load(golden_path)
    assert np.array_equal(pid, g["pointid"]), "point ordering / mesh identity changed"
    if not np.allclose(pts, g["points"], rtol=DEFORM_RTOL, atol=DEFORM_ATOL):
        d = np.linalg.norm(pts - g["points"], axis=1)
        i = int(np.argmax(d))
        raise AssertionError(
            f"deformed coordinates differ beyond rtol={DEFORM_RTOL} atol={DEFORM_ATOL}: "
            f"node PointID={int(pid[i])} actual={pts[i]} golden={g['points'][i]} "
            f"|delta|={d[i]:.3e} m"
        )
