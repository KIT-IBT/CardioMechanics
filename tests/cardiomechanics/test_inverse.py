"""Inverse-problem (active-stress estimator) regression test.

Full round-trip on the ellipsoid example: a forward run applies a known
DoubleHill active tension and deforms the mesh; its surface is extracted as the
target; the estimator then recovers the per-element active stress that
reproduces that target. Both runs are shortened to StopTime=0.25 (through the
activation onset) to keep the test tractable; the estimator is direct/serial to
match the golden. We freeze our own recovered ActiveStress and deformation as
the golden (there is no external reference: the old code's target-generating
ActiveStress path was intentionally removed in the modern fork).
"""
import shutil
import sys
from pathlib import Path

import numpy as np
import pytest

from helpers.compare import read_vtu_cell_field, read_vtu_points
from helpers.run import run_binary

REPO_ROOT = Path(__file__).resolve().parents[2]
SRC = REPO_ROOT / "examples" / "inverseEllipsoid"
GOLDEN_DIR = Path(__file__).parent / "golden"
TOOLS_PYTHON = REPO_ROOT / "tools" / "python"

SIM_LENGTH = 0.25          # forward + inverse StopTime; export dt is 1e-2
LAST = 25                  # Inverse.<LAST>.vtu at t=0.25
STRESS_RTOL, STRESS_ATOL = 1e-4, 1e-2      # ActiveStress in Pa
DEFORM_RTOL, DEFORM_ATOL = 1e-4, 1e-8      # coordinates in m


@pytest.fixture(scope="module")
def inverse_vtu_dir(binary, cm_env, tmp_path_factory):
    """Stage the example, run forward -> target extraction -> inverse, serially."""
    pytest.importorskip("meshio")
    wd = tmp_path_factory.mktemp("inverse")
    shutil.copytree(SRC / "geometry", wd / "geometry")
    for xml in ("forwardT4.xml", "inverseT4.xml"):
        text = (SRC / xml).read_text().replace(
            "<StopTime>1</StopTime>", f"<StopTime>{SIM_LENGTH}</StopTime>")
        (wd / xml).write_text(text)
    (wd / "Results").mkdir()

    run_binary(binary("CardioMechanics"), ["-settings", "forwardT4.xml"],
               cwd=wd, env=cm_env, timeout=300)

    sys.path.insert(0, str(TOOLS_PYTHON))
    import CreateTargetSurfaces as cts
    cts.create_target_surfaces(
        wd / "Results" / "forward_vtu",
        wd / "geometry" / "meshT4IP.node", wd / "geometry" / "meshT4IP.sur",
        wd / "TargetSurfaces", dt=0.01, list_path=wd / "TargetSurfaces.list",
        tool=binary("ExtractSurfaceNodesFromVTU"))

    run_binary(binary("CardioMechanics"), ["-settings", "inverseT4.xml"],
               cwd=wd, env=cm_env, timeout=1200)
    return wd / "Results" / "Inverse_vtu"


def test_inverse_active_stress(inverse_vtu_dir, update_golden):
    cid, stress = read_vtu_cell_field(inverse_vtu_dir / f"Inverse.{LAST}.vtu", "ActiveStress")
    golden_path = GOLDEN_DIR / "inverse_activestress.npz"
    if update_golden:
        GOLDEN_DIR.mkdir(exist_ok=True)
        np.savez_compressed(golden_path, cellid=cid, stress=stress)
        pytest.skip(f"updated golden {golden_path.name}")
    assert golden_path.is_file(), f"missing golden {golden_path}; run with --update-golden"
    g = np.load(golden_path)
    assert np.array_equal(cid, g["cellid"]), "cell ordering / mesh identity changed"
    if not np.allclose(stress, g["stress"], rtol=STRESS_RTOL, atol=STRESS_ATOL):
        d = np.abs(stress - g["stress"])
        i = int(np.argmax(d))
        raise AssertionError(
            f"recovered ActiveStress differs beyond rtol={STRESS_RTOL} atol={STRESS_ATOL}: "
            f"cell CellID={int(cid[i])} actual={stress[i]:.6e} golden={g['stress'][i]:.6e} "
            f"|delta|={d[i]:.3e} Pa")


def test_inverse_deformation(inverse_vtu_dir, update_golden):
    pid, pts = read_vtu_points(inverse_vtu_dir / f"Inverse.{LAST}.vtu")
    golden_path = GOLDEN_DIR / "inverse_deformation.npz"
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
            f"|delta|={d[i]:.3e} m")
