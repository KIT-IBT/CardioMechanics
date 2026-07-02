from pathlib import Path

import numpy as np
import pytest

from helpers.compare import compare_columns, read_golden, read_vtu_points, write_golden
from helpers.run import run_binary

GOLDEN_DIR = Path(__file__).parent / "golden"
NP = 4
DEFORM_RTOL, DEFORM_ATOL = 1e-4, 1e-8
# Sensor traces are written with 1e-6 absolute precision and the coupled EP
# solve (with mechano-electric feedback) is reproducible only to a few 1e-6
# run-to-run — near solver tolerance. An absolute floor of 2e-5 absorbs that.
SENSOR_RTOL, SENSOR_ATOL = 1e-3, 2e-5

pytestmark = [pytest.mark.mpi, pytest.mark.slow]


@pytest.fixture(scope="module")
def em01_em_out(em01_root, em01_sim_length, cm_env, binary):
    """Run the coupled electromechanics (NewmarkBeta solver + acCELLerate
    plugin, Land17 tension) once on the EM01 cube, shortened to em01_sim_length,
    and return the Results directory."""
    settings = em01_root / "settings"
    xml = (settings / "M_1mm.xml").read_text().replace(
        "<StopTime>1.0</StopTime>", f"<StopTime>{em01_sim_length}</StopTime>"
    )
    (settings / "M_short.xml").write_text(xml)
    run_binary(
        binary("CardioMechanics"),
        ["-settings", "M_short.xml"],
        cwd=settings,
        env=cm_env,
        np=NP,
        timeout=1200,
    )
    return em01_root / "Results"


def _last_vtu(vtu_dir):
    return max(vtu_dir.glob("Cube.*.vtu"), key=lambda p: int(p.stem.split(".")[1]))


def test_em01_deformation(em01_em_out, update_golden):
    pid, pts = read_vtu_points(_last_vtu(em01_em_out / "Cube_vtu"))
    golden_path = GOLDEN_DIR / "em01_em_deformation.npz"
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
            f"node PointID={int(pid[i])} |delta|={d[i]:.3e} m"
        )


def test_em01_coupled_sensors(em01_em_out, update_golden):
    """P8 Vm/Cai from the plugin EP solve (with mechano-electric feedback), a
    distinct signal from the standalone acCELLerate run."""
    from helpers.compare import read_headerless

    ep = em01_em_out / "EP"
    t_vm, vm = read_headerless(ep / "P8_Vm.txt", ["t", "Vm"])[1].T
    t_cai, cai = read_headerless(ep / "P8_Cai.txt", ["t", "Cai"])[1].T
    assert np.allclose(t_vm, t_cai), "Vm and Cai sensor time bases differ"
    actual = (["t", "Vm", "Cai"], np.column_stack([t_vm, vm, cai]))

    golden_path = GOLDEN_DIR / "em01_em_sensors.csv"
    if update_golden:
        GOLDEN_DIR.mkdir(exist_ok=True)
        write_golden(golden_path, *actual)
        pytest.skip(f"updated golden {golden_path.name}")
    assert golden_path.is_file(), f"missing golden {golden_path}; run with --update-golden"
    compare_columns(actual, read_golden(golden_path), rtol=SENSOR_RTOL, atol=SENSOR_ATOL)
