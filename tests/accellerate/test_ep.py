from pathlib import Path

import numpy as np
import pytest

from helpers.compare import compare_columns, read_golden, read_headerless, write_golden
from helpers.run import run_binary

GOLDEN = Path(__file__).parent / "golden" / "em01_sensors.csv"
NP = 4
RTOL, ATOL = 1e-4, 1e-8

pytestmark = pytest.mark.mpi


def _read_sensors(results_ep):
    """Merge the P8 Vm and Cai sensor traces into one (names, data) table."""
    t_vm, vm = read_headerless(results_ep / "P8_Vm.txt", ["t", "Vm"])[1].T
    t_cai, cai = read_headerless(results_ep / "P8_Cai.txt", ["t", "Cai"])[1].T
    assert np.allclose(t_vm, t_cai), "Vm and Cai sensor time bases differ"
    return ["t", "Vm", "Cai"], np.column_stack([t_vm, vm, cai])


@pytest.mark.slow
def test_em01_ep_sensors(em01_root, cm_env, binary, update_golden):
    """Standalone acCELLerate EP solve on the EM01 cube; compare the P8 sensor
    traces (Vm, Cai) that capture the wavefront reaching the far corner."""
    run_binary(
        binary("acCELLerate"),
        ["ep_0.25mm.aclt"],
        cwd=em01_root / "settings",
        env=cm_env,
        np=NP,
        timeout=600,
    )
    actual = _read_sensors(em01_root / "Results" / "EP")

    if update_golden:
        GOLDEN.parent.mkdir(exist_ok=True)
        write_golden(GOLDEN, *actual)
        pytest.skip(f"updated golden {GOLDEN.name}")

    assert GOLDEN.is_file(), f"missing golden {GOLDEN}; run with --update-golden"
    compare_columns(actual, read_golden(GOLDEN), rtol=RTOL, atol=ATOL)
