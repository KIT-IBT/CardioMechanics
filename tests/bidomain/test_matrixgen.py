import json
from pathlib import Path

import pytest

from helpers.compare import compare_scalars, read_petsc_mat, read_petsc_vec

GOLDEN = Path(__file__).parent / "golden" / "cube_0.25mm_invariants.json"
RTOL = 1e-9  # matrix assembly is deterministic; this is cross-build margin


def test_matrix_invariants(em01_root, update_golden):
    """BidomainMatrixGenerator output (run once in the em01_root fixture):
    compare structural/numeric invariants of the stiffness matrix, mass matrix
    and material vector against a committed golden."""
    geo = em01_root / "geoFiles"
    actual = {
        "stiffness": read_petsc_mat(geo / "cube_0.25mm.mat"),
        "mass": read_petsc_mat(geo / "cube_0.25mm.mass.mat"),
        "material": read_petsc_vec(geo / "cube_0.25mm.vec"),
    }
    if update_golden:
        GOLDEN.parent.mkdir(exist_ok=True)
        GOLDEN.write_text(json.dumps(actual, indent=2, sort_keys=True))
        pytest.skip(f"updated golden {GOLDEN.name}")

    assert GOLDEN.is_file(), f"missing golden {GOLDEN}; run with --update-golden"
    golden = json.loads(GOLDEN.read_text())
    for key in golden:
        compare_scalars(actual[key], golden[key], rtol=RTOL)
