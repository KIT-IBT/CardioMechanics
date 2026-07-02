import os
import shutil
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]

# EM01 short-run milestones (full 1.0s beat is ~30+ min; these capture the
# physics of interest): EP wavefront reaches the far sensor P8 at ~42ms, first
# clear mechanical deformation appears by ~30ms. 0.05s covers both with buffer.
EM01_SIM_LENGTH = 0.05


def pytest_addoption(parser):
    parser.addoption(
        "--update-golden",
        action="store_true",
        default=False,
        help="Regenerate golden references instead of asserting against them.",
    )


@pytest.fixture(scope="session")
def update_golden(request):
    return request.config.getoption("--update-golden")


@pytest.fixture(scope="session")
def repo_root():
    return REPO_ROOT


@pytest.fixture(scope="session")
def cm_env():
    """Environment for every binary invocation.

    kaRootDir must be set: the binaries resolve their bundled data files relative
    to it at runtime and segfault (getenv returns NULL) when it is unset. The
    CMake presets only set kaRootDir at build time, so a plain shell does not
    have it. OMP_NUM_THREADS=1 keeps serial runs deterministic.
    """
    env = os.environ.copy()
    env["kaRootDir"] = str(REPO_ROOT)
    env["OMP_NUM_THREADS"] = "1"
    return env


def _find_binary(name):
    override = os.environ.get("CM_BIN_DIR")
    if override:
        cand = Path(override) / name
        if cand.is_file():
            return cand
        raise FileNotFoundError(f"CM_BIN_DIR={override} set but {name} not found there")
    matches = [
        p
        for p in REPO_ROOT.glob(f"_build/*/bin/**/{name}")
        if p.is_file() and os.access(p, os.X_OK)
    ]
    if not matches:
        raise FileNotFoundError(
            f"{name} not found. Searched $CM_BIN_DIR and "
            f"{REPO_ROOT}/_build/*/bin/**/{name} — build the project first."
        )
    return max(matches, key=lambda p: p.stat().st_mtime)


@pytest.fixture(scope="session")
def binary():
    """Return a resolver: binary('CellModelTest') -> Path to the newest build."""
    return _find_binary


@pytest.fixture(scope="session")
def em01_sim_length():
    return EM01_SIM_LENGTH


@pytest.fixture(scope="session")
def em01_root(binary, cm_env, tmp_path_factory):
    """Stage the EM01 example into an isolated tree and run the FE-matrix
    preprocessing once. Shared by the BidomainMatrixGenerator, acCELLerate and
    full-EM (CardioMechanics M_1mm) tests, which all consume its output.

    Requires mpirun for the downstream parallel runs; the preprocessing itself
    is serial. The staged ep_0.25mm.aclt is shortened to EM01_SIM_LENGTH.
    """
    from helpers.run import run_binary

    if shutil.which("mpirun") is None:
        pytest.skip("mpirun not found")

    src = REPO_ROOT / "examples" / "EM01"
    root = tmp_path_factory.mktemp("em01")
    for sub in ["settings", "tetgen", "cellmodelFiles", "stimFiles", "sensorFiles", "materialFiles"]:
        shutil.copytree(src / sub, root / sub)
    (root / "geoFiles").mkdir()
    shutil.copy(src / "geoFiles" / "cube_0.25mm.vtu", root / "geoFiles")
    (root / "Results" / "EP").mkdir(parents=True)

    # Preprocessing: assemble mono-domain stiffness/mass matrices and material vec.
    run_binary(
        binary("BidomainMatrixGenerator"),
        ["cube_0.25mm.vtu", "-o", "cube_0.25mm", "-mono", "../materialFiles/materialIntra.def"],
        cwd=root / "geoFiles",
        env=cm_env,
        timeout=300,
    )

    # Shorten the EP project length in place (affects both the standalone
    # acCELLerate run and the plugin used by the coupled EM run).
    aclt = root / "settings" / "ep_0.25mm.aclt"
    aclt.write_text(aclt.read_text().replace("CalcLength 1.0", f"CalcLength {EM01_SIM_LENGTH}"))
    return root
