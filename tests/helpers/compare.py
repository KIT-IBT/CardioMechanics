import re

import numpy as np

# A header column is `<index>:<name>`; the name may contain spaces (e.g.
# "Ki [M]"), so split on the index markers rather than on whitespace.
_HEADER_COL = re.compile(r"\d+:(.*?)(?=\s+\d+:|\s*$)")


def read_trace(path):
    """Parse a CellModelTest -outfile trace.

    Line 1 is a header of `idx:name` tokens; the remaining lines are
    whitespace-separated floats. Returns (names, data) with data shape
    (rows, cols).
    """
    path = str(path)
    with open(path) as f:
        header = f.readline()
    names = _HEADER_COL.findall(header)
    data = np.loadtxt(path, skiprows=1)
    if data.ndim == 1:
        data = data.reshape(1, -1)
    assert data.shape[1] == len(names), (
        f"{path}: {len(names)} header columns but {data.shape[1]} data columns"
    )
    return names, data


def read_table(path):
    """Parse a whitespace-delimited table whose first line is plain column names.

    Used for outputs like CardioMechanics Pressure.dat
    (`time pressure130 volume130`). Returns (names, data).
    """
    path = str(path)
    with open(path) as f:
        names = f.readline().split()
    data = np.loadtxt(path, skiprows=1)
    if data.ndim == 1:
        data = data.reshape(1, -1)
    assert data.shape[1] == len(names), (
        f"{path}: {len(names)} header columns but {data.shape[1]} data columns"
    )
    return names, data


def read_headerless(path, names):
    """Load a whitespace numeric table that has no header line (e.g. acCELLerate
    sensor files: `t value` per row). Returns (names, data)."""
    data = np.loadtxt(str(path))
    if data.ndim == 1:
        data = data.reshape(1, -1)
    assert data.shape[1] == len(names), (
        f"{path}: expected {len(names)} columns but found {data.shape[1]}"
    )
    return list(names), data


def read_petsc_mat(path):
    """Structural/numeric invariants of a PETSc binary sparse matrix.

    Format (big-endian): int32 classid(1211216), rows, cols, nnz;
    int32 row_lengths[rows]; int32 col_indices[nnz]; float64 values[nnz].
    Returns dict(rows, cols, nnz, fro, vsum) — avoids committing the (large)
    matrix or a PETSc-version-fragile byte comparison.
    """
    with open(str(path), "rb") as f:
        classid, rows, cols, nnz = (int(x) for x in np.fromfile(f, dtype=">i4", count=4))
        assert classid == 1211216, f"{path}: not a PETSc matrix (classid {classid})"
        rowlen = np.fromfile(f, dtype=">i4", count=rows)
        np.fromfile(f, dtype=">i4", count=nnz)  # column indices (skip)
        val = np.fromfile(f, dtype=">f8", count=nnz)
    assert int(rowlen.sum()) == nnz and len(val) == nnz, f"{path}: truncated/inconsistent"
    return {"rows": rows, "cols": cols, "nnz": nnz,
            "fro": float(np.sqrt((val ** 2).sum())), "vsum": float(val.sum())}


def read_petsc_vec(path):
    """Invariants of a PETSc binary vector: dict(n, l2, vsum)."""
    with open(str(path), "rb") as f:
        classid, n = (int(x) for x in np.fromfile(f, dtype=">i4", count=2))
        assert classid == 1211214, f"{path}: not a PETSc vector (classid {classid})"
        val = np.fromfile(f, dtype=">f8", count=n)
    return {"n": n, "l2": float(np.sqrt((val ** 2).sum())), "vsum": float(val.sum())}


def compare_scalars(actual, golden, rtol, atol=0.0):
    """Assert two dicts of scalars match. Integer-valued keys (dims, counts)
    must match exactly; floats compare within tolerance."""
    assert set(actual) == set(golden), f"keys differ: {set(actual)} vs {set(golden)}"
    for k in golden:
        a, g = actual[k], golden[k]
        if isinstance(g, int):
            assert a == g, f"'{k}' changed: actual={a} golden={g}"
        elif not np.isclose(a, g, rtol=rtol, atol=atol):
            raise AssertionError(
                f"'{k}' differs beyond rtol={rtol} atol={atol}: actual={a:.8e} golden={g:.8e}"
            )


def read_vtu_points(path):
    """Deformed point coordinates from a VTU, ordered by PointID.

    CardioMechanics stores no explicit displacement field; the mesh point
    coordinates are the deformed geometry. Sorting by PointID makes the order
    independent of MPI partitioning. Returns (pointid, points) with points
    shape (N, 3). Requires meshio.
    """
    import meshio

    m = meshio.read(str(path))
    pts = m.points
    pid = m.point_data.get("PointID")
    if pid is None:
        return np.arange(len(pts)), pts
    pid = np.asarray(pid).ravel()
    order = np.argsort(pid, kind="stable")
    return pid[order], pts[order]


def write_golden(path, names, data):
    # CellModelTest writes ~7 significant figures; %.8e keeps that faithfully
    # without storing re-expansion noise, and stays well above the compare rtol.
    np.savetxt(str(path), data, fmt="%.8e", delimiter=",", header=",".join(names), comments="")


def read_golden(path):
    path = str(path)
    with open(path) as f:
        names = f.readline().strip().split(",")
    data = np.loadtxt(path, delimiter=",", skiprows=1)
    if data.ndim == 1:
        data = data.reshape(1, -1)
    return names, data


def compare_columns(actual, golden, rtol, atol=0.0):
    """Assert two (names, data) traces match within tolerance, column by column.

    Raises AssertionError naming the offending column, row, and delta.
    """
    a_names, a = actual
    g_names, g = golden
    assert a_names == g_names, f"column names differ:\n actual={a_names}\n golden={g_names}"
    assert a.shape == g.shape, f"shape mismatch: actual={a.shape} golden={g.shape}"
    for j, name in enumerate(g_names):
        if not np.allclose(a[:, j], g[:, j], rtol=rtol, atol=atol, equal_nan=True):
            diff = np.abs(a[:, j] - g[:, j])
            i = int(np.nanargmax(diff))
            raise AssertionError(
                f"column '{name}' differs beyond rtol={rtol} atol={atol}: "
                f"row {i} actual={a[i, j]:.6e} golden={g[i, j]:.6e} "
                f"|delta|={diff[i]:.3e}"
            )
