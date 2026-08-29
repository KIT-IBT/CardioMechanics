"""Unit tests for VTK2tetgen's fiber/sheet orthonormalization math.

createONS and NormFiberSheetNormal are pure (numpy only); they build the local
material basis written to the .bases file. A wrong basis silently corrupts the
mechanics, so orthonormality and right-handedness are worth pinning.
"""
import numpy as np
import pytest

import VTK2tetgen as vt


def _unpack(m):
    m = np.asarray(m, dtype=float)
    return m[0:3], m[3:6], m[6:9]


FIBERS = [
    [1, 0, 0], [0, 1, 0], [0, 0, 1],   # axis-aligned: exercises each seed branch
    [1, 1, 0], [2, 3, 6], [-1, 2, -2], [0.1, 0.0, 0.0],
]


@pytest.mark.parametrize("f", FIBERS)
def test_createons_orthonormal_righthanded(f):
    fn, s, sn = _unpack(vt.createONS(f))
    # unit length
    for v in (fn, s, sn):
        assert np.isclose(np.linalg.norm(v), 1.0)
    # mutually orthogonal
    assert np.isclose(np.dot(fn, s), 0.0, atol=1e-12)
    assert np.isclose(np.dot(fn, sn), 0.0, atol=1e-12)
    assert np.isclose(np.dot(s, sn), 0.0, atol=1e-12)
    # fiber direction preserved (just normalized)
    assert np.allclose(fn, np.asarray(f, float) / np.linalg.norm(f))
    # right-handed: f x s == sn and det == +1
    assert np.allclose(np.cross(fn, s), sn, atol=1e-12)
    assert np.isclose(np.linalg.det(np.array([fn, s, sn])), 1.0)


def test_createons_matches_identity_for_x_axis():
    assert np.allclose(vt.createONS([1, 0, 0]), [1, 0, 0, 0, 1, 0, 0, 0, 1])


def test_normfibersheetnormal_normalizes_and_preserves_direction():
    f, s, sn = [2, 0, 0], [0, 3, 0], [0, 0, 4]
    fn, sfn, snn = _unpack(vt.NormFiberSheetNormal(f, s, sn))
    assert np.allclose(fn, [1, 0, 0])
    assert np.allclose(sfn, [0, 1, 0])
    assert np.allclose(snn, [0, 0, 1])


def test_normfibersheetnormal_keeps_nonorthogonal_input():
    # It only normalizes; a non-orthogonal input stays non-orthogonal (no re-basing).
    fn, s, sn = _unpack(vt.NormFiberSheetNormal([1, 1, 0], [0, 1, 0], [0, 0, 1]))
    assert np.isclose(np.linalg.norm(fn), 1.0)
    assert not np.isclose(np.dot(fn, s), 0.0)
