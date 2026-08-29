"""Unit tests for CreateTargetSurfaces' pure transform logic.

The tool extracts surface node coordinates per forward timestep and writes them
as the binary .dat format CBDataFromFile reads, plus a time-indexed .list. The
extraction itself is an external binary; here we pin the deterministic glue:
timestep indexing, text parsing, binary packing, and time normalization.
"""
import struct

import pytest

import CreateTargetSurfaces as cts


@pytest.mark.parametrize("name,expected", [
    ("Inverse.25.vtu", 25),
    ("forward.0.vtu", 0),
    ("/some/dir/forward.7.vtu", 7),
    ("prefix.123.vtu", 123),
])
def test_vtu_index(name, expected):
    assert cts.vtu_index(name) == expected


def test_parse_node_table_skips_count_and_blanks():
    text = "3\n1 2 3\n4 5 6\n\n7 8 9\n"
    assert cts.parse_node_table(text) == [(1.0, 2.0, 3.0), (4.0, 5.0, 6.0), (7.0, 8.0, 9.0)]


def test_pack_surface_layout_and_roundtrip():
    coords = [(1.0, 2.0, 3.0), (4.0, 5.0, 6.0)]
    blob = cts.pack_surface(coords)
    count = struct.unpack_from("i", blob, 0)[0]
    assert count == 6  # 3 * number of nodes
    values = struct.unpack_from("6d", blob, struct.calcsize("i"))
    assert list(values) == [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
    assert len(blob) == struct.calcsize("i") + 6 * struct.calcsize("d")


def test_pack_surface_empty():
    blob = cts.pack_surface([])
    assert struct.unpack_from("i", blob, 0)[0] == 0
    assert len(blob) == struct.calcsize("i")


def test_list_entries_sorted_by_time():
    assert cts.list_entries([2, 0, 1], 0.01) == [(0.0, 0), (0.01, 1), (0.02, 2)]


def test_list_entries_uses_dt_scaling():
    (t0, i0), (t1, i1) = cts.list_entries([0, 10], 0.05)
    assert (t0, i0) == (0.0, 0)
    assert t1 == pytest.approx(0.5) and i1 == 10
