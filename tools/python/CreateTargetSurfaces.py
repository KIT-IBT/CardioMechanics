#!/usr/bin/env python

"""Generate inverse-problem target surfaces from a CardioMechanics forward run.

For every forward-output VTU the requested surface nodes are extracted (via the
ExtractSurfaceNodesFromVTU tool) and written as a binary `.dat` file; a
time-indexed `.list` file mapping each time to its `.dat` is then emitted for
the PointsCtrl plugin to consume. Times are `index * dt` seconds, sorted
ascending.

The `.dat` binary layout matches what CBDataFromFile reads: a native-endian
int32 value count (= 3 * number of nodes) followed by that many native-endian
float64 coordinates, node by node as x, y, z.
"""

import argparse
import glob
import os
import struct
import subprocess
from pathlib import Path


def vtu_index(path):
    """Timestep index N from a `<prefix>.<N>.vtu` filename."""
    return int(Path(path).name.split(".")[1])


def parse_node_table(text):
    """Parse ExtractSurfaceNodesFromVTU text output into a list of (x, y, z).

    Line 0 is the node count; the remaining non-empty lines are `x y z`.
    """
    lines = [ln for ln in text.splitlines() if ln.strip()]
    return [tuple(map(float, ln.split())) for ln in lines[1:]]


def pack_surface(coords):
    """Pack (x, y, z) rows into the CBDataFromFile binary payload (bytes)."""
    flat = [v for xyz in coords for v in xyz]
    return struct.pack("i", len(flat)) + struct.pack(f"{len(flat)}d", *flat)


def list_entries(indices, dt):
    """Map timestep indices to (time, index) pairs sorted ascending by time."""
    return sorted(((i * dt, i) for i in indices), key=lambda e: e[0])


def create_target_surfaces(results_dir, node_file, sur_file, output_dir, dt,
                           list_path, surfaces="1 2",
                           tool="ExtractSurfaceNodesFromVTU"):
    """Extract, pack and index every VTU in results_dir. Returns the .list text."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    list_path = Path(list_path)

    dats = {}
    for vtu in glob.glob(os.path.join(results_dir, "*.vtu")):
        dat = output_dir / (os.path.basename(vtu) + ".dat")
        subprocess.run([str(tool), vtu, str(node_file), str(sur_file),
                        surfaces, str(dat)], check=True)
        dat.write_bytes(pack_surface(parse_node_table(dat.read_text())))
        dats[vtu_index(vtu)] = dat

    lines = [f"{t} {os.path.relpath(dats[i], list_path.parent)}\n"
             for t, i in list_entries(dats.keys(), dt)]
    text = "".join(lines)
    list_path.write_text(text)
    return text


def parse():
    parser = argparse.ArgumentParser(
        prog="CreateTargetSurfaces",
        description="Generate inverse-problem target surfaces from a forward run.")
    parser.add_argument("-results", required=True,
                        help="Directory of forward-output VTU files.")
    parser.add_argument("-node", required=True, help="Mesh .node file.")
    parser.add_argument("-sur", required=True, help="Mesh .sur file.")
    parser.add_argument("-output", required=True,
                        help="Output directory for the target .dat files.")
    parser.add_argument("-list", required=True, help="Output .list file path.")
    parser.add_argument("-dt", type=float, required=True,
                        help="Time between exported VTUs in seconds.")
    parser.add_argument("-surfaces", default="1 2",
                        help="Surface indices to extract (default: '1 2').")
    parser.add_argument("-tool", default="ExtractSurfaceNodesFromVTU",
                        help="Path to the ExtractSurfaceNodesFromVTU binary.")
    return parser.parse_args()


if __name__ == "__main__":
    args = parse()
    create_target_surfaces(args.results, args.node, args.sur, args.output,
                           args.dt, args.list, args.surfaces, args.tool)
