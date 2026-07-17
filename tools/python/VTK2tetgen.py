#!/usr/bin/env python

'''
convert a vtk file from gmsh to .node and .ele files for CardioMechanics
Lukas Baron, Wed Aug 10 2016
Features:
- .node takes fixation from vtk "Fixation" field
- .ele evaluates vtk "Material" field
- .ele is T4 or T10 depending on the vtk Tetrahedron type
- .bases uses vtk "Fiber" field
- .sur has "0 0" as vtk/surface numbers

Tobias Gerach, Tue Dec 11 2018
Merged the two versions of VTK2tetgen.py. Specifically, I moved the fixMaterial option from Robin Andlauer's version
into this one and added proper argument parsing to the script.
This is now the one and only version...the other one was removed.

Tobias Gerach, Tue Apr 9 2019
Implemented that fiber, sheet, and sheet normal direction are used from the vtk if available.
If only the fiber direction is given, sheet and sheet normal directions are set according to scalar and cross product rules.
Fiber, sheet, and sheet normal vectors will all be normalized correctly.
The .bases file will now be written in T10 format if VTU is a T10 mesh.
'''

import argparse
import os

import numpy as np

# CardioMechanics ids start from 1, vtk ids start from 0.
START_FROM_ONE = True


def nodeIsFixed(data, array, point, fixMaterial):
    """Return 7 if any cell incident to `point` has a material in `fixMaterial`."""
    import vtk

    pointCellList = vtk.vtkIdList()
    data.GetPointCells(point, pointCellList)
    for i in range(pointCellList.GetNumberOfIds()):
        if array.GetValue(pointCellList.GetId(i)) in fixMaterial:
            return 7
    return 0


def createONS(f):
    """Orthonormal fiber/sheet/sheet-normal basis from a fiber vector alone.

    The sheet direction is seeded from the Cartesian axis least aligned with the
    fiber, then made orthonormal to it; the sheet normal completes a right-handed
    system. Returns the nine components [f, s, sn].
    """
    f = f / np.linalg.norm(f)
    a, b, c = abs(np.dot(f, [1, 0, 0])), abs(np.dot(f, [0, 1, 0])), abs(np.dot(f, [0, 0, 1]))
    if a <= b and a <= c:
        s = np.array([1, 0, 0])
    elif b <= a and b <= c:
        s = np.array([0, 1, 0])
    else:
        s = np.array([0, 0, 1])
    sn = np.cross(f, s)
    sn = sn / np.linalg.norm(sn)
    s = np.cross(sn, f)
    s = s / np.linalg.norm(s)
    return [f[0], f[1], f[2], s[0], s[1], s[2], sn[0], sn[1], sn[2]]


def NormFiberSheetNormal(f, s, sn):
    """Normalize a given fiber/sheet/sheet-normal triple. Returns [f, s, sn]."""
    f = f / np.linalg.norm(f)
    s = s / np.linalg.norm(s)
    sn = sn / np.linalg.norm(sn)
    return [f[0], f[1], f[2], s[0], s[1], s[2], sn[0], sn[1], sn[2]]


def parse():
    parser = argparse.ArgumentParser(
        description='Convert a vtk file to .node and .ele files for CardioMechanics.')
    parser.add_argument('mesh', help='VTK input file.')
    parser.add_argument('-outfile', type=str, help='Set an alternative output filename PREFIX if desired.')
    parser.add_argument('-scale', type=int, help='Scale the VTK input file by a factor SCALE.')
    parser.add_argument('-fixMaterial', nargs='+', type=int,
                        help='Fix nodes if they are connected to cells that match the given mask in "Material". Example: "-fixMaterial 32 162"')
    return parser.parse_args()


def _write(path, lines):
    """Write tetgen-style output: the trailing summary line first, then the rows."""
    with open(path, "w") as fh:
        fh.write(lines.pop() + "\n")
        for line in lines:
            fh.write(line + "\n")


def main():
    import vtk

    args = parse()
    filename = args.mesh
    unitScaling = float(args.scale) if args.scale else 1

    if filename.lower().endswith('vtu'):  # CardioMechanics output files
        reader = vtk.vtkXMLUnstructuredGridReader()
    else:
        reader = vtk.vtkUnstructuredGridReader()
    reader.SetFileName(filename)
    reader.Update()
    data = reader.GetOutput()

    matArray = data.GetCellData().GetArray('Material') if args.fixMaterial else None

    hasFixation = data.GetPointData().GetArray("Fixation") is not None
    if hasFixation:
        print("Node fixation found")
    hasMaterial = data.GetCellData().GetArray("Material") is not None
    if hasMaterial:
        print("Material information found")

    print("Nodes:", data.GetNumberOfPoints())
    nodeLines = []
    for pid in range(data.GetNumberOfPoints()):
        coords = unitScaling * np.array(data.GetPoint(pid))
        fixation = 0
        if args.fixMaterial:
            fixation = nodeIsFixed(data, matArray, pid, args.fixMaterial)
        if hasFixation:
            fixation = int(data.GetPointData().GetArray("Fixation").GetTuple1(pid))
        idx = pid + 1 if START_FROM_ONE else pid
        nodeLines.append(f"{idx} {coords[0]} {coords[1]} {coords[2]} {fixation}")
    nodeLines.append(f"{data.GetNumberOfPoints()} 3 1 0")

    hasT10 = False
    print("Elements:", data.GetNumberOfCells())
    eleLines = []
    cid = 0
    for vtkcid in range(data.GetNumberOfCells()):
        cell = data.GetCell(vtkcid)
        pids = cell.GetPointIds()
        npts = cell.GetNumberOfPoints()
        if npts == 10 and not hasT10:
            hasT10 = True
            print("  T10 elements found")
        if npts in (4, 10):
            idx = cid + 1 if START_FROM_ONE else cid
            parts = [str(idx)]
            for pid in range(pids.GetNumberOfIds()):
                parts.append(str(pids.GetId(pid) + 1 if START_FROM_ONE else pids.GetId(pid)))
            mat = int(data.GetCellData().GetArray("Material").GetTuple1(vtkcid)) if hasMaterial else 0
            parts.append(str(mat))
            eleLines.append(" ".join(parts))
            cid += 1
    eleLines.append(f"{cid} 10 1" if hasT10 else f"{cid} 4 1")

    fibersArrayName = ""
    sheetArrayName = ""
    sheetnormalArrayName = ""
    hasFibers = hasSheets = hasSheetNormals = False
    for i in range(data.GetCellData().GetNumberOfArrays()):
        name = data.GetCellData().GetArray(i).GetName()
        if name in ("Fiber", "DifferenceVector"):
            fibersArrayName, hasFibers = name, True
        if name == "Sheet":
            sheetArrayName, hasSheets = name, True
        if name == "Sheetnormal":
            sheetnormalArrayName, hasSheetNormals = name, True

    basesLines = []
    NumQP = 1
    if hasFibers:
        print("Bases:", data.GetNumberOfCells())
        print("  Fiber information found in array", "\"" + fibersArrayName + "\"")
        cid = 0
        for vtkcid in range(data.GetNumberOfCells()):
            cell = data.GetCell(vtkcid)
            npts = cell.GetNumberOfPoints()
            if npts == 10:
                hasT10 = True
            if npts in (4, 10):
                idx = cid + 1 if START_FROM_ONE else cid
                f = data.GetCellData().GetArray(fibersArrayName)
                if hasSheets and hasSheetNormals:
                    s = data.GetCellData().GetArray(sheetArrayName)
                    sn = data.GetCellData().GetArray(sheetnormalArrayName)
                    m = NormFiberSheetNormal(f.GetTuple3(vtkcid), s.GetTuple3(vtkcid), sn.GetTuple3(vtkcid))
                else:
                    m = createONS(f.GetTuple3(vtkcid))
                NumQP = 5 if hasT10 else 1
                parts = [str(idx)]
                parts.extend(str(m[k]) for _ in range(NumQP) for k in range(9))
                basesLines.append(" ".join(parts))
                cid += 1
        basesLines.append(f"{cid} {NumQP}")

    print("Surfaces:", data.GetNumberOfCells())
    # same code as for elements, except number of nodes per cell
    surLines = []
    scid = 0
    for vtkscid in range(data.GetNumberOfCells()):
        cell = data.GetCell(vtkscid)
        pids = cell.GetPointIds()
        npts = cell.GetNumberOfPoints()
        if npts in (3, 6):
            idx = scid + 1 if START_FROM_ONE else scid
            parts = [str(idx)]
            if START_FROM_ONE:
                for pid in range(pids.GetNumberOfIds()):
                    parts.append(str(pids.GetId(pid) + 1))
                mat = int(data.GetCellData().GetArray("Material").GetTuple1(vtkscid)) if hasMaterial else 0
                parts.append(f"{mat} {mat}")
            else:
                for pid in range(pids.GetNumberOfIds()):
                    parts.append(str(pids.GetId(pid)))
                parts.append("10 1")
            surLines.append(" ".join(parts))
            scid += 1
    surLines.append(f"{scid} 3 2")

    outFileName = args.outfile if args.outfile else os.path.splitext(os.path.basename(filename))[0]
    _write(f"{outFileName}.node", nodeLines)
    _write(f"{outFileName}.ele", eleLines)
    _write(f"{outFileName}.sur", surLines)
    if hasFibers:
        _write(f"{outFileName}.bases", basesLines)


if __name__ == "__main__":
    main()
