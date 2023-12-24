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

import vtk
import sys, os
import argparse
import numpy as np

parser = argparse.ArgumentParser(description = 'Convert a vtk file to .node and .ele files for CardioMechanics.')
parser.add_argument('mesh', help = 'VTK input file.')
parser.add_argument('-outfile', type = str, help = 'Set an alternative output filename PREFIX if desired.')
parser.add_argument('-scale', type = int, help = 'Scale the VTK input file by a factor SCALE.')
parser.add_argument('-fixMaterial', nargs = '+', type = int, help = 'Fix nodes if they are connected to cells that match the given mask in "Material". Example: "-fixMaterial 32 162"')
args = parser.parse_args()

# CardioMechanics ids start from 1
# vtk ids start from 0

start_from_one = True

filename=args.mesh
unitScaling = 1
if args.scale:
        unitScaling=float(args.scale)

def nodeIsFixed(data,array,Point, fixMaterial):
    pointCellList = vtk.vtkIdList()
    data.GetPointCells(Point,pointCellList)
    for iPointCell in xrange(pointCellList.GetNumberOfIds()):
        if array.GetValue(pointCellList.GetId(iPointCell)) in fixMaterial: 
            return 7
    return 0

# creates an orientation matrix from a given fiber vector
def createONS(f):
    f = f/np.linalg.norm(f)
    [a, b, c] = [ abs(np.dot(f,[1,0,0])), abs(np.dot(f,[0,1,0])), abs(np.dot(f,[0,0,1])) ]
    s = np.array([0,0,0])
    sn = np.array([0,0,0])
    if a<=b and a<=c:
            s = np.array([1,0,0])
    elif b<=a and b<=c:
            s = np.array([0,1,0])
    else:
            s = np.array([0,0,1])
    sn = np.cross(f,s)
    sn = sn/np.linalg.norm(sn)
    s = np.cross(sn,f)
    s = s/np.linalg.norm(s)
    return [ f[0], f[1], f[2], s[0], s[1], s[2], sn[0], sn[1], sn[2] ]

def NormFiberSheetNormal(f,s,sn):
    f = f/np.linalg.norm(f)
    s = s/np.linalg.norm(s)
    sn = sn/np.linalg.norm(sn)
    return [ f[0], f[1], f[2], s[0], s[1], s[2], sn[0], sn[1], sn[2] ]


if (filename.lower().endswith('vtu')):  # CardioMechanics output files
    reader = vtk.vtkXMLUnstructuredGridReader()
else:
    reader = vtk.vtkUnstructuredGridReader()
reader.SetFileName(filename)
reader.Update()

data = reader.GetOutput()

if args.fixMaterial:
    matArray = data.GetCellData().GetArray('Material')

hasFixation = False
if data.GetPointData().GetArray("Fixation") != None:
    print("Node fixation found")
    hasFixation = True

hasMaterial = False
if data.GetCellData().GetArray("Material") != None:
    print("Material information found")
    hasMaterial = True

print("Nodes:", data.GetNumberOfPoints())
nodeLines = []
for pid in range( 0,data.GetNumberOfPoints() ):
    coords = unitScaling * np.array(data.GetPoint(pid))
    fixation = 0
    if args.fixMaterial:
        fixation = nodeIsFixed(data, matArray, pid, args.fixMaterial)
    if hasFixation:
        fixation = int(data.GetPointData().GetArray("Fixation").GetTuple1(pid))
    if start_from_one:
        nodeLines.append( str(pid+1)+" "+str(coords[0])+" "+str(coords[1])+" "+str(coords[2])+" "+str(fixation) )
    else:
        nodeLines.append( str(pid)+" "+str(coords[0])+" "+str(coords[1])+" "+str(coords[2])+" "+str(fixation) )
nodeLines.append( str(data.GetNumberOfPoints())+" 3 1 0" )

hasT10 = False

print("Elements:", data.GetNumberOfCells())
eleLines=[]
cid = 0
for vtkcid in range( 0, data.GetNumberOfCells() ):
    eleLine=""
    cell = data.GetCell(vtkcid)
    pids = cell.GetPointIds()
    if (cell.GetNumberOfPoints()==10 and hasT10==False):
        hasT10 = True
        print("  T10 elements found")
    if (cell.GetNumberOfPoints() == 4 or cell.GetNumberOfPoints()==10):
        if start_from_one:
            eleLine = eleLine + str(cid+1)+" "
            for pid in range( 0, pids.GetNumberOfIds() ):
                eleLine = eleLine + str(pids.GetId(pid)+1)+" "
        else:
            eleLine = eleLine + str(cid)+" "
            for pid in range( 0, pids.GetNumberOfIds() ):
                eleLine = eleLine + str(pids.GetId(pid))+" "
        mat = 0
        if hasMaterial:
            mat = int(data.GetCellData().GetArray("Material").GetTuple1(vtkcid))
        eleLine = eleLine + str(mat)
        cid = cid + 1
        eleLines.append(eleLine)
if not hasT10:
    eleLines.append(str(cid)+" 4 1")
else:
    eleLines.append(str(cid)+" 10 1")

fibersArrayName = ""
hasFibers = False
hasSheets = False
hasSheetNormals = False
basesLines=[]
for i in range( 0,data.GetCellData().GetNumberOfArrays() ):
    a = data.GetCellData().GetArray(i)
    if a.GetName() == "Fiber" or a.GetName() == "DifferenceVector":
        fibersArrayName = a.GetName()
        hasFibers = True
    if a.GetName() == "Sheet":
        sheetArrayName = a.GetName()
        hasSheets = True
    if a.GetName() == "Sheetnormal":
        sheetnormalArrayName = a.GetName()
        hasSheetNormals = True
if hasFibers:
    print( "Bases:", data.GetNumberOfCells())
    print( "  Fiber information found in array", "\""+fibersArrayName+"\"")
    cid = 0
    for vtkcid in range( 0, data.GetNumberOfCells() ):
        basesLine=""
        cell = data.GetCell(vtkcid)
        if (cell.GetNumberOfPoints()==10):
            hasT10 = True
        if (cell.GetNumberOfPoints() == 4 or cell.GetNumberOfPoints()==10):
            if start_from_one:
                basesLine = basesLine + str(cid+1)+" "
            else:
                basesLine = basesLine + str(cid)+" "
            f = data.GetCellData().GetArray(fibersArrayName)
            if hasSheets and hasSheetNormals:
                s = data.GetCellData().GetArray(sheetArrayName)
                sn = data.GetCellData().GetArray(sheetnormalArrayName)
                m = NormFiberSheetNormal(f.GetTuple3(vtkcid),s.GetTuple3(vtkcid),sn.GetTuple3(vtkcid))
            else:
                m = createONS(f.GetTuple3(vtkcid))
            if hasT10:
                NumQP = 5
            else:
                NumQP = 1
            for j in range(NumQP):
                for i in range(0,9):
                    basesLine = basesLine + " " + str(m[i])
            cid = cid + 1
            basesLines.append(basesLine)
    basesLines.append(str(cid)+" "+ str(NumQP))

print("Surfaces:", data.GetNumberOfCells())
# same code as for elements, except number of nodes per cell
surLines=[]
scid = 0
for vtkscid in range( 0, data.GetNumberOfCells() ):
    surLine = ""
    cell = data.GetCell(vtkscid)
    pids = cell.GetPointIds()
    if (cell.GetNumberOfPoints() == 3 or cell.GetNumberOfPoints() == 6):
        if start_from_one:
            surLine = surLine + str(scid+1)+" "
            for pid in range( 0, pids.GetNumberOfIds() ):
                surLine = surLine + str(pids.GetId(pid)+1)+" "
            mat = 0
            if hasMaterial:
                    mat = int(data.GetCellData().GetArray("Material").GetTuple1(vtkscid))
            surLine = surLine + str(mat) + " " + str(mat)
        else:
            surLine = surLine + str(scid)+" "
            for pid in range( 0, pids.GetNumberOfIds() ):
                surLine = surLine + str(pids.GetId(pid))+" "
            surLine = surLine + "10 1"
        scid = scid + 1
        surLines.append(surLine)
surLines.append(str(scid)+" 3 2")

if args.outfile:
    OutFileName = args.outfile
else:
    OutFileName = os.path.splitext(os.path.basename(filename))[0]

f=open("{}.node".format(OutFileName),"w")
f.write( nodeLines.pop()+"\n" )
for l in nodeLines:
        f.write(l+"\n")
f.close()

f=open("{}.ele".format(OutFileName),"w")
f.write( eleLines.pop()+"\n" )
for l in eleLines:
        f.write(l+"\n")
f.close()

f=open("{}.sur".format(OutFileName),"w")
f.write( surLines.pop()+"\n" )
for l in surLines:
        f.write(l+"\n")
f.close()

if hasFibers:
        f=open("{}.bases".format(OutFileName),"w")
        f.write( basesLines.pop()+"\n" )
        for l in basesLines:
                f.write(l+"\n")
        f.close()
