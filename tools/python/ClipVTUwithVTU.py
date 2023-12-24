#!/opt/local/bin/python3

'''
Author: tg319
Date: 15.06.2021

Clip vtk/vtu/vtp time series using a predefined clip.
'''
import sys
import numpy as np
import vtk
import argparse
import xml.etree.cElementTree as ET


def read_VTX_XML(filename):
    if filename.endswith('vtk'):
        reader = vtk.vtkXMLDataSetReader()
    elif filename.endswith('vtu'):
        reader = vtk.vtkXMLUnstructuredGridReader()
    elif filename.endswith('vtp'):
        reader = vtk.vtkXMLPolyDataReader()
    else:
        IOError("Unknown FileEnding ")
    reader.SetFileName(filename)
    reader.Update()
    return reader.GetOutput()


def write_VTX_XML(filename, data):
    if filename.endswith('vtk'):
        writer = vtk.vtkXMLUnstructuredGridWriter()
    elif filename.endswith('vtu'):
        writer = vtk.vtkXMLUnstructuredGridWriter()
    elif filename.endswith('vtp'):
        writer = vtk.vtkXMLPolyDataWriter()
    writer.SetFileName(filename)
    writer.SetInputData(data)
    writer.Write()

def parse():
    parser = argparse.ArgumentParser(description='ClipVTUwithVTU clips the data and keeps the same elements.')
    parser.add_argument('mesh', help = 'vtk/vtu/vtp data you want to clip.')
    parser.add_argument('clip', help = 'clipped version of <mesh> created in e.g. ParaView invtu/vtk/vtp.')
    parser.add_argument('-outfile', type = str, default = 'ClipResult', help = 'prefix of the results file.')
    parser.add_argument('-start', type = int, default = 0, help = 'start value for the time series.')
    parser.add_argument('-stop', type = int, default = 100, help = 'stop value for the time series.')
    parser.add_argument('-step', type = int, default = 1, help = 'increment of the time series.')

    args = parser.parse_args()

    return args

def main():

    args = parse()

    heart_VTU_prefix = args.mesh
    clip_vtu_fn = args.clip
    out_prefix_fn = args.outfile
    start = args.start
    stop = args.stop
    i_step = args.step

    clip_VTU = read_VTX_XML(clip_vtu_fn)

    id_list = vtk.vtkIdTypeArray()
    for cid in range(clip_VTU.GetNumberOfCells()):
        cell_id = int(clip_VTU.GetCellData().GetArray("CellID").GetComponent(cid, 0))
        id_list.InsertNextValue(cell_id)

    for i in (range(start, stop, i_step)):
        print('{:03.2f} %'.format(100*(i-start)/(stop-start)), end="\r")
        heart_VTU = read_VTX_XML(heart_VTU_prefix + "." + str(i) + ".vtu")
        ori_id_name = "OriginalIDs"
        ori_ids = vtk.vtkIdFilter()
        ori_ids.SetInputData(heart_VTU)
        ori_ids.PointIdsOn()
        ori_ids.CellIdsOn()
        ori_ids.FieldDataOn()
        ori_ids.SetIdsArrayName(ori_id_name)

        SelectionNode = vtk.vtkSelectionNode()
        SelectionNode.SetFieldType(vtk.vtkSelectionNode.CELL)
        SelectionNode.SetContentType(vtk.vtkSelectionNode.INDICES)
        SelectionNode.SetSelectionList(id_list)
        selection = vtk.vtkSelection()
        selection.AddNode(SelectionNode)

        ExtractSelection = vtk.vtkExtractSelection()
        ExtractSelection.SetInputData(0, heart_VTU)
        ExtractSelection.SetInputData(1, selection)
        ExtractSelection.SetInputConnection(ori_ids.GetOutputPort())
        ExtractSelection.Update()
        selected = vtk.vtkUnstructuredGrid()
        selected.ShallowCopy(ExtractSelection.GetOutput())
        write_VTX_XML(out_prefix_fn + "." + str(i-start) + ".vtu", selected)

if __name__ == '__main__':
    main()
