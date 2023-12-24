#!/usr/bin/env python

'''
Author: Tobias Gerach
Date: 17.11.2020

TuneCV utility for tetrahedral meshes in use with acCELLerate.
Based on the method proposed in

Costa CM, Hoetzl E, Rocha BM, Prassl AJ, Plank G.
Automatic Parameterization Strategy for Cardiac Electrophysiology Simulations.
Comput Cardiol (2010). 2013;40:373-376.

Dependencies:
Gmsh (currently with hardcoded location /Volumes/bordeaux/ServerApps/Gmsh.app/Contents/MacOS/gmsh)
VTK
acCELLerate (Lattice CellModel XCELLentTools LatticePETSc)
'''

import argparse
import os
import vtk
from vtk.util import numpy_support as VN
import subprocess
import logging


def parse():
    parser = argparse.ArgumentParser(
        description='TuneCV utility for tetrahedral meshes in use with acCELLerate.')

    # mesh settings
    parser.add_argument('-resolution',
                        type=float, default=0.5,
                        help='mesh resolution in mm (default is 0.5 mm)')
    parser.add_argument('-length',
                        type=float, default=20.0,
                        help='length of strand in mm (default is 20.0 mm)')

    # ionic model settings
    parser.add_argument('model',
                        type=str,
                        choices=['BeelerReuter', 'TenTusscher2', 'CourtemancheEtAl', 'OHaraRudy',
                                 'HimenoEtAl', 'KoivumaekiEtAl', 'GrandiEtAlVentricle', 'GrandiEtAlAtrium'],
                        help='Ionic model')
    parser.add_argument('amplitude',
                        type=float,
                        help='Stimulus amplitude')

    parser.add_argument('duration',
                        type=float,
                        help='Stimulus duration in s')

    # conductivity/velocity settings
    parser.add_argument('-velocity',
                        type=float, default=0.67,
                        help='desired conduction velocity in m/s (default: 0.67 m/s)')
    parser.add_argument('-gi',
                        type=float, default=0.174,
                        help='initial value intracellular conductivity in S/m (default: 0.174 S/m)')
    parser.add_argument('-ge',
                        type=float, default=0.625,
                        help='initial value extracellular conductivity in S/m (default: 0.625 S/m)')
    parser.add_argument('-gm',
                        type=float, default=0.136,
                        help='initial value monodomain conductivity (harmonic mean of gi/ge) in S/m (default: 0.136 S/m)')
    parser.add_argument('-beta',
                        type=float, default=140000,
                        help='Membrane surface-to-volume ratio (default: 140000 m^-1)')
    parser.add_argument('-cm',
                        type=float, default=0.01,
                        help='Membrane capacitance per unit area (default: 0.01 Fm^-2)')
    parser.add_argument('-sourceModel',
                        default='mono',
                        choices=['mono'],
                        help='pick type of electrical source model (default is monodomain)')
    parser.add_argument('-tol',
                        type=float, default=0.01,
                        help='Percentage of error in conduction velocity (default: 0.01)')
    parser.add_argument('-converge',
                        default=False, action='store_true',
                        help='Iterate until converged velocity (default: False)')
    parser.add_argument('-maxit',
                        type=int, default=20,
                        help='Maximum number of conduction velocity iterations (default: 20)')

    # numerical settings
    parser.add_argument('-dt',
                        type=float, default=1.e-5,
                        help='Integration time step on seconds (default: 1e-5 s)')

    parser.add_argument('-theta',
                        type=float, default=0.5,
                        choices=[0, 0.5, 1],
                        help='Choose time stepping method 0: explicit, 0.5: CN, 1: implicit (default: 0.5)')

    parser.add_argument('-stol',
                        type=float, default=1.e-6,
                        help='Solver tolerance (default: 10^-6)')

    parser.add_argument('-lumping',
                        type=bool, default=False,
                        choices=[True, False],
                        help='Use mass lumping (default: False)')

    parser.add_argument('-activationThreshold',
                        type=float, default=-0.02,
                        help='Threshold to record LAT (default: -0.02 V)')

    # output settings
    parser.add_argument('-log',
                        default='tuneCV.log',
                        help='Generate/append to a log file.')

    parser.add_argument('-gmshexec',
                        type=str, default='/Applications/ServerApps/Gmsh.app/Contents/MacOS/gmsh',
                        help='Use this option to define an alternative Gmsh path.')

    return parser


def jobID(args):
    """
    Generate name of top level output directory.
    """
    tpl = '{}_vel_{}_dx_{}'
    return tpl.format(args.model, args.velocity, args.resolution)


def readVTK(filename):
    logging.info('Reading file:   {}'.format(filename))
    reader = vtk.vtkUnstructuredGridReader()
    if (filename.lower().endswith('vtu')):  # CardioMechanics output files
        reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(filename)
    reader.Update()
    return reader.GetOutput()


def writeVTK(data, filename, datatype='vtu'):
    logging.info('Writing file:   {}'.format(filename))
    if datatype == 'vtu':
        writer = vtk.vtkXMLUnstructuredGridWriter()
    else:
        writer = vtk.vtkUnstructuredGridWriter()
    writer.SetFileName(filename)
    writer.SetInputData(data)
    writer.Write()


def addCellArray(ArrayName, value, mesh):
    noCells = mesh.GetNumberOfCells()
    Array = vtk.vtkDoubleArray()
    Array.SetName(ArrayName)
    Array.SetNumberOfComponents(1)
    for iCell in range(noCells):
        Array.InsertNextValue(value)
    mesh.GetCellData().AddArray(Array)

    return mesh


def Init():
    """
    Initialize arguments and directory
    """
    parser = parse()
    args = parser.parse_args()

    if os.path.isdir(jobID(args)):
        pass
    else:
        try:
            os.mkdir(jobID(args))
        except:
            raise OSError("Could not read or create top level directory!")

    print("-- -- -- -- -- -- -- -- Initialize tuneCV -- -- -- -- -- -- -- --")
    print("Parameterlist:")
    for arg in vars(args):
        print("{:12s}\t{}".format(arg, getattr(args, arg)))
    print("-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --")

    return args


def CreateMesh(args):
    """
    Create the mesh using Gmsh
    """
    # create .geo file
    f = open("{}/mesh.geo".format(jobID(args)), "w+")
    f.write('SetFactory("OpenCASCADE");\n')
    f.write('lc = {};\n'.format(args.resolution))
    f.write('x = {};\n'.format(args.length))
    f.write('y = lc;\nz = lc;\n')
    f.write(
        'Point(1) = {0.0, 0.0, 0.0, lc};\nPoint(2) = {0.0, y, 0.0, lc};\n//+\nLine(1) = {1, 2};\n')
    f.write('//+\nExtrude {0, 0, z} { Curve{1}; Layers{z/lc}; }\n')
    f.write('//+\nExtrude {x, 0, 0} { Surface{1}; Layers{x/lc}; }\n')
    f.write('//+\nPhysical Volume("tets") = {1};\n')
    f.close()

    # run Gmsh to create .vtk
    subprocess.call([args.gmshexec, "{}/mesh.geo".format(jobID(args)), "-3",
                     "-format", "vtk", "-v", "0", "-o", "{}/mesh.vtk".format(jobID(args))])

    # create material and fiber arrays
    mesh = readVTK("{}/mesh.vtk".format(jobID(args)))
    mesh = addCellArray('Material', 31, mesh)
    mesh = addCellArray('Phi', 0, mesh)
    mesh = addCellArray('Theta', 127, mesh)

    writeVTK(mesh, "{}/mesh.vtu".format(jobID(args)), 'vtu')


def CreateCellModelFile(args):
    """
    Create .cmf file with chosen model
    """
    f = open("{}/CellModelFile.cmf".format(jobID(args)), "w+")
    f.write('31 {} 0\n'.format(args.model))
    f.close()


def CreateLogFile(args):
    """
    Create .log file for iterative run
    """
    f = open("{}/{}".format(jobID(args), args.log), "w+")
    f.write('Iteration\tConductivity (S/m)\tCV (m/s)\tResidual')
    f.close()


def CreateConditionFile(args):
    """
    Create .cnd and .acnd files
    """
    f = open("{}/ConditionFile.cnd".format(jobID(args)), "w+")
    f.write('-1.0-0.0 0.0-{} 0.0-{} {} 0 {} 0 Ii\n'.format(args.resolution,
                                                           args.resolution, args.amplitude, args.duration))
    f.close()

    subprocess.call(['GeometricCnd2Acnd', "{}/mesh.vtu".format(jobID(args)), "{}/ConditionFile.cnd".format(jobID(args)),
                     "{}/AcCELLerateConditionFile.acnd".format(jobID(args))], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def CreateMaterialFile(sigma, iteration, args):
    """
    Create .def file with the conductivity sigma
    """
    f = open("{}/Iter_{}/materialIntra.def".format(jobID(args), iteration), "w+")
    f.write('[Slab 31 31]\n\n')
    f.write('Kappa\t0\t{}\n'.format(sigma))
    f.write('AnisotropyX\t1.0\nAnisotropyY\t1.0\nAnisotropyZ\t1.0\n')
    f.close()

    # Set up system and mass matrix
    subprocess.call(['BidomainMatrixGenerator', "{}/mesh.vtu".format(jobID(args)), "-o", "{}/Iter_{}/mesh".format(jobID(args), iteration),
                     '-mono', "{}/Iter_{}/materialIntra.def".format(jobID(args), iteration)], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def CreateACLTFile(args, iteration):
    """
    Create .aclt file
    """
    f = open("{}/Iter_{}/settings.aclt".format(jobID(args), iteration), "w+")
    f.write('Condition {}/AcCELLerateConditionFile.acnd\n'.format(jobID(args)))
    # Set to at least 2x the theoretical propagation time
    CalcLength = round(
        args.dt * round((2 * (args.length * 1e-3 / args.velocity)) / args.dt), 7)
    f.write('CalcLength {}\n'.format(CalcLength))
    f.write('DTCell {}\n'.format(args.dt))
    f.write('DTIntra {}\n'.format(args.dt))
    f.write('CellModelFile {}/CellModelFile.cmf\n'.format(jobID(args)))
    f.write('Domain {}\n'.format(args.sourceModel))
    f.write('ThetaIntra {}\n'.format(args.theta))
    f.write('MembraneCapacitance {}\n'.format(args.cm))
    f.write('SurfaceToVolume {}\n'.format(args.beta))
    f.write('ActivationThreshold {}\n'.format(args.activationThreshold))
    f.write('Results AT\n')
    f.write('Material {}/Iter_{}/mesh.vec\n'.format(jobID(args), iteration))
    f.write('MatrixIntra {}/Iter_{}/mesh.mat\n'.format(jobID(args), iteration))
    f.write('MassMatrixIntra {}/Iter_{}/mesh.mass.mat\n'.format(jobID(args), iteration))
    f.write('Resprefix {}/Iter_{}/r\n'.format(jobID(args), iteration))
    f.write('Protocol {}/Iter_{}/Sim.log\n'.format(jobID(args), iteration))
    f.write('NoExport')
    f.close()


def StartSimulation(sigma, iteration, args):
    """
    Wrapper function to start a single simulation
    """
    # Create result directory
    if os.path.isdir("{}/Iter_{}".format(jobID(args), iteration)):
        pass
    else:
        try:
            os.mkdir("{}/Iter_{}".format(jobID(args), iteration))
        except:
            raise OSError("Could not create result directory!")

    CreateMaterialFile(sigma, iteration, args)
    CreateACLTFile(args, iteration)

    # Run acCELLerate
    subprocess.call(['acCELLerate', "{}/Iter_{}/settings.aclt".format(jobID(args),
                                                                      iteration)], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def EvaluateCV(args, iteration):
    """
    Calculate CV estimate from simulation
    """
    # Convert PETScVec to VTK
    subprocess.call(['PETScVec2VTK', "{}/Iter_{}/r_AT.vec".format(jobID(args), iteration), "{}/mesh.vtu".format(jobID(args)),
                     "{}/Iter_{}/AT.vtu".format(jobID(args), iteration)], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    # Process activation times
    ActivationTimes = readVTK(
        "{}/Iter_{}/AT.vtu".format(jobID(args), iteration))
    LAT = VN.vtk_to_numpy(
        ActivationTimes.GetPointData().GetAbstractArray("Results"))
    firstActivation = min(LAT)
    lastActivation = max(LAT)
    if (firstActivation < 0 and lastActivation > 0):
        return -1
    elif (firstActivation < 0 and lastActivation < 0):
        return -2
    propagationTime = lastActivation - firstActivation
    propagationDistance = args.length * 1.e-3

    return propagationDistance / propagationTime


def UpdateConductivity(sigma, cv, args):
    """
    Calculate new sigma_m
    """
    sigma *= (args.velocity / cv)**2

    return sigma


def run():
    # Initialize tuneCV
    args = Init()
    # Init logfile
    CreateLogFile(args)
    # Prepare Mesh
    CreateMesh(args)
    # Create static acCELLerate settings
    CreateCellModelFile(args)
    CreateConditionFile(args)
    # Initial CV estimation
    iteration = 0
    sigma = [args.gm]
    CV = []
    StartSimulation(sigma[-1], iteration, args)
    print("Initial estimate:")
    CV_estimate = EvaluateCV(args, iteration)
    if CV_estimate == -1:
        raise RuntimeError(
            "Simulation time is too short. Increase initial sigma to get a better estimate.")
    elif CV_estimate == -2:
        raise RuntimeError(
            "Tissue was not activated. For high sigma values, try to increase stimulus amplitude.")
    else:
        CV.append(CV_estimate)
    residual = abs(CV[-1] - args.velocity) / args.velocity
    print(
        "sigma_m = {:.4f} S/m\tvelocity = {:.4f} m/s".format(sigma[-1], CV[-1]))
    # Tune CV
    if args.converge:
        print("Start Iterating...")
        while (iteration < args.maxit and residual > args.tol):
            # update .log file
            f = open("{}/{}".format(jobID(args), args.log), "a")
            f.write('\n{}\t{:.4f}\t{:.4f}\t{:.4f}'.format(
                iteration, sigma[-1], CV[-1], residual))
            f.close()
            iteration += 1
            sigma.append(UpdateConductivity(sigma[-1], CV[-1], args))
            StartSimulation(sigma[-1], iteration, args)
            CV_estimate = EvaluateCV(args, iteration)
            if CV_estimate == -1:
                raise RuntimeError(
                    "Simulation time is too short. Increase initial sigma to get a better estimate.")
            elif CV_estimate == -2:
                raise RuntimeError(
                    "Tissue was not activated. For high sigma values, try to increase stimulus amplitude.")
            else:
                CV.append(CV_estimate)
            print(
                "sigma_m = {:.4f} S/m\tvelocity = {:.4f} m/s".format(sigma[-1], CV[-1]))
            residual = abs(CV[-1] - args.velocity) / args.velocity
        # finish optimization with update of log file
        print("DONE. Residual lower than tolerance.")
        f = open("{}/{}".format(jobID(args), args.log), "a")
        f.write('\n{}\t{:.4f}\t{:.4f}\t{:.4f}'.format(
            iteration, sigma[-1], CV[-1], residual))
        f.close()


if __name__ == '__main__':
    run()
