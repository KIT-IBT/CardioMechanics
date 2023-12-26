#!/usr/bin/env python

import argparse
import os


def parse():
    parser = argparse.ArgumentParser(
        prog='createSimDir',
        description='Little helper to create a basic simulation directory for electromechanical simulations.\n')

    parser.add_argument('-material',
                        type=int, default=[0], nargs='*',
                        help='List of material numbers (default: 0)')

    parser.add_argument('-model',
                        type=str, default=['ElphyDummy'], nargs='*',
                        choices=['BeelerReuter', 'TenTusscher2', 'CourtemancheEtAl', 'OHaraRudy',
                                 'HimenoEtAl', 'KoivumaekiEtAl', 'GrandiEtAlVentricle', 'GrandiEtAlAtrium'],
                        help='List of ionic models (default: ElphyDummy)')

    parser.add_argument('-sigma',
                        type=float, default=[0.136], nargs='*',
                        help='List of monodomain conductivities in S/m (default: 0.136 S/m)')

    parser.add_argument('-beta',
                        type=float, default=140000,
                        help='Membrane surface-to-volume ratio (default: 140000 m^-1)')

    parser.add_argument('-cm',
                        type=float, default=0.01,
                        help='Membrane capacitance per unit area (default: 0.01 Fm^-2)')

    parser.add_argument('-time',
                        type=float, default=0.8,
                        help='Simulation time in s (default: 0.8 s)')

    parser.add_argument('-dt',
                        type=float, default=1.e-5,
                        help='Integration time step in seconds (default: 1e-5 s)')

    parser.add_argument('-theta',
                        type=float, default=0.5,
                        choices=[0, 0.5, 1],
                        help='Choose time stepping method 0: explicit, 0.5: CN, 1: implicit (default: 0.5)')

    parser.add_argument('-sourceModel',
                        default='mono',
                        choices=['mono'],
                        help='pick type of electrical source model (default is monodomain)')

    parser.add_argument('-activationThreshold',
                        type=float, default=-0.02,
                        help='Threshold to record LAT (default: -0.02 V)')

    # output settings
    parser.add_argument('-dir',
                        type=str, default='SimDir',
                        help='Name of top level directory (default: SimDir)')

    return parser


def CreateCellModelFile(args):
    """
    Create .cmf file with chosen models
    """
    path = os.path.join(args.dir, 'cellmodelFiles')
    if os.path.isdir(path):
        pass
    else:
        os.mkdir(path)
    f = open("{}/CellModelFile.cmf".format(path), "w+")
    for idx, mat in enumerate(args.material):
        f.write('{} {} ForceDummy\n'.format(mat, args.model[idx]))
    f.close()


def CreateConditionFile(args):
    """
    Create .acnd files
    """
    path = os.path.join(args.dir, 'stimFiles')
    if os.path.isdir(path):
        pass
    else:
        os.mkdir(path)
    f = open("{}/stimulationFile.acnd".format(path), "w+")
    f.write('#pointID amplitude cycleLength duration activationTime Ii|Ie|If|Vm|Ve|Vf\n')
    f.close()


def CreateSensorFile(args):
    """
    Create a dummy sensor file 
    """
    path = os.path.join(args.dir, 'sensorFiles')
    if os.path.isdir(path):
        pass
    else:
        os.mkdir(path)
    f = open("{}/Sensor.txt".format(path), "w+")
    f.write('#pointID ExportPath tBeginSave dt VAR\n')
    f.close()


def CreateMaterialFile(args):
    """
    Create .def file with the conductivity sigma
    """
    path = os.path.join(args.dir, 'materialFiles')
    if os.path.isdir(path):
        pass
    else:
        os.mkdir(path)
    f = open("{}/materialIntra.def".format(path), "w+")
    for idx, mat in enumerate(args.material):
        f.write('[Tissue_{} {} {}]\n\n'.format(idx, mat, mat))
        f.write('Kappa\t0\t{}\n'.format(args.sigma[idx]))
        f.write('AnisotropyX\t1.0\nAnisotropyY\t1.0\nAnisotropyZ\t1.0\n\n')
    f.close()


def CreateACLTFile(args):
    """
    Create .aclt file
    """
    path = os.path.join(args.dir, 'settings')
    if os.path.isdir(path):
        pass
    else:
        os.mkdir(path)
    f = open("{}/settings.aclt".format(path), "w+")
    f.write('Condition ../stimFiles/stimulationFile.acnd\n')
    f.write('CalcLength {}\n'.format(args.time))
    f.write('DTCell {}\n'.format(args.dt))
    f.write('DTIntra {}\n'.format(args.dt))
    f.write('CellModelFile ../cellmodelFiles/CellModelFile.cmf\n')
    f.write('Domain {}\n'.format(args.sourceModel))
    f.write('ThetaIntra {}\n'.format(args.theta))
    f.write('MembraneCapacitance {}\n'.format(args.cm))
    f.write('SurfaceToVolume {}\n'.format(args.beta))
    f.write('ActivationThreshold {}\n'.format(args.activationThreshold))
    f.write('Material ../geoFiles/mesh.vec\n')
    f.write('MatrixIntra ../geoFiles/mesh.mat\n')
    f.write('MassMatrixIntra ../geoFiles/mesh.mass.mat\n')
    f.write('Results Vm AT\n')
    f.write('Resprefix ../Results/EP/r\n')
    f.write('NoExport')
    f.close()


def Init():
    """
    Initialize arguments and directory
    """
    parser = parse()
    args = parser.parse_args()

    if os.path.isdir(args.dir):
        pass
    else:
        try:
            os.mkdir(args.dir)
        except OSError as error:
            print(error)

    path = os.path.join(args.dir, 'geoFiles')
    if os.path.isdir(path):
        pass
    else:
        os.mkdir(path)

    path = os.path.join(args.dir, 'tetgen')
    if os.path.isdir(path):
        pass
    else:
        os.mkdir(path)

    print("-- -- -- -- -- -- Creating Simulation Directory -- -- -- -- -- --")
    print('The following structure of the directory is expected.')
    print('Files marked with * are not generated by this script.')
    print('Make sure to adapt all autogenerated files to your simulation.\n')
    print('{}'.format(args.dir))
    print('|- cellmodelFiles')
    print('|  |- CellModelFile.cmf')
    print('|- stimFiles')
    print('|  |- stimulationFile.acnd')
    print('|- materialFiles')
    print('|  |- materialIntra.def')
    print('|- sensorFiles')
    print('|  |- Sensor.txt')
    print('|- tetgen')
    print('|  |- *.node')
    print('|  |- *.ele')
    print('|  |- *.bases')
    print('|  |- *.sur')
    print('|- geoFiles')
    print('|  |- *.mass.mat')
    print('|  |- *.mass.mat.info')
    print('|  |- *.mat')
    print('|  |- *.mat.info')
    print('|  |- *.vec')
    print('|  |- *.vec.info')
    print('|  |- *.vtu')
    print('|- settings')
    print('   |- settings.aclt')
    print('   |- *.xml')
    print("-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --")

    return args


def run():
    args = Init()
    CreateCellModelFile(args)
    CreateConditionFile(args)
    CreateMaterialFile(args)
    CreateACLTFile(args)
    CreateSensorFile(args)


if __name__ == '__main__':
    run()
