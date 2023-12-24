#!/opt/local/bin/python

'''
Author:         Tobias Gerach
Date:           01.10.18
Description:    This script automatically converts all *.VEC files in a given directory to *.VTU format
                using PETScVec2VTK. Afterwards it creates an XML file that can be used to visualize all
                converted files in, e.g. ParaView.
-------------------------------------------------------------------------------------------------------
Changelog:
Version 1.0 tg319
'''

from optparse import OptionParser
import subprocess
import glob
import os, time, sys
import xml.etree.cElementTree as ET
import re

def convert(options):
    print("Currently trying to convert files...")

    # get all *.vec files in the target directory
    print("Search for *.vec files in directory './{}'".format(options.directory))
    files = sorted(glob.glob('{}/*_Vm.vec'.format(options.directory)))
    print("{} files found".format(len(files)))

    # setting regular expressions
    timestamp = re.compile(r'_\d\d\d\d.\d\d\d\d\d\d_')
    directory = re.compile(r'{}/'.format(options.directory))

    for file in files:
        filename = re.split(directory, file)[-1]
        res_prefix = re.split(timestamp, filename)[0]
        res_suffix = os.path.splitext(re.split(timestamp, filename)[-1])[0]

    # loop over all *.vec files
    print("Begin to convert all files to VTU...")
    for i,file in enumerate(files, start = 1):
        # target_name = os.path.splitext(file)[0] + '.vtu'
        target_name = "{}/Vm.{}.vtu".format(options.directory,i)
        print(os.path.splitext(file)[0])
        if options.clean_flag:
            subprocess.call([ 'PETScVec2VTK', file, '{}'.format(options.geometry), target_name, "--arrayName", "Vm", "--delete" ])
        else:
            subprocess.call([ 'PETScVec2VTK', file, '{}'.format(options.geometry), target_name, "--arrayName", "Vm" ])

    print("Conversion finished successfully.\n")

def save_as_pvd(options):
    print("Create an XML file containing the simulation results\n")
    print("Search for *.vtu files in directory './{}'".format(options.directory))
    # filenames = sorted(glob.glob('./{}/*.vtu'.format(options.directory)))
    filenames = sorted(glob.glob('{}/*_Vm.vec'.format(options.directory)))
    print("{} files found".format(len(filenames)))

    # setting regular expressions
    timestamp = re.compile(r'\d\d\d\d.\d\d\d\d\d\d')

    VTKFile = ET.Element("VTKFile", type = "Collection", version = "0.1", byte_order = "LittleEndian", compressor = "vtkZLibDataCompressor")
    Collection = ET.SubElement(VTKFile, "Collection")

    for i,filename in enumerate(filenames, start = 1):
        file = os.path.split(filename)[1]
        timestep =  float(timestamp.findall(file)[0])
        file = "Vm.{}.vtu".format(i)
        ET.SubElement(Collection, "DataSet", timestep = "{}".format(timestep), part = "0", group = "", file = "{}".format(file))

    tree = ET.ElementTree(VTKFile)
    tree.write('./{0}/{1}.pvd'.format(options.directory, options.outfile))
    if os.path.isfile('./{0}/{1}.pvd'.format(options.directory, options.outfile)):
        print("Created XML file successfully")
    else:
        raise OSError('File could not be saved')


def main():
    parser = OptionParser(  usage = "Usage: %prog [options]",
                            description = "%prog automatically converts all *.vec files in the provided directory to *.vtu format and creates a *.pvd file containing all results ready to be viewed in ParaView.",
                            version = "%prog 1.0"   )
    parser.add_option(  "-d", "--directory",
                        action = "store", dest = "directory", metavar = "PATH", type = "string",
                        default = None,
                        help = "Directory containing simulation results"    )
    parser.add_option(  "-g", "--geometry",
                        action = "store", dest = "geometry", metavar = "GEOFILE", type = "string",
                        default = None,
                        help = "Geometry file (*.vtu) used for the simulation"  )
    parser.add_option(  "-o", "--outfile",
                        action = "store", dest = "outfile", metavar = "FILENAME", type = "string",
                        default = "results",
                        help = "Filename of the resulting *.pvd file"   )
    parser.add_option(  "-c", "--clean",
                        action = "store_true", dest = "clean_flag", default = False,
                        help = "Choose wether or not to delete converted *.vec files"   )

    (options, args) = parser.parse_args()

    # Check for major parsing errors
    if options.geometry == None:
        parser.error("Geometry file was not specified.")
    if options.directory == None:
        parser.error("Directory was not specified.")
    if os.path.isdir(options.directory):
        pass
    else: 
        parser.error("Directory does not exist.")
    if os.path.isfile(options.geometry):
        pass
    else: 
        parser.error("Geometry file does not exist.")

    print("Parsing finished")

    # Run conversion and save PVD file
    convert(options)
    save_as_pvd(options)


if __name__ == '__main__':
    main()