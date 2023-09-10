//

//  VTK2PETScVec.cpp
//  LatticePETSc
//
/*! \file VTK2PETScVec.cpp
   \brief Generate PETSc vector from VTK surface or volumetric mesh

   \author Axel Loewe, IBT, KIT
 */


#include <PETScLSE.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkPointSet.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkDoubleArray.h>
#include <unordered_set>
#include <petscviewer.h>
#include <iostream>
#include <cstdio>

template<class DataType>

void PETScVec2VTK(DataType *pointSet, std::string vecName, std::string arrayName) {
#if KADEBUG
  fprintf(stderr, "PETScVec2VTK<ValTyp>::PETScVec2VTK()\n");
#endif  // if KADEBUG

  PetscErrorCode ierr;

  Vec v;
  PetscViewer fd;
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD, vecName.c_str(), FILE_MODE_READ, &fd); CHKERRQ(ierr);
  VecCreate(PETSC_COMM_WORLD, &v);
  ierr = VecLoad(v, fd); CHKERRQ(ierr);
  PetscInt nPoints, x, y, nPointsVec;
  std::string origVtkFilename;
  kaSharedMatrixN<double> m;  // Geometry matrix
  // LoadLatInfo(vecName.c_str(), m, x, y, nPoints);
  nPoints = pointSet->GetNumberOfPoints();
  VecGetSize(v, &nPointsVec);
  ierr = PetscViewerDestroy(&fd); CHKERRQ(ierr);

  // if(pointSet->GetNumberOfPoints() != nPoints)
  //    throw kaBaseException("Number of points %u does not match VTK file %u.\nOriginal VTK file: %s.", nPoints,
  // pointSet->GetNumberOfPoints(), origVtkFilename.c_str());
  /*if(pointSet->GetNumberOfCells() != nCells)
      throw kaBaseException("Number of cells %u does not match VTK file %u.\nOriginal VTK file: %s.", nCells,
         pointSet->GetNumberOfCells(), origVtkFilename.c_str());*/
  if (nPointsVec != nPoints) {
    throw kaBaseException("Number of vector elements %u does not match number of points %u.\nOriginal VTK file: %s.",
                          nPointsVec, nPoints, origVtkFilename.c_str());
  }

  // PETScVec is always double. Maybe give option to user to store result in using smaller data types
  vtkSmartPointer<vtkDoubleArray> results = vtkSmartPointer<vtkDoubleArray>::New();
  results->SetName(arrayName.c_str());
  results->SetNumberOfComponents(1);
  results->SetNumberOfValues(nPoints);
  PetscScalar *pv;
  ierr = VecGetArray(v, &pv); CHKERRQ(ierr);
  for (unsigned int i = 0; i < nPoints; i++)
    results->SetValue(i, *pv++);

  pointSet->GetPointData()->AddArray(results);
  ierr = VecRestoreArray(v, &pv); CHKERRQ(ierr);
  ierr = VecDestroy(&v); CHKERRQ(ierr);
}  // PETScVec2VTK

int main(int argc, char **argv) {
#if KADEBUG
  fprintf(stderr, "main\n");
#endif  // if KADEBUG
  PetscErrorCode ierr;
  PetscInitialize(&argc, &argv, (char *)0, "");
  bool doDelete = false;

  if (argc < 4) {
    cerr << argv[0] << " <vector filename> <VTK filename (.vtp/.vtu) in>  <VTK filename (.vtp/.vtu) out>" << endl;
    cerr << "\t[--arrayName]\tName of the VTK array to be appended (default: Results)" << endl;
    cerr << "\t[--delete]\tDelte PETScVec after conversion" << endl;
    exit(-1);
  }
  try {
    std::string arrayName = "Results";
    for (int i = 4; i < argc; i++)
      if (strcmp(argv[i], "--arrayName") == 0) {
        arrayName = argv[++i];
      } else {
        if (strcmp(argv[i], "--delete") == 0)
          doDelete = true;
        else
          throw kaBaseException("Unknown option %s", argv[i]);
      }

    std::string filenameVtkIn  = argv[2];
    std::string filenameVtkOut = argv[3];
    std::string extension      = filenameVtkIn.substr(filenameVtkIn.length()-3);
    std::string extensionOut   = filenameVtkOut.substr(filenameVtkOut.length()-3);
    if (extension.compare(extensionOut) != 0) {
      throw kaBaseException("Extension of input VTK file (%s) does not match that of output file (%s).", argv[2],
                            argv[3]);
    }
    std::string filenameVec = argv[1];
    if (extension.compare("vtp") == 0) {
      vtkSmartPointer<vtkXMLPolyDataReader> reader =
        vtkSmartPointer<vtkXMLPolyDataReader>::New();
      reader->SetFileName(filenameVtkIn.c_str());
      reader->Update();
      vtkPolyData *pointSet = reader->GetOutput();
      PETScVec2VTK<vtkPolyData>(pointSet, filenameVec, arrayName);
      vtkSmartPointer<vtkXMLPolyDataWriter> writer =
        vtkSmartPointer<vtkXMLPolyDataWriter>::New();
      writer->SetFileName(filenameVtkOut.c_str());
#if VTK_MAJOR_VERSION > 5
      writer->SetInputData(pointSet);
#else  // if VTK_MAJOR_VERSION > 5
      writer->SetInput(pointSet);
#endif  // if VTK_MAJOR_VERSION > 5
      writer->Write();
    } else if (extension.compare("vtu") == 0) {
      vtkSmartPointer<vtkXMLUnstructuredGridReader> reader =
        vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
      reader->SetFileName(filenameVtkIn.c_str());
      reader->Update();
      vtkUnstructuredGrid *pointSet = reader->GetOutput();
      PETScVec2VTK<vtkUnstructuredGrid>(pointSet, filenameVec, arrayName);
      vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer =
        vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
      writer->SetFileName(filenameVtkOut.c_str());
#if VTK_MAJOR_VERSION > 5
      writer->SetInputData(pointSet);
#else  // if VTK_MAJOR_VERSION > 5
      writer->SetInput(pointSet);
#endif  // if VTK_MAJOR_VERSION > 5
      writer->Write();
    } else {
      throw kaBaseException("Unknown or unsupported type of input file %s", argv[1]);
    }
    if (doDelete) {
      remove(filenameVec.c_str());
      remove((filenameVec+".header").c_str());
      remove((filenameVec+".info").c_str());
    }
  } catch (kaBaseException &e) {
    cerr << argv[0] << " Error: " << e << endl;
    ierr = PetscFinalize(); CHKERRQ(ierr);
    exit(-1);
  }
  ierr = PetscFinalize(); CHKERRQ(ierr);
  return 0;
}  // main
