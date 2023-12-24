/*
 * File: VTKIOHelper.h
 *
 * Institute of Biomedical Engineering, 
 * Karlsruhe Institute of Technology (KIT)
 * https://www.ibt.kit.edu
 * 
 * Repository: https://github.com/KIT-IBT/CardioMechanics
 *
 * License: GPL-3.0 (See accompanying file LICENSE or visit https://www.gnu.org/licenses/gpl-3.0.html)
 *
 */


#ifndef VTKIOHELPER_H
#define VTKIOHELPER_H

#include "Verbosity.h"
#include "ProgressBar.h"
#include "vtkLatticeSG.h"

#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkImageData.h>
#include <vtkPolyData.h>
#include <vtkDataSet.h>
#include <vtkXMLWriter.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkXMLDataReader.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkXMLImageDataReader.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkDataSetReader.h>

#include <string>

using namespace std;

class VTKIOHelper {
 public:
  template<class T> static void Write(vtkSmartPointer<T>, string);
  template<class T> static vtkSmartPointer<T> Read(string);
  static vtkSmartPointer<vtkDataSet> Read(std::string);

 private:
  static void write(vtkSmartPointer<vtkXMLWriter>, vtkSmartPointer<vtkDataSet>, string);

  template<class T> static vtkSmartPointer<T> read(vtkSmartPointer<vtkXMLDataReader>, string);
};

vtkSmartPointer<vtkDataSet> VTKIOHelper::Read(std::string fn) {
  std::string ext = fn.substr(fn.rfind('.')+1);

  vtkSmartPointer<vtkXMLDataReader> r;
  if (ext == "vtu") {
    r = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
  } else if (ext == "vtp") {
    r = vtkSmartPointer<vtkXMLPolyDataReader>::New();
  } else if (ext == "vti") {
    r = vtkSmartPointer<vtkXMLImageDataReader>::New();
  } else if (ext == "vtk") {
    vtkSmartPointer<vtkDataSetReader> rd = vtkSmartPointer<vtkDataSetReader>::New();
    rd->SetFileName(fn.c_str());
    rd->Update();
    return rd->GetOutput();
  } else {return NULL;}

  r->SetFileName(fn.c_str());
  r->Update();
  return r->GetOutputAsDataSet();
}

template<class T>

void VTKIOHelper::Write(vtkSmartPointer<T> mesh, string filename) {
  /* This is a default template should never get called (or even compiled)
   * If you encountered an error here, check your template class T
   * or add a new specialized template for your T below */
  write(vtkSmartPointer<vtkXMLWriter>::New(), mesh, filename);
}

template<>  /* inline because specialized template */

inline void VTKIOHelper::Write(vtkSmartPointer<vtkImageData> mesh, string filename) {
  write(vtkSmartPointer<vtkXMLImageDataWriter>::New(), mesh, filename);
}

template<>  /* inline because specialized template */

inline void VTKIOHelper::Write(vtkSmartPointer<vtkUnstructuredGrid> mesh, string filename) {
  if (filename.substr(filename.size() - 3) == "vtk") {  /* Convert lagcy vtk OTF */
    Verbosity::Stream(0) << "WARNING: Legacy VTK file " << filename;
    filename.replace(filename.size() - 3, 3, "vtu");
    Verbosity::Stream(0) << " saved as " << filename << std::endl;
  }
  write(vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New(), mesh, filename);
}

template<>  /* more general specialized template */

inline void VTKIOHelper::Write(vtkSmartPointer<vtkDataSet> ds, string file) {
  /* TODO: add missing vtk data types */
  switch (ds->GetDataObjectType()) {
    case VTK_UNSTRUCTURED_GRID:
      VTKIOHelper::Write<vtkUnstructuredGrid>(vtkUnstructuredGrid::SafeDownCast(ds), file);
      break;
    case VTK_IMAGE_DATA:
      VTKIOHelper::Write<vtkImageData>(vtkImageData::SafeDownCast(ds), file);
      break;

    // case VTK_POLY_DATA:
    //    VTKIOHelper::Write<vtkPolyData>(vtkPolyData::SafeDownCast(ds), file);
    //    break;
    default:
      break;
  }
}

template<class T>
vtkSmartPointer<T> VTKIOHelper::Read(string name) {
  /* This is a default template should never get called (or even compiled)
   * If you encountered an error here, check your template class T
   * or add a new specialized template for your T below */
  vtkSmartPointer<vtkXMLDataReader> reader = vtkSmartPointer<vtkXMLDataReader>::New();
  return read<T>(reader, name);
}

template<>  /* inline because specialized template */
inline vtkSmartPointer<vtkImageData> VTKIOHelper::Read(string name) {
  vtkSmartPointer<vtkXMLImageDataReader> reader = vtkSmartPointer<vtkXMLImageDataReader>::New();
  return read<vtkImageData>(reader, name);
}

template<>  /* inline because specialized template */
inline vtkSmartPointer<vtkUnstructuredGrid> VTKIOHelper::Read(string name) {
  vtkSmartPointer<vtkXMLUnstructuredGridReader> reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
  return read<vtkUnstructuredGrid>(reader, name);
}

// template<> /* a more general template */
// inline vtkSmartPointer<vtkDataSet> VTKIOHelper::Read(string filename) {
//    if (filename.rfind(".vtu") != string::npos)
//        return VTKIOHelper::Read<vtkUnstructuredGrid>(filename);
//    else if (filename.rfind(".vti") != string::npos)
//        return VTKIOHelper::Read<vtkImageData>(filename);
//    else if (filename.rfind("lat") != string::npos) {
//        vtkSmartPointer<vtkLatticeSG> l = vtkSmartPointer<vtkLatticeSG>::New();
//        l->SetLattice(filename.c_str());
//        return l;
//    }
//    else
//        throw "VTKIOHelper::Read(string): unexpected data type\n";
// }

template<class T>
vtkSmartPointer<T> VTKIOHelper::read(vtkSmartPointer<vtkXMLDataReader> reader, string name) {
  Verbosity::Stream(1) << "Reading " << name << " ... " << flush;
  reader->SetFileName(name.c_str());
  reader->Update();
  reader->GetOutputAsDataSet()->Register(reader);
  Verbosity::Stream(1) << "done." << endl;
  return T::SafeDownCast(reader->GetOutputAsDataSet());
}

void VTKIOHelper::write(vtkSmartPointer<vtkXMLWriter> writer, vtkSmartPointer<vtkDataSet> mesh, string name) {
  Verbosity::Stream(0) << "Writing " << name << " ... " << flush;
  writer->SetFileName(name.c_str());
#if VTK_MAJOR_VERSION > 5
  writer->SetInputData(mesh);
#else  // if VTK_MAJOR_VERSION > 5
  writer->SetInput(mesh);
#endif  // if VTK_MAJOR_VERSION > 5
  writer->Write();
  Verbosity::Stream(0) << "done." << endl;
}

#endif  // ifndef VTKIOHELPER_H
