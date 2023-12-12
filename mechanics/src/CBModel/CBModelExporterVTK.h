/*
 * File: CBModelExporterVTK.h
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


#ifndef CB_MODEL_EXPORTER_VTK
#define CB_MODEL_EXPORTER_VTK

#include <sys/types.h>
#include <sys/stat.h>

#include <map>
#include <set>
#include <vector>
#include <string>
#include <sstream>

#include <vtkCellArray.h>
#include <vtkPoints.h>
#include <vtkSmartPointer.h>
#include <vtkLine.h>
#include <vtkTriangle.h>
#include <vtkPolyData.h>
#include <vtkDataSetMapper.h>
#include <vtkPlane.h>
#include <vtkTetra.h>
#include <vtkQuadraticTetra.h>
#include <vtkQuadraticTriangle.h>
#include <vtkUnstructuredGrid.h>
#include <vtkIntArray.h>
#include <vtkDoubleArray.h>
#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkProperty.h>
#include <vtkXMLUnstructuredGridWriter.h>

#include "Matrix3.h"

#include "CBModelExporter.h"


#include "ParameterMap.h"

class CBModelExporterVTK : public CBModelExporter
{
public:
    CBModelExporterVTK(CBModel* model, ParameterMap* parameters) : CBModelExporter(model, parameters)
    {
        InitPvdFile();
        InitTetgen();
    }
    ~CBModelExporterVTK(){}
    bool ExportModel();
    vtkSmartPointer<vtkUnstructuredGrid> GetMesh();
protected:
private:
    void InitPvdFile();
    CBModelExporterVTK();
    typedef CBModelExporter   Base;
    
    void InitTetgen();
    void ExportTetgenNode();
    void ExportTetgenBases();
};

#endif
