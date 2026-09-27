/*
 * File: CBacCELLerate.h
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


#pragma once

#include <vtkCellLocator.h>
#include <vtkCellArray.h>
#include <vtkPoints.h>
#include <vtkSmartPointer.h>
#include <vtkLine.h>
#include <vtkTriangle.h>
#include <vtkPolyData.h>
#include <vtkDataSetMapper.h>
#include <vtkPlane.h>
#include <vtkTetra.h>
#include <vtkUnstructuredGrid.h>
#include <vtkIntArray.h>
#include <vtkDoubleArray.h>
#include <vtkDataArray.h>
#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkProperty.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkMeshQuality.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkFillHolesFilter.h>
#include <vtkPolyDataNormals.h>
#include <vtkFeatureEdges.h>
#include <vtkCleanPolyData.h>
#include <vtkPolyDataConnectivityFilter.h>
#include <vtkAppendPolyData.h>
#include <vtkStripper.h>
#include <vtkTriangleFilter.h>
#include <vtkPointLocator.h>
#include <vtkSelectEnclosedPoints.h>
#include <vtkMath.h>
#include <vtkLine.h>

#include "CBacCELLerate.h"
#include "CBDataCtrl.h"
#include "Matrix4.h"
#include "CBSolver.h"
#include <kaPoint.h>
#include "CBSolverPlugin.h"
#include "CBDataFromFile.h"
#include "CBDataCtrl.h"
#include "acCELLerate.h"
#include "acltCellModel.h"
#include "acltTime.h"
#include <PETScLSE.h>
#include <Material.h>
#include <array>
#include <vector>

class CBacCELLerate : public CBSolverPlugin {
public:
    ~CBacCELLerate() {
        VecScatterDestroy(&localNodesScatter_);
        VecDestroy(&localNodes_);
    }
    
    std::string GetName() override { return std::string("acCELLerate"); }
    
    void Init() override;
    void Prepare() override;
    void Apply(TFloat time) override;
    void Export(TFloat) override;
    void StepBack() override;
    
    CBStatus GetStatus() override { return status_; }
    
protected:
private:
    void UpdateNodes(bool useReferenceNodes = false);
    void UpdateStretch();
    void UpdateVelocity();
    void AssembleMatrix();
    void InitPvdFile();
    void InitMapping();
    void InitParameters();
    void InitMesh();
    void InitPetscVec();
    void InitLocalMaps();
    void CalcShapeFunctionDeriv();
    void CalcDeformationTensor(vtkIdType cellID, const Vector3<TFloat> *nodesCoords, Matrix3<TFloat> &deformationTensor);
    
    Matrix3<TFloat> GetBasisAtCell(vtkIdType cellID) {return Q_[cellID];}
    
    void GetFibersFromMechanicsMesh();
    void ApplySpatialSortPCA();
    void DetermineElementRanges();
    
    
    TFloat stopTime_ = 0;
    float stepBackTime_ = 0;
    bool stepBack_    = false;
    CBStatus status_ = CBStatus::WAITING;
    
    /// Plugin parameters
    int NumQP_;
    float offsetTime_ = 0;
    bool export_ = true;
    bool constStretchRate_ = true;
    std::vector<int> NumConele_;
    std::string MEF_;
    std::string MaterialFileName_;
    std::string resPreFix_;
    std::string resFolder_;
    std::string pvdFilename_;
    std::string accprojectFile_;
    std::vector<TInt> materialCoupling_;
    std::vector<TInt> priorityVector_;
    
    /// Global vectors
    Vec stretchVecF_;
    Vec velocityVec_;
    Vec accNodes_;
    Vec stepbackForce_;
    Vec timestepForce_;

    /// Global Matrices
    Mat sysMatrix_ = NULL;
    Mat massMatrix_ = NULL;
    
    /// acCELLerate instance
    acCELLerate *act_ = NULL;
    
    /// List of mech solid Elements
    std::vector<CBElementSolid *> solidElements_;
    
    /// acMesh_ point expressed in the shape functions of the solid element containing (or closest to) it
    struct PointMapping {
        PetscInt point;
        CBElementSolid *element;
        Vector4<TFloat> shapeFun;
    };

    /// Points owned by this process, sorted by point index. Each point is owned by exactly one process
    /// so that UpdateNodes/UpdateStretch set every entry of the global vectors once.
    std::vector<PointMapping> pointMappings_;

    /// Quadrature point of a solid element expressed in the shape functions of the closest acMesh_ cell
    struct QPMapping {
        PetscInt points[4];
        Vector4<TFloat> shapeFun;
    };

    /// Indexed by position in solidElements_ * NumQP_ + quadrature point; unused for uncoupled materials
    std::vector<QPMapping> qpMappings_;

    /// Coupled solid element as sent to process 0 for the mapping
    struct MappingElement {
        PetscInt eIdx;      // position in solidElements_ of the sending process
        double coords[12];  // vertices [m]
    };

    /// acMesh_ point mapped to the solid element at position eIdx in solidElements_ of the owning process
    struct MappingPoint {
        PetscInt point;
        PetscInt eIdx;
        double shapeFun[4];
    };

    /// Quadrature point mapped to an acMesh_ cell
    struct MappingQP {
        PetscInt points[4];
        double shapeFun[4];
    };

    /// Computes on process 0 the mapping of the coupled solid elements of all processes; elements of process r
    /// are elements[elementDispls[r]:elementDispls[r+1]]. points are grouped by owning process, qps follow elements.
    void MapElements(const std::vector<MappingElement> &elements, const std::vector<int> &elementDispls,
                     std::vector<MappingPoint> &points, std::vector<int> &pointCounts, std::vector<MappingQP> &qps);

    /// Per local acMesh_ cell: index 0-11 contains dN/dX; index 12 contains tet volume
    std::vector<std::array<double, 13>> dNdX_;

    /// Per local acMesh_ cell: reference basis
    std::vector<Matrix3<TFloat>> Q_;

    /// Per local acMesh_ cell: vertices as indices into localPoints_
    std::vector<std::array<PetscInt, 4>> localCells_;

    /// Per local acMesh_ cell: material
    std::vector<int> localMaterials_;

    /// Global indices of the acMesh_ points used by local cells, sorted
    std::vector<PetscInt> localPoints_;

    /// Coordinates [mm] of localPoints_
    std::vector<Vector3<TFloat>> localCoords_;

    /// Gathers the coordinates of localPoints_ from accNodes_
    VecScatter localNodesScatter_ = NULL;
    Vec localNodes_ = NULL;

    /// for node permutation using pca
    std::vector<TInt> backwardMapping_;
    std::vector<TInt> forwardMapping_;
    
    /// parallel element layout
    std::vector<PetscInt> elementRanges_;
    PetscInt localElementsFrom_ = 0;
    PetscInt localElementsTo_ = 0;
    PetscInt numLocalElements_ = 0;
    
    TInt mpirank_;
    TInt mpisize_;
    
    PetscInt Istart, Iend;
    
    /// Mesh properties
    vtkIdType nPoints_;
    vtkIdType nCells_;
    std::vector<CBElementSolid *> elements_;

    /// Full acMesh_, only held by process 0 during Init(); afterwards each process keeps its local cells in localCells_
    vtkSmartPointer<vtkUnstructuredGrid> acMesh_;
    vtkSmartPointer<vtkDoubleArray> acMeshMaterials_;
    vtkSmartPointer<vtkDataArray> acMeshFiberValues_;
    vtkSmartPointer<vtkDataArray> acMeshSheetValues_;
    vtkSmartPointer<vtkDataArray> acMeshNormalValues_;
};
