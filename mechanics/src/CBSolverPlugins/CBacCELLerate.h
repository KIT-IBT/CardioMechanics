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
#include <map>

class CBacCELLerate : public CBSolverPlugin {
public:
    ~CBacCELLerate() {}
    
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
    Matrix3<TFloat> GetInitialBasisAtCell(vtkIdType cellID);
    
    Matrix3<TFloat> GetBasisAtCell(vtkIdType cellID) {return Q_[cellID];}
    
    void GetFibersFromMechanicsMesh();
    void ApplySpatialSortPCA();
    void DetermineElementRanges();
    
    
    TFloat stopTime_ = 0;
    float stepBackTime_ = 0;
    bool stepBack_    = false;
    int counter_ = 0;
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
    
    /// Global vectors
    Vec stretchVecF_;
    Vec velocityVec_;
    Vec accNodes_;
    Vec stepbackForce_;
    Vec timestepForce_;
    Vec ClosestEleShapeFun_;
    Vec deformation_;
    
    /// Global Matrices
    Mat sysMatrix_ = NULL;
    Mat massMatrix_ = NULL;
    
    /// acCELLerate instance
    acCELLerate *act_ = NULL;
    
    /// List of mech solid Elements
    std::vector<CBElementSolid *> solidElements_;
    
    /// List that maps accNodes_ to solidElements_
    std::map<PetscInt, CBElementSolid *> eleList_;
    
    /// List that maps accNodes_ to shape functions within closest solidElement
    /// 1: numPoints*5
    /// #1 is distance to ele centroid; #2-5: Shape fun
    std::map<PetscInt, vector<double>> nearC_;
    
    /// List that maps Gauss point of solidElements_ to accNodes_
    std::map<PetscInt, std::vector<PetscInt>> nearP_;
    
    /// List that maps Gauss point to shape functions within closest acc ele
    std::map<PetscInt, Vector4<TFloat>> ForceSF_;
    
    /// List that maps acMesh_ elements to the respective dN/dX
    /// index 0-11 contains dNdX; index 12 contains tet volume
    std::map<vtkIdType, vector<double>> dNdX_;
    
    /// List that maps acMesh_ elements to current deformation tensor
    std::map<vtkIdType, Matrix3<TFloat>> F_;
    
    /// List that maps acMesh_ elements to the reference basis
    std::map<PetscInt, Matrix3<TFloat>> Q_;
    
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
    vtkSmartPointer<vtkUnstructuredGrid> acMesh_;
    vtkSmartPointer<vtkDoubleArray> acMeshMaterials_;
    vtkSmartPointer<vtkDoubleArray> acMeshFiberValues_;
    vtkSmartPointer<vtkDoubleArray> acMeshSheetValues_;
    vtkSmartPointer<vtkDoubleArray> acMeshNormalValues_;
};
