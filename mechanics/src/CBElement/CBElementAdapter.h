/*
 *  CBElementAdapter.h
 *  CardioMechanics
 *
 *  Created by Thomas Fritz on 15.06.10.
 *  Copyright 2010 Institute for Biomedical Engineering, Karlsruhe Institute of Technology (KIT). All rights reserved.
 *
 */

#ifndef CB_ELEMENT_ADAPTER_H
#define CB_ELEMENT_ADAPTER_H

#include <vector>

#include <petscksp.h>

#include "Matrix3.h"

#include "DCType.h"

class CBSolver;
class CBElement;

using namespace math_pack;


class CBElementAdapter
{
public:
    CBElementAdapter() : activeStressTensor_(0), numActiveStressTensorComponents_(0), activeStressTensorComponentsIndices_(0), activeStressTensorOwnershipLow_(0), nodesComponentsBoundaryConditionsLocal_(0), nodesComponentsBoundaryConditionsGlobal_(0), areBoundaryConditionsActive_(true)
    {}
    
    virtual ~CBElementAdapter(){fdEpsilon_ = 0;
        solver_    = 0; }
    
    virtual void Init(){}
    void SetSolver(CBSolver* solver){solver_ = solver; }
    CBSolver* GetSolver(){return(solver_); }
    
    TInt GetNumberOfNodes();
    TInt GetNumberOfElements();
    
    Vec GetNodes(){return nodes_;} 
    
    void AddToForces(Vec f){VecAXPY(nodalForces_,1,f);}
    void AddToNodalForcesJacobianDiagonal(Vec d){MatDiagonalSet(nodalForcesJacobian_,d,ADD_VALUES);}
    
    void GetActiveStressTensor(PetscInt elementIndex, Matrix3<PetscScalar>& stressTensor);
    
    //! numNodesCoords = amount of nodes components indices = 3* amount of nodes since each node has 3 components.
    //!
    //! nodesCoordsIndices = indices of the nodes components. e.g. If a node has the local index 90 ist components have the components indices 3*90+0, 3*90+1, 3*90+2
    //! see info about node and component numbering !!!
    //! nodesCoords = Array with size sizeof(T)*numNodesCoords
    //! Example: 4 nodes with local indices 10 15 30 90 -> nodesCoordsIndices = {30,31,32,45,46,47,90,91,92,270,271,271}
    
    void GetNodesCoords(PetscInt numNodesCoords, const PetscInt* nodesCoordsIndices, PetscScalar* nodesCoords);
    void GetGlobalNodesCoords(PetscInt numNodesCoords, const PetscInt* nodesCoordsIndices, PetscScalar* nodesCoords);
    void GetRefNodesCoords(PetscInt numNodesCoords, const PetscInt* nodesCoordsIndices, PetscScalar* refNodesCoords);
    
    void SetNodesCoords(PetscInt numNodesCoords, const PetscInt* nodesCoordsIndices, PetscScalar* nodesCoords);
    void SetGlobalNodesCoords(PetscInt numNodesCoords, const PetscInt* nodesCoordsIndices, PetscScalar* nodesCoords);
    void SetRefNodesCoords(PetscInt numNodesCoords, const PetscInt* nodesCoordsIndices, PetscScalar* nodesCoords);
    
    
    void AddNodalForce(PetscInt i, const Vector3<PetscScalar>* force);
    void InsertNodalForce(PetscInt i, const Vector3<PetscScalar>* force);
    
    //! The indexing of InsertNodalForcesComponents is different to InsertNodalForces since it uses nodes Indices and the other used nodes forces components indices !!!!!
    //!
    //! numNodalForcesComponents: number of vector components to set = 3 * nodes per element
    //! nodalForcesComponentsIndices: indices of the vector components to set
    //! nodalForcesComponents: components of the forces vectors to set.
    //!
    //! Example: The element has 4 nodes and 4 force vectors for each node -> number of force components = 4 * 3 = 12,
    //!
    //! nodalForcesComponentsIndices = {3 * nodesIndices[0], 3 * nodesIndices[0]+1, 3 * nodesIndices[0]+2,
    //! 3 * nodesIndices[1] ,3 * nodesIndices[1]+1, 3 * nodesIndices[1]+2,
    //! 3 * nodesIndices[2] ,3 * nodesIndices[2]+1, 3 * nodesIndices[2]+2,
    //! 3 * nodesIndices[3] ,3 * nodesIndices[3]+1, 3 * nodesIndices[3]+2}
    
    void AddNodalForcesComponents(PetscInt numNodalForcesComponents, const PetscInt* nodalForcesComponentsIndices, const PetscScalar* nodalForcesComponents);
    
    void InsertNodalForcesComponents(PetscInt numNodalForcesComponents, const PetscInt* nodalForcesComponentsIndices, const PetscScalar* nodalForcesComponents);
    
    void GetNodalForcesComponents(PetscInt numNodalForcesComponents, const PetscInt* nodalForcesComponentsIndices, PetscScalar* nodalForcesComponents);
    
    void AddNodalForcesComponentsGlobal(PetscInt numNodalForcesComponents, const PetscInt* nodalForcesComponentsIndices, const PetscScalar* nodalForcesComponents);
    
    void InsertNodalForcesComponentsGlobal(PetscInt numNodalForcesComponents, const PetscInt* nodalForcesComponentsIndices, const PetscScalar* nodalForcesComponents);
    
    void AddNodalForcesJacobianEntries(PetscInt numRows, const PetscInt* rowsIndices, PetscInt numCols, const PetscInt* colsIndices, const PetscScalar* nodalForcesJacobianEntries);
    void AddNodalForcesActiveStressJacobianEntries(PetscInt numRows, const PetscInt* rowsIndices, PetscInt numCols, const PetscInt* colsIndices, const PetscScalar* nodalForcesJacobianEntries);
    void AddNodalForcesActiveStressTensorAndFiberOrientationJacobianEntries(PetscInt numRows, const PetscInt* rowsIndices, PetscInt numCols, const PetscInt* colsIndices, const PetscScalar* nodalForcesJacobianEntries);
    
    void InsertNodalForcesJacobianEntries(PetscInt numRows, const PetscInt* rowsIndices, PetscInt numCols, const PetscInt* colsIndices, const PetscScalar* nodalForcesJacobianEntries);
    
    void AddNodalForcesJacobianEntriesGlobal(const PetscInt numRows, const PetscInt* rowsIndices, const PetscInt numCols, const PetscInt* colsIndices, const PetscScalar* nodalForcesJacobianEntries);
    
    void InsertNodalForcesJacobianEntriesGlobal(const PetscInt numRows, const PetscInt* rowsIndices, const PetscInt numCols, const PetscInt* colsIndices, const PetscScalar* nodalForcesJacobianEntries);
    
    void SetMassMatrixEntries(PetscInt numRows, const PetscInt* rowsIndices, PetscInt numCols, const PetscInt* colsIndices, PetscScalar* massMatrixEntries);
    
    void SetStiffnessMatrixEntries(PetscInt numRows, const PetscInt* rowsIndices, PetscInt numCols, const PetscInt* colsIndices, PetscScalar* stiffnessMatrixEntries);
    
    void AddLaplacianEntriesGlobal(PetscInt numRows, const PetscInt* rowsIndices, PetscInt numCols, const PetscInt* colsIndices, PetscScalar* laplacianEntries);
    void GetNodesComponentsBoundaryConditions(PetscInt numNodesCoords, const PetscInt* nodesCoordsIndices, bool* nodesComponentsBoundaryConditions);
    void GetNodesComponentsBoundaryConditionsGlobal(PetscInt numNodesCoords, const PetscInt* nodesCoordsIndices, bool* globalNodesComponentsBoundaryConditions);
    void ApplyLocalToGlobalMapping(PetscInt* nodesCoordsIndices, PetscInt numIndices){ISLocalToGlobalMappingApply(nodesIndicesMapping_, numIndices, nodesCoordsIndices, nodesCoordsIndices); }
    void ApplyGlobalToLocalMapping(PetscInt* nodesCoordsIndices, PetscInt numIndices) { ISGlobalToLocalMappingApply(nodesIndicesMapping_, IS_GTOLM_MASK, numIndices, nodesCoordsIndices, &numIndices, nodesCoordsIndices); } // index is -1 if not on local process
    void ApplyGlobalToLocalMappingNonGhosted(PetscInt* nodesCoordsIndices, PetscInt numIndices) { ISGlobalToLocalMappingApply(nodesIndicesMappingNonGhosted_, IS_GTOLM_MASK, numIndices, nodesCoordsIndices, &numIndices, nodesCoordsIndices); } // index is -1 if not on local process
    
    
    void LinkNodes(Vec nodes){nodes_ = nodes;}
    void LinkGlobalNodes(Vec nodes){globalNodes_ = nodes;}
    void LinkRefNodes(Vec refNodes){refNodes_ = refNodes;}
    void LinkNodesComponentsBoundaryConditionsGlobal(bool* nodesComponentsBoundaryConditionsGlobal); // Sets nodesComponentsBoundaryConditionsLocal automatically
    void LinkNodalForces(Vec nodalForces){nodalForces_ = nodalForces; }
    void LinkMassMatrix(Mat massMatrix){massMatrix_ = massMatrix; }
    void LinkStiffnessMatrix(Mat stiffnessMatrix){stiffnessMatrix_ = stiffnessMatrix; }
    void LinkNodalForcesJacobian(Mat nodalForcesJacobian){nodalForcesJacobian_ = nodalForcesJacobian; }
    void LinkNodalForcesActiveStressJacobian(Mat nodalForcesActiveStressJacobian){nodalForcesActiveStressJacobian_ = nodalForcesActiveStressJacobian; }
    void LinkNodalForcesActiveStressTensorAndFiberOrientationJacobian(Mat nodalForcesActiveStressTensorAndFiberOrientationJacobian){nodalForcesActiveStressTensorAndFiberOrientationJacobian_ = nodalForcesActiveStressTensorAndFiberOrientationJacobian; }
    void LinkActiveStress(Vec activeStressTensor, PetscInt numActiveStressTensorComponents, PetscInt* activeStressTensorComponentsIndices);
    void LinkLaplacian(Mat laplacian){laplacian_ = laplacian; }
    
    void LinkLocalToGlobalMapping(ISLocalToGlobalMapping nodesIndicesMapping){nodesIndicesMapping_ = nodesIndicesMapping; }
    void LinkLocalToGlobalMappingNonGhosted(ISLocalToGlobalMapping nodesIndicesMapping){nodesIndicesMappingNonGhosted_ = nodesIndicesMapping; }
    
    void ActivateBoundaryConditions(){areBoundaryConditionsActive_=true;}
    void DeactivateBoundaryConditions(){areBoundaryConditionsActive_=false;}
    bool ActiveBoundaryConditions(){return areBoundaryConditionsActive_;}
    PetscInt GetLocalNodesFrom();
    PetscInt GetLocalNodesTo();
    
    void SetFiniteDifferencesEpsilon(TFloat fdEpsilon){fdEpsilon_ = fdEpsilon; }
    TFloat GetFiniteDifferencesEpsilon(){return(fdEpsilon_); }
    const std::vector<CBElement*>& GetElementVector();
    
protected:
private:
    
    TFloat    fdEpsilon_;
    CBSolver* solver_;
    
    Vec                    nodes_;
    Vec                    globalNodes_;
    Vec                    refNodes_;
    
    Vec                    nodalForces_;
    Mat                    nodalForcesJacobian_;
    Mat                    nodalForcesActiveStressJacobian_;
    Mat                    nodalForcesActiveStressTensorAndFiberOrientationJacobian_;
    
    Mat                    massMatrix_;
    Mat                    stiffnessMatrix_;
    Mat                    laplacian_;
    
    Vec                    activeStressTensor_;
    PetscInt               numActiveStressTensorComponents_;
    PetscInt*              activeStressTensorComponentsIndices_;
    PetscInt               activeStressTensorOwnershipLow_;
    bool*                  nodesComponentsBoundaryConditionsLocal_;
    bool*                  nodesComponentsBoundaryConditionsGlobal_;
    bool                   areBoundaryConditionsActive_;
    ISLocalToGlobalMapping nodesIndicesMapping_;
    ISLocalToGlobalMapping nodesIndicesMappingNonGhosted_;
};

#endif
