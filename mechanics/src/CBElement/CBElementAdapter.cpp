/*
 * File: CBElementAdapter.cpp
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


#include "CBSolver.h"
#include "CBElementAdapter.h"



TInt CBElementAdapter::GetNumberOfNodes()
{
    return(solver_->GetNumberOfNodes());
}


TInt CBElementAdapter::GetNumberOfElements()
{
    return(solver_->GetNumberOfElements());
}

void CBElementAdapter::LinkNodesComponentsBoundaryConditionsGlobal(bool* nodesComponentsBoundaryConditionsGlobal)
{
    nodesComponentsBoundaryConditionsGlobal_ = nodesComponentsBoundaryConditionsGlobal;
    PetscInt  numTotalNodes = solver_->GetNumberOfLocalNodes() + solver_->GetNumberOfGhostNodes();
    if(nodesComponentsBoundaryConditionsLocal_ == 0)
        nodesComponentsBoundaryConditionsLocal_ = new bool[3 * numTotalNodes];
    
    for(PetscInt i = 0; i <  numTotalNodes; i++)
    {
        int n[3]={3*i,3*i+1,3*i+2};
        ISLocalToGlobalMappingApply(nodesIndicesMapping_, 3, n, n);
        nodesComponentsBoundaryConditionsLocal_[3*i]   = nodesComponentsBoundaryConditionsGlobal[n[0]];
        nodesComponentsBoundaryConditionsLocal_[3*i+1] = nodesComponentsBoundaryConditionsGlobal[n[1]];
        nodesComponentsBoundaryConditionsLocal_[3*i+2] = nodesComponentsBoundaryConditionsGlobal[n[2]];
    }
    
    areBoundaryConditionsActive_=true;
    solver_->SetNodesComponentsBoundaryConditionsGlobal(nodesComponentsBoundaryConditionsGlobal);
}

const std::vector<CBElement*>& CBElementAdapter::GetElementVector() {
    return solver_->GetElementVector();
}

void CBElementAdapter::LinkActiveStress(Vec activeStressTensor, PetscInt numActiveStressTensorComponents, PetscInt* activeStressTensorComponentsIndices)
{
    activeStressTensor_                  = activeStressTensor;
    numActiveStressTensorComponents_     = numActiveStressTensorComponents;
    activeStressTensorComponentsIndices_ = activeStressTensorComponentsIndices;
    PetscInt dummy;
    if(activeStressTensor_)
        VecGetOwnershipRange(activeStressTensor_, &activeStressTensorOwnershipLow_, &dummy);
}


// --------------------------- Coordinates -----------------------------

void CBElementAdapter::GetNodesCoords(PetscInt numNodesCoords, const PetscInt* nodesCoordsIndices, PetscScalar* nodesCoords)
{
    VecGetValues(nodes_, numNodesCoords, nodesCoordsIndices, nodesCoords);
}

void CBElementAdapter::GetGlobalNodesCoords(PetscInt numNodesCoords, const PetscInt* nodesCoordsIndices, PetscScalar* nodesCoords)
{
    VecGetValues(globalNodes_, numNodesCoords, nodesCoordsIndices, nodesCoords);
}

void CBElementAdapter::GetRefNodesCoords(PetscInt numNodesCoords, const PetscInt* nodesCoordsIndices, PetscScalar* refNodesCoords)
{
    if(refNodesCoords)
        VecGetValues(refNodes_, numNodesCoords, nodesCoordsIndices, refNodesCoords);
    else
        throw std::runtime_error("void CBElementAdapter::GetRefNodesCoords(PetscInt numNodesCoords, const PetscInt* nodesCoordsIndices, PetscScalar* refNodesCoords): refNodesCoords is zero");
}

void CBElementAdapter::SetNodesCoords(PetscInt numNodesCoords, const PetscInt* nodesCoordsIndices, PetscScalar* nodesCoords)
{
    VecSetValues(nodes_, numNodesCoords, nodesCoordsIndices, nodesCoords, INSERT_VALUES);
}

void CBElementAdapter::SetGlobalNodesCoords(PetscInt numNodesCoords, const PetscInt* nodesCoordsIndices, PetscScalar* nodesCoords)
{
    VecSetValues(globalNodes_, numNodesCoords, nodesCoordsIndices, nodesCoords, INSERT_VALUES);
}



// --------------------------- Boundary Conditions ---------------------
void CBElementAdapter::GetNodesComponentsBoundaryConditions(PetscInt numNodesCoords, const PetscInt* nodesCoordsIndices, bool* nodesComponentsBoundaryConditionsLocal)
{
    for(PetscInt i = 0; i < numNodesCoords; i++)
    {
        if(areBoundaryConditionsActive_)
            nodesComponentsBoundaryConditionsLocal[i] = nodesComponentsBoundaryConditionsLocal_[nodesCoordsIndices[i]];
        else
            nodesComponentsBoundaryConditionsLocal[i] = false;
    }
}

void CBElementAdapter::GetNodesComponentsBoundaryConditionsGlobal(PetscInt numNodesComponentsIndices, const PetscInt* nodesComponentsIndices, bool* nodesComponentsBoundaryConditionsGlobal)
{
    for(PetscInt i = 0; i < numNodesComponentsIndices; i++)
    {
        if(areBoundaryConditionsActive_)
            nodesComponentsBoundaryConditionsGlobal[i] = nodesComponentsBoundaryConditionsGlobal_[nodesComponentsIndices[i]];
        else
            nodesComponentsBoundaryConditionsGlobal[i] = false;
    }
}



// --------------------------- Nodal Forces ----------------------------

void CBElementAdapter::AddNodalForce(PetscInt i, const Vector3<PetscScalar>* force)
{
    PetscInt pos[3] = {3 * i, 3 * i + 1, 3 * i + 2};
    VecSetValuesLocal(nodalForces_, 3, pos, force->GetArray(), ADD_VALUES);
}

void CBElementAdapter::InsertNodalForce(PetscInt i, const Vector3<PetscScalar>* force)
{
    PetscInt pos[3] = {3 * i, 3 * i + 1, 3 * i + 2};
    VecSetValuesLocal(nodalForces_, 3, pos, force->GetArray(), INSERT_VALUES);
}


void CBElementAdapter::AddNodalForcesComponents(PetscInt numNodalForcesComponents, const PetscInt* nodalForcesComponentsIndices, const PetscScalar* nodalForcesComponents)
{
    VecSetValuesLocal(nodalForces_, numNodalForcesComponents, nodalForcesComponentsIndices, nodalForcesComponents, ADD_VALUES);
}

void CBElementAdapter::InsertNodalForcesComponents(PetscInt numNodalForcesComponents, const PetscInt* nodalForcesComponentsIndices, const PetscScalar* nodalForcesComponents)
{
    VecSetValuesLocal(nodalForces_, numNodalForcesComponents, nodalForcesComponentsIndices, nodalForcesComponents, INSERT_VALUES);
}

void CBElementAdapter::GetNodalForcesComponents(PetscInt numNodalForcesComponents, const PetscInt* nodalForcesComponentsIndices, PetscScalar* nodalForcesComponents)
{
    VecGetValues(nodalForces_, numNodalForcesComponents, nodalForcesComponentsIndices, nodalForcesComponents);
}

void CBElementAdapter::AddNodalForcesComponentsGlobal(PetscInt numNodalForcesComponents, const PetscInt* nodalForcesComponentsIndices, const PetscScalar* nodalForcesComponents)
{
    VecSetValues(nodalForces_, numNodalForcesComponents, nodalForcesComponentsIndices, nodalForcesComponents, ADD_VALUES);
}


// --------------------------- Nodal Forces Jacobian -------------------
void CBElementAdapter::AddNodalForcesJacobianEntries(PetscInt numRows, const PetscInt* rowsIndices, PetscInt numCols, const PetscInt* colsIndices, const PetscScalar* nodalForcesJacobianEntries)
{
    MatSetValuesLocal(nodalForcesJacobian_, numRows, rowsIndices, numCols, colsIndices, nodalForcesJacobianEntries, ADD_VALUES);
}

void CBElementAdapter::AddNodalForcesActiveStressJacobianEntries(PetscInt numRows, const PetscInt* rowsIndices, PetscInt numCols, const PetscInt* colsIndices, const PetscScalar* nodalForcesJacobianEntries)
{
    MatSetValues(nodalForcesActiveStressJacobian_, numRows, rowsIndices, numCols, colsIndices, nodalForcesJacobianEntries, ADD_VALUES);
}

void CBElementAdapter::AddNodalForcesActiveStressTensorAndFiberOrientationJacobianEntries(PetscInt numRows, const PetscInt* rowsIndices, PetscInt numCols, const PetscInt* colsIndices, const PetscScalar* nodalForcesJacobianEntries)
{
    MatSetValues(nodalForcesActiveStressTensorAndFiberOrientationJacobian_, numRows, rowsIndices, numCols, colsIndices, nodalForcesJacobianEntries, ADD_VALUES);
}
void CBElementAdapter::AddNodalForcesJacobianEntriesGlobal(const PetscInt numRows, const PetscInt* rowsIndices, const PetscInt numCols, const PetscInt* colsIndices, const PetscScalar* nodalForcesJacobianEntries)
{
    MatSetValues(nodalForcesJacobian_, numRows, rowsIndices, numCols, colsIndices, nodalForcesJacobianEntries, ADD_VALUES);
}

void CBElementAdapter::InsertNodalForcesComponentsGlobal(PetscInt numNodalForcesComponents, const PetscInt* nodalForcesComponentsIndices, const PetscScalar* nodalForcesComponents)
{
    VecSetValues(nodalForces_, numNodalForcesComponents, nodalForcesComponentsIndices, nodalForcesComponents, INSERT_VALUES);
}

void CBElementAdapter::InsertNodalForcesJacobianEntries(PetscInt numRows, const PetscInt* rowsIndices, PetscInt numCols, const PetscInt* colsIndices, const PetscScalar* nodalForcesJacobianEntries)
{
    MatSetValuesLocal(nodalForcesJacobian_, numRows, rowsIndices, numCols, colsIndices, nodalForcesJacobianEntries, INSERT_VALUES);
}

void CBElementAdapter::InsertNodalForcesJacobianEntriesGlobal(const PetscInt numRows, const PetscInt* rowsIndices, const PetscInt numCols, const PetscInt* colsIndices, const PetscScalar* nodalForcesJacobianEntries)
{
    MatSetValues(nodalForcesJacobian_, numRows, rowsIndices, numCols, colsIndices, nodalForcesJacobianEntries, INSERT_VALUES);
}


// ---------------------------- Mass Matrix ----------------------------
void CBElementAdapter::SetMassMatrixEntries(PetscInt numRows, const PetscInt* rowsIndices, PetscInt numCols, const PetscInt* colsIndices, PetscScalar* massMatrixEntries)
{
    MatSetValuesLocal(massMatrix_, numRows, rowsIndices, numCols, colsIndices, massMatrixEntries, ADD_VALUES);
}

// ---------------------------- Stiffness Matrix ----------------------------
void CBElementAdapter::SetStiffnessMatrixEntries(PetscInt numRows, const PetscInt* rowsIndices, PetscInt numCols, const PetscInt* colsIndices, PetscScalar* stiffnessMatrixEntries)
{
    // Assemble the global Stiffness Matrix K, which is the sum of all element stiffness matrices Ke over all elements
    MatSetValuesLocal(stiffnessMatrix_, numRows, rowsIndices, numCols, colsIndices, stiffnessMatrixEntries, ADD_VALUES);
}

// ---------------------------- Laplacian Matrix ----------------------------
void CBElementAdapter::AddLaplacianEntriesGlobal(PetscInt numRows, const PetscInt* rowsIndices, PetscInt numCols, const PetscInt* colsIndices, PetscScalar* laplacianEntries)
{
    MatSetValues(laplacian_, numRows, rowsIndices, numCols, colsIndices, laplacianEntries, ADD_VALUES);
}

// ---------------------------- Active Stress Tensor -----------------
void CBElementAdapter::GetActiveStressTensor(PetscInt elementIndex, Matrix3<PetscScalar>& stressTensor)
{
    PetscScalar activeStressTensorComponents[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
    if(activeStressTensor_)
    {
        PetscInt indices[9];
        for(unsigned int i = 0; i < numActiveStressTensorComponents_; i++)
            indices[i] = numActiveStressTensorComponents_ * (elementIndex + activeStressTensorOwnershipLow_) + activeStressTensorComponentsIndices_[i];
        VecGetValues(activeStressTensor_, numActiveStressTensorComponents_, indices, activeStressTensorComponents);
    }
    stressTensor.SetArray(activeStressTensorComponents);
}
