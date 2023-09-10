/*
 *  CBElementSolidT4.cpp
 *  CardioMechanics
 *
 *  Created by Thomas Fritz on 13.09.11.
 *  Copyright 2011 Institute for Biomedical Engineering, Karlsruhe Institute of Technology (KIT). All rights reserved.
 *
 */

#include "CBElementSolidT4.h"
#include "CBElementAdapter.h"
#include "CBSolver.h"

void CBElementSolidT4::SetNodeIndex(unsigned int i, TInt j) {
    if (i > 3)
        throw std::runtime_error("CBElementSolidT4<T>::SetNode(unsigned int i,TInt j) -> i = out of range");
    else
        nodesIndices_[i] = j;
}

CBElementSolidT4::CBElementSolidT4(CBElementSolidT4 &other) : CBElementSolid(other) {
    nodesIndices_ = other.nodesIndices_;
    dNdX_ = other.dNdX_;
    basisAtQuadraturePoint_ = other.basisAtQuadraturePoint_;
}

CBElement *CBElementSolidT4::New() {
    return new CBElementSolidT4;
}

CBElement *CBElementSolidT4::Clone() {
    return new CBElementSolidT4(*this);
}

TInt CBElementSolidT4::GetNodeIndex(unsigned int i) {
    if (i > 3)
        throw std::runtime_error("CBElementSolidT4<T>::SetNode(unsigned int i,TInt j) -> i = out of range");
    else
        return nodesIndices_[i];
}

void CBElementSolidT4::SetBasisAtQuadraturePoint(int i, const Matrix3<TFloat> &basis) {
    if (i == 0) {
        basisAtQuadraturePoint_ = basis;
    } else {
        throw std::runtime_error("void CBElementSolidT4::void SetBasisAtQuadraturePoint(int i, const Matrix3<TFloat>& basis): i out of range");
    }
}

Matrix3<TFloat> *CBElementSolidT4::GetBasisAtQuadraturePoint(int i) {
    {
        if (i == 0) {
            return &basisAtQuadraturePoint_;
        } else {
            throw std::runtime_error("Matrix3<TFloat>* CBElementSolidT4::void GetBasisAtQuadraturePoint(int i, const Matrix3<TFloat>& basis): i out of range");
        }
    }
}

void CBElementSolidT4::GetNodesCoordsIndices(int *nodesCoordsIndices) {
    for (unsigned int i = 0; i < 4; i++) {
        nodesCoordsIndices[3*i]   = 3*nodesIndices_[i];
        nodesCoordsIndices[3*i+1] = 3*nodesIndices_[i]+1;
        nodesCoordsIndices[3*i+2] = 3*nodesIndices_[i]+2;
    }
}

void CBElementSolidT4::CalcShapeFunctionDerivatives(TFloat l1, TFloat l2, TFloat l3, TFloat l4, TFloat *dNdX,
                                                    bool useReferenceNodes) {
    TFloat nodesCoords[12];
    TInt   nodesCoordsIndices[12];
    
    
    GetNodesCoordsIndices(nodesCoordsIndices);
    
    if (useReferenceNodes)
        Base::adapter_->GetRefNodesCoords(12, nodesCoordsIndices, nodesCoords);
    else
        Base::adapter_->GetNodesCoords(12, nodesCoordsIndices, nodesCoords);
    
    
    TFloat z43 = (nodesCoords[11] - nodesCoords[8]);
    TFloat z42 = (nodesCoords[11] - nodesCoords[5]);
    TFloat z41 = (nodesCoords[11] - nodesCoords[2]);
    TFloat z32 = (nodesCoords[8] - nodesCoords[5]);
    TFloat z31 = (nodesCoords[8] - nodesCoords[2]);
    TFloat z21 = (nodesCoords[5] - nodesCoords[2]);
    
    TFloat y43 = (nodesCoords[10] - nodesCoords[7]);
    TFloat y42 = (nodesCoords[10] - nodesCoords[4]);
    TFloat y41 = (nodesCoords[10] - nodesCoords[1]);
    TFloat y32 = (nodesCoords[7] - nodesCoords[4]);
    TFloat y31 = (nodesCoords[7] - nodesCoords[1]);
    TFloat y21 = (nodesCoords[4] - nodesCoords[1]);
    
    
    TFloat x41 = (nodesCoords[9] - nodesCoords[0]);
    TFloat x31 = (nodesCoords[6] - nodesCoords[0]);
    TFloat x21 = (nodesCoords[3] - nodesCoords[0]);
    
    detJ_ = x21 * (y31 * z41 - y41 * z31) + y21 * (x41 * z31 - x31 * z41) + z21 * (x31 * y41 - x41 * y31);
    
    dNdX[0] = 1.0 / detJ_ *(nodesCoords[4] * z43 - nodesCoords[7] * z42 + nodesCoords[10] * z32);
    dNdX[3] = 1.0 / detJ_ *(-nodesCoords[1] * z43 + nodesCoords[7] * z41 - nodesCoords[10] * z31);
    dNdX[6] = 1.0 / detJ_ *(nodesCoords[1] * z42 - nodesCoords[4] * z41 + nodesCoords[10] * z21);
    dNdX[9] = 1.0 / detJ_ *(-nodesCoords[1] * z32 + nodesCoords[4] * z31 - nodesCoords[7] * z21);
    
    dNdX[1] = 1.0 / detJ_ *(-nodesCoords[3] * z43 + nodesCoords[6] * z42 - nodesCoords[9] * z32);
    dNdX[4] = 1.0 / detJ_ *(nodesCoords[0] * z43 - nodesCoords[6] * z41 + nodesCoords[9] * z31);
    dNdX[7] = 1.0 / detJ_ *(-nodesCoords[0] * z42 + nodesCoords[3] * z41 - nodesCoords[9] * z21);
    dNdX[10] = 1.0 / detJ_ *(nodesCoords[0] * z32 - nodesCoords[3] * z31 + nodesCoords[6] * z21);
    
    dNdX[2] = 1.0 / detJ_ *(nodesCoords[3] * y43 - nodesCoords[6] * y42 + nodesCoords[9] * y32);
    dNdX[5] = 1.0 / detJ_ *(-nodesCoords[0] * y43 + nodesCoords[6] * y41 - nodesCoords[9] * y31);
    dNdX[8] = 1.0 / detJ_ *(nodesCoords[0] * y42 - nodesCoords[3] * y41 + nodesCoords[9] * y21);
    dNdX[11] = 1.0 / detJ_ *(-nodesCoords[0] * y32 + nodesCoords[3] * y31 - nodesCoords[6] * y21);
} // CBElementSolidT4::CalcShapeFunctionDerivatives

void CBElementSolidT4::CalcShapeFunctionsDerivatives() {
    // the result is independent from l1,l2,l3,l4 since dNdX are constant, however for consistency i have chosen the centroid of the element
    CalcShapeFunctionDerivatives(0.25, 0.25, 0.25, 0.25, dNdX_.data());
    Ancestor::initialVolume_ = GetVolume();
}

void CBElementSolidT4::CalcDeformationTensorWithLocalBasis(const TFloat *nodesCoords,
                                                           Matrix3<TFloat> &deformationTensor) {
    TFloat *f = deformationTensor.GetArray();
    
    for (unsigned int i = 0; i < 3; i++) {
        f[i]   = dNdX_[i] * nodesCoords[0] + dNdX_[i+3] * nodesCoords[3] + dNdX_[i+6] * nodesCoords[6] + dNdX_[i+9] *
        nodesCoords[9];
        f[3+i] = dNdX_[i] * nodesCoords[1] + dNdX_[i+3] * nodesCoords[4] + dNdX_[i+6] * nodesCoords[7] + dNdX_[i+9] *
        nodesCoords[10];
        f[6+i] = dNdX_[i] * nodesCoords[2] + dNdX_[i+3] * nodesCoords[5] + dNdX_[i+6] * nodesCoords[8] + dNdX_[i+9] *
        nodesCoords[11];
    }
    deformationTensor = GetBasisAtQuadraturePoint(0)->GetTranspose() * deformationTensor *
    GetBasisAtQuadraturePoint(0)->GetInverse().GetTranspose();
}

// TODO: remove activeStress from function arguments, communication should now rather happen over the TensionModel (e.g. CBFileTension)
CBStatus CBElementSolidT4::CalcNodalForcesHelperFunction(const TFloat *nodesCoords, const bool *boundaryConditions, TFloat *forces) {
    CBStatus rc;
    
    Matrix3<TFloat> deformationTensor;
    Matrix3<TFloat> as     = {0, 0, 0, 0, 0, 0, 0, 0, 0};
    Matrix3<TFloat> stress = {0, 0, 0, 0, 0, 0, 0, 0, 0};
    
    CalcDeformationTensorWithLocalBasis(nodesCoords, deformationTensor);
    
    // compute passive stress component
    rc = Base::material_->GetConstitutiveModel()->CalcPK2Stress(deformationTensor, stress);
    
    // print corrupt element messages from each processor
    if (rc == CBStatus::CORRUPT_ELEMENT) {
        std::cout << "SolidT4::CalcNodalForcesHelperFunction(): Element with index " << index_ <<
        " is corrupt. Det = " << deformationTensor.Det() << std::endl;
    }
    
    if (rc != CBStatus::SUCCESS)
        return rc;
    
    TFloat time = adapter_->GetSolver()->GetTiming().GetCurrentTime();
    
    // compute active stress
    as = Base::tensionModel_->CalcActiveStress(deformationTensor, time);
    stress += as;
    
    // convert PK2 stress into nominal stress with respect to the local coordinate system aligned with the fibres
    stress = stress * deformationTensor.GetTranspose();
    
    // transform nominal stress back into the global coordinate system
    stress = GetBasisAtQuadraturePoint(0)->GetInverse().GetTranspose() * stress * GetBasisAtQuadraturePoint(0)->GetTranspose();
    
    // Calculate force contribution from all nodes with 1st order quadrature rule:
    // W_i = 1
    // tetCoords = {a,a,a,a} (a = 0.25)
    // f_i = V * f_j(a,a,a,a)
    for (int k = 0; k <= 9; k++) { // 9 = 3 * (number of nodes - 1)
        if (boundaryConditions[k] == 0) {
            forces[k] = Ancestor::initialVolume_ *
            (stress.Get(0, 0) * dNdX_[k] + stress.Get(1, 0) * dNdX_[k+1] + stress.Get(2, 0) * dNdX_[k+2]);
        } else {
            forces[k] = 0;
        }
        
        k++;
        if (boundaryConditions[k] == 0) {
            forces[k] = Ancestor::initialVolume_ *
            (stress.Get(0, 1) * dNdX_[k-1] + stress.Get(1, 1) * dNdX_[k] + stress.Get(2, 1) * dNdX_[k+1]);
        } else {
            forces[k] = 0;
        }
        
        k++;
        if (boundaryConditions[k] == 0) {
            forces[k] = Ancestor::initialVolume_ *
            (stress.Get(0, 2) * dNdX_[k-2] + stress.Get(1, 2) * dNdX_[k-1] + stress.Get(2, 2) * dNdX_[k]);
        } else {
            forces[k] = 0;
        }
    }
    
    return rc;
} // CBElementSolidT4::CalcNodalForcesHelperFunction

CBStatus CBElementSolidT4::CalcNodalForces() {
    CBStatus rc;
    TFloat       nodesCoords[12];
    bool         boundaryConditions[12];
    TInt         nodesCoordsIndices[12];
    
    Matrix3<TFloat> deformationTensor;
    
    GetNodesCoordsIndices(nodesCoordsIndices);
    Base::adapter_->GetNodesCoords(12, nodesCoordsIndices, nodesCoords);
    Ancestor::adapter_->GetNodesComponentsBoundaryConditions(12, nodesCoordsIndices, boundaryConditions);
    
    TFloat forces[12];
    
    rc = CalcNodalForcesHelperFunction(nodesCoords, boundaryConditions, forces);
    
    Base::adapter_->AddNodalForcesComponents(12, nodesCoordsIndices, forces);
    return rc;
} // CBElementSolidT4::CalcNodalForces

CBStatus CBElementSolidT4::GetDeformationTensor(Matrix3<TFloat> &f) {
    TFloat nodesCoords[12];
    bool   boundaryConditions[12];
    TInt   nodesCoordsIndices[12];
    
    GetNodesCoordsIndices(nodesCoordsIndices);
    
    Base::adapter_->GetNodesCoords(12, nodesCoordsIndices, nodesCoords);
    
    Ancestor::adapter_->GetNodesComponentsBoundaryConditions(12, nodesCoordsIndices, boundaryConditions);
    
    CalcDeformationTensorWithLocalBasis(nodesCoords, f);
    return CBStatus::SUCCESS;
}

TFloat *CBElementSolidT4::GetShapeFunctionsDerivatives() {
    return dNdX_.data();
}

TFloat CBElementSolidT4::GetDeformationEnergy() {
    TFloat e;
    Matrix3<TFloat> f;
    
    GetDeformationTensor(f);
    Base::material_->GetConstitutiveModel()->CalcEnergy(f, e);
    return initialVolume_ * e;
}

Matrix3<TFloat> CBElementSolidT4::GetPK2Stress() {
    Matrix3<TFloat> deformationTensor;
    
    GetDeformationTensor(deformationTensor);
    
    Matrix3<TFloat> pk2Stress = {0, 0, 0, 0, 0, 0, 0, 0, 0};
    CBStatus rc = Base::material_->GetConstitutiveModel()->CalcPK2Stress(deformationTensor, pk2Stress);
    
    return pk2Stress;
}

CBStatus CBElementSolidT4::GetCauchyStress(Matrix3<TFloat> &cauchyStress) {
    Matrix3<TFloat> deformationTensor;
    
    GetDeformationTensor(deformationTensor);
    
    Matrix3<TFloat> pk2Stress = {0, 0, 0, 0, 0, 0, 0, 0, 0};
    CBStatus rc = Base::material_->GetConstitutiveModel()->CalcPK2Stress(deformationTensor, pk2Stress);
    if (rc != CBStatus::SUCCESS)
        return rc;
    
    Matrix3<TFloat> a;
    
    // Base::adapter_->GetActiveStressTensor(localIndex_, a);
    TFloat time = adapter_->GetSolver()->GetTiming().GetCurrentTime();
    a = GetTensionModel()->CalcActiveStress(deformationTensor, time);
    
    pk2Stress += a;
    
    cauchyStress = 1.0 / deformationTensor.Det() *  deformationTensor * pk2Stress * deformationTensor.GetTranspose();
    
    return rc;
}

CBStatus CBElementSolidT4::CalcNodalForcesJacobian() {
    CBStatus rc;
    TFloat   nodesCoords[12];
    TFloat   nodeCoord;
    bool     boundaryConditions[12];
    TInt     nodesCoordsIndices[12];
    TInt     indices[12];
    TFloat   f1[12];
    TFloat   f2[12];
    TFloat   forcesJacobian[12*12]; // may be initialized to zero, but not necessary, as entries for fixated nodes are not copied (negative indices)
    
    GetNodesCoordsIndices(nodesCoordsIndices);
    Base::adapter_->GetNodesComponentsBoundaryConditions(12, nodesCoordsIndices, boundaryConditions);
    memcpy(indices, nodesCoordsIndices, 12*sizeof(TInt));
    
    Base::adapter_->GetNodesCoords(12, nodesCoordsIndices, nodesCoords);
    
    TFloat epsilon = Base::adapter_->GetFiniteDifferencesEpsilon();
    TFloat epsilon2 = 2.0*epsilon;
    
    // Calculate forces jacobian
    for (int i = 0; i < 12; i++) {
        if (boundaryConditions[i]) {
            indices[i] = -1; // negative indices will be ignored by MatSetValues
            continue;
        }
        
        /*
         * Finite Differences :
         * j = 0 -> node i + basisInverse * (epsilon,0,0),
         * j = 1 -> node i + basisInverse * (0,epsilon,0),
         * j = 2 -> node i + basisInverse * (0,0,epsilon)
         */
        
        nodeCoord = nodesCoords[i];
        
        nodesCoords[i] = nodeCoord + epsilon;
        rc = CalcNodalForcesHelperFunction(nodesCoords, boundaryConditions, f1);
        if (rc != CBStatus::SUCCESS)
            return rc;
        
        nodesCoords[i] = nodeCoord - epsilon;
        rc = CalcNodalForcesHelperFunction(nodesCoords, boundaryConditions, f2);
        
        nodesCoords[i] = nodeCoord;
        
        for (int k = 0; k <= 9; k++) { // 9 = 3 * (number of nodes - 1)
            forcesJacobian[12*k + i] = (f1[k] - f2[k]) / epsilon2;
            k++;
            forcesJacobian[12*k + i] = (f1[k] - f2[k]) / epsilon2;
            k++;
            forcesJacobian[12*k + i] = (f1[k] - f2[k]) / epsilon2;
        }
    }
    
    Base::adapter_->AddNodalForcesJacobianEntries(12, indices, 12, indices, forcesJacobian);
    return CBStatus::SUCCESS;
} // CBElementSolidT4::CalcNodalForcesJacobian

CBStatus CBElementSolidT4::CalcNodalForcesAndJacobian() {
    CBStatus rc;
    
    rc = CalcNodalForces();
    
    if (rc != CBStatus::SUCCESS)
        return rc;
    
    rc = CalcNodalForcesJacobian();
    return rc;
}

CBStatus CBElementSolidT4::CalcConsistentMassMatrix() {
    TFloat nodesCoords[12];
    bool   boundaryConditions[12];
    TInt   nodesCoordsIndices[12];
    
    GetNodesCoordsIndices(nodesCoordsIndices);
    
    Base::adapter_->GetNodesCoords(12, nodesCoordsIndices, nodesCoords);
    Base::adapter_->GetNodesComponentsBoundaryConditions(12, nodesCoordsIndices, boundaryConditions);
    
    /*  Consistent Mass Matrix for a 4-Node Tetrahedron
     *
     M =    pV/20 *
     | 2 0 0 | 1 0 0 | 1 0 0 | 1 0 0 |
     | 0 2 0 | 0 1 0 | 0 1 0 | 0 1 0 |
     | 0 0 2 | 0 0 1 | 0 0 1 | 0 0 1 |
     --------| ------| ------| -------
     | 1 0 0 | 2 0 0 | 1 0 0 | 1 0 0 |
     | 0 1 0 | 0 2 0 | 0 1 0 | 0 1 0 |
     | 0 0 1 | 0 0 2 | 0 0 1 | 0 0 1 |
     --------| ------| ------| -------
     | 1 0 0 | 1 0 0 | 2 0 0 | 1 0 0 |
     | 0 1 0 | 0 1 0 | 0 2 0 | 0 1 0 |
     | 0 0 1 | 0 0 1 | 0 0 2 | 0 0 1 |
     --------| ------| ------| -------
     | 1 0 0 | 1 0 0 | 1 0 0 | 2 0 0 |
     | 0 1 0 | 0 1 0 | 0 1 0 | 0 2 0 |
     | 0 0 1 | 0 0 1 | 0 0 1 | 0 0 2 |
     */
    
    TFloat c = (Base::material_->GetMassDensity() * GetVolume()) / 20;
    
    TFloat massMatrixEntries[144];
    
    TFloat r1[16] = { 2, 1, 1, 1,
        1, 2, 1, 1,
        1, 1, 2, 1,
        1, 1, 1, 2 };
    
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            for (int k = 0; k < 3; k++) {
                for (int l = 0; l < 3; l++) {
                    if ((k == l) && (boundaryConditions[3*i+k] == 0)) {
                        massMatrixEntries[12*(3*i+k)+(3*j+l)] = r1[4*i+j] * c;
                    } else {
                        massMatrixEntries[12*(3*i+k)+(3*j+l)] = 0;
                    }
                }
            }
        }
    }
    
    Base::adapter_->SetMassMatrixEntries(12, nodesCoordsIndices, 12, nodesCoordsIndices, massMatrixEntries);
    return CBStatus::SUCCESS;
} // CBElementSolidT4::CalcConsistentMassMatrix

CBStatus CBElementSolidT4::CalcLumpedMassMatrix() {
    TInt nodesCoordsIndices[12];
    
    GetNodesCoordsIndices(nodesCoordsIndices);
    
    /*  Lumped Mass Matrix for a 4-Node Tetrahedron
     *
     M =    pV/20 *
     | 5 0 0 0 0 0 0 0 0 0 0 0 |
     | 0 5 0 0 0 0 0 0 0 0 0 0 |
     | 0 0 5 0 0 0 0 0 0 0 0 0 |
     | 0 0 0 5 0 0 0 0 0 0 0 0 |
     | 0 0 0 0 5 0 0 0 0 0 0 0 |
     | 0 0 0 0 0 5 0 0 0 0 0 0 |
     | 0 0 0 0 0 0 5 0 0 0 0 0 |
     | 0 0 0 0 0 0 0 5 0 0 0 0 |
     | 0 0 0 0 0 0 0 0 5 0 0 0 |
     | 0 0 0 0 0 0 0 0 0 5 0 0 |
     | 0 0 0 0 0 0 0 0 0 0 5 0 |
     | 0 0 0 0 0 0 0 0 0 0 0 5 |
     */
    
    TFloat c = (Base::material_->GetMassDensity() * GetVolume()) / 4;
    TFloat massMatrixEntries[144];
    
    for (int i = 0; i < 12; i++)
        for (int j = 0; j < 12; j++)
            if (i == j)
                massMatrixEntries[12*i+j] = c;
            else
                massMatrixEntries[12*i+j] = 0;
    
    Base::adapter_->SetMassMatrixEntries(12, nodesCoordsIndices, 12, nodesCoordsIndices, massMatrixEntries);
    return CBStatus::SUCCESS;
}

bool CBElementSolidT4::IsElementInverted() {
    TFloat          nodesCoords[12];
    TInt            nodesCoordsIndices[12];
    Matrix3<TFloat> deformationTensor;
    
    GetNodesCoordsIndices(nodesCoordsIndices);
    Base::adapter_->GetNodesCoords(12, nodesCoordsIndices, nodesCoords);
    CalcDeformationTensorWithLocalBasis(nodesCoords, deformationTensor);
    
    if (deformationTensor.Det() < 0)
        return true;
    else
        return false;
}

TFloat CBElementSolidT4::GetVolume() {
    return detJ_ / 6.;
}

void CBElementSolidT4::CheckNodeSorting() {
    if (Base::adapter_ != 0) {
        TFloat nodesCoords[12];
        TInt   nodesCoordsIndices[12];
        
        GetNodesCoordsIndices(nodesCoordsIndices);
        
        Base::adapter_->GetNodesCoords(12, nodesCoordsIndices, nodesCoords);
        
        TFloat z41 = (nodesCoords[11] - nodesCoords[2]);
        TFloat z31 = (nodesCoords[8] - nodesCoords[2]);
        TFloat z21 = (nodesCoords[5] - nodesCoords[2]);
        TFloat y41 = (nodesCoords[10] - nodesCoords[1]);
        TFloat y31 = (nodesCoords[7] - nodesCoords[1]);
        TFloat y21 = (nodesCoords[4] - nodesCoords[1]);
        TFloat x41 = (nodesCoords[9] - nodesCoords[0]);
        TFloat x31 = (nodesCoords[6] - nodesCoords[0]);
        TFloat x21 = (nodesCoords[3] - nodesCoords[0]);
        
        TFloat v = x21 * (y31 * z41 - y41 * z31) + y21 * (x41 * z31 - x31 * z41) + z21 * (x31 * y41 - x41 * y31);
        
        if (v < 0) {
            TInt a = nodesIndices_[1];
            nodesIndices_[1] = nodesIndices_[2];
            nodesIndices_[2] = a;
        }
    } else {
        throw std::runtime_error("Adapter has to be set before running void CBElementSolidT4::CheckNodeSorting()");
    }
} // CBElementSolidT4::CheckNodeSorting

CBStatus CBElementSolidT4::CalculateLaplacian() {
    TInt nodesCoordsIndices[12];
    
    GetNodesCoordsIndices(nodesCoordsIndices);
    TFloat          nodesCoords[12];
    Base::adapter_->GetNodesCoords(12, nodesCoordsIndices, nodesCoords);
    TFloat laplacian[16];
    
    UpdateShapeFunctions();
    CalcShapeFunctionsDerivatives();
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            // Calculate nabla Ni * nabla Nj
            
            laplacian[4*i+j] = initialVolume_ *
            (dNdX_[3*i]*dNdX_[3*j]  + dNdX_[3*i+1]*dNdX_[3*j+1] + dNdX_[3*i+2]*dNdX_[3*j+2]);
    
    Base::adapter_->AddLaplacianEntriesGlobal(4, nodesIndices_.data(), 4, nodesIndices_.data(), laplacian);
    return CBStatus::SUCCESS;
}

CBStatus CBElementSolidT4::GetDeformationTensorAtQuadraturePoints(Matrix3<TFloat> *f) {
    for (int i = 0; i < 4; i++) {
        CBElementSolidT4::GetDeformationTensor(f[i]);
    }
    
    return CBStatus::SUCCESS;
}
