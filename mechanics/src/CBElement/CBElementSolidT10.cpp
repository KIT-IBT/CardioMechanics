/*
 *
 *  CBElementSolidT10.cpp
 *  CardioMechanics
 *
 *  Created by Thomas Fritz on 13.09.11.
 *  Copyright 2011 Institute for Biomedical Engineering, Karlsruhe Institute of Technology (KIT). All rights reserved.
 *
 */

#include "CBElementSolidT10.h"
#include "CBElementAdapter.h"
#include "Matrix4.h"
#include "CBData.h"
#include <algorithm>
#include "DCCtrl.h"

#include "CBSolver.h"
#include "CBData.h"

CBElementSolidT10::CBElementSolidT10(CBElementSolidT10 &other) : CBElementSolid(other) {
    nodesIndices_ = other.nodesIndices_;
    dNdXW_ = other.dNdXW_;
    dNdXt4_ = other.dNdXt4_;
    detJ_ = other.detJ_;
    basisAtQuadraturePoint_ = other.basisAtQuadraturePoint_;
}

CBElement *CBElementSolidT10::New() {
    return new CBElementSolidT10;
}

CBElement *CBElementSolidT10::Clone() {
    return new CBElementSolidT10(*this);
}

void CBElementSolidT10::UpdateShapeFunctions() {
    CalcShapeFunctionDerivativesAtQuadraturePoints();
    CalcShapeFunctionDerivativesAtCentroid();
    CalcT4ShapeFunctionsDerivatives();
}

void CBElementSolidT10::SetNodeIndex(unsigned int i, TInt j) {
    if (i > 10)
        throw std::runtime_error("CBElementSolidT10<T>::SetNode(unsigned int i,TInt j) -> i out of range");
    else
        nodesIndices_[i] = j;
}

TInt CBElementSolidT10::GetNodeIndex(unsigned int i) {
    if (i > 10)
        throw std::runtime_error("CBElementSolidT10<T>::SetNode(unsigned int i,TInt j) -> i out of range");
    else
        return nodesIndices_[i];
}

void CBElementSolidT10::SetBasisAtQuadraturePoint(int i, const Matrix3<TFloat> &basis) {
    if ((i >= 0) && (i < 5)) {
        basisAtQuadraturePoint_[i] = basis;
    } else {
        throw std::runtime_error(
                                 "void CBElementSolidT10::SetBasisAtQuadraturePoint(int i, const Matrix3<TFloat>& basis): This element has only 4 quadrature points, i is out of range");
    }
}

Matrix3<TFloat> *CBElementSolidT10::GetBasisAtQuadraturePoint(int i) {
    if ((i >= 0) && (i < 5)) {
        return &basisAtQuadraturePoint_[i];
    } else {
        throw std::runtime_error(
                                 "Matrix3<TFloat>* CBElementSolidT10::GetBasisAtQuadraturePoint(int i): This element has only 4 quadrature points, i is out of range");
    }
}

void CBElementSolidT10::CalcShapeFunctionDerivativesAtQuadraturePoints() {
    // RM: A bit irritating as Bases are stored with the centroid at [0] in bases - but correct!
    
    CalcShapeFunctionDerivatives((5+3*sqrt(5))/20, (5-sqrt(5))/20, (5-sqrt(5))/20, (5-sqrt(5))/20, &(dNdXW_[0]));
    CalcShapeFunctionDerivatives((5-sqrt(5))/20, (5+3*sqrt(5))/20, (5-sqrt(5))/20, (5-sqrt(5))/20, &(dNdXW_[30]));
    CalcShapeFunctionDerivatives((5-sqrt(5))/20, (5-sqrt(5))/20, (5+3*sqrt(5))/20, (5-sqrt(5))/20, &(dNdXW_[60]));
    CalcShapeFunctionDerivatives((5-sqrt(5))/20, (5-sqrt(5))/20, (5-sqrt(5))/20, (5+3*sqrt(5))/20, &(dNdXW_[90]));
    CalcShapeFunctionDerivatives(0.25, 0.25,  0.25, 0.25, &(dNdXW_[120]));
    CalcT4ShapeFunctionsDerivatives();
    Ancestor::initialVolume_ = GetVolume();
}

void CBElementSolidT10::CalcShapeFunctionDerivativesAtCentroid()
{}

void CBElementSolidT10::CalcShapeFunctionDerivatives(TFloat l1, TFloat l2, TFloat l3, TFloat l4, TFloat *dNdX,
                                                     bool useReferenceNodes) {
    
    TFloat nodesCoords[30];
    TInt   nodesCoordsIndices[30];
    
    for (unsigned int i = 0; i < 10; i++) {
        nodesCoordsIndices[3*i]   = 3*nodesIndices_[i];
        nodesCoordsIndices[3*i+1] = 3*nodesIndices_[i]+1;
        nodesCoordsIndices[3*i+2] = 3*nodesIndices_[i]+2;
    }
    
    if (useReferenceNodes) {
        Base::adapter_->GetRefNodesCoords(30, nodesCoordsIndices, nodesCoords);
    } else {
        Base::adapter_->GetNodesCoords(30, nodesCoordsIndices, nodesCoords);
    }
    TFloat x1 = nodesCoords[0];
    TFloat y1 = nodesCoords[1];
    TFloat z1 = nodesCoords[2];
    
    TFloat x2 = nodesCoords[3];
    TFloat y2 = nodesCoords[4];
    TFloat z2 = nodesCoords[5];
    
    TFloat x3 = nodesCoords[6];
    TFloat y3 = nodesCoords[7];
    TFloat z3 = nodesCoords[8];
    
    TFloat x4 = nodesCoords[9];
    TFloat y4 = nodesCoords[10];
    TFloat z4 = nodesCoords[11];
    
    TFloat x5 = nodesCoords[12];
    TFloat y5 = nodesCoords[13];
    TFloat z5 = nodesCoords[14];
    
    TFloat x6 = nodesCoords[15];
    TFloat y6 = nodesCoords[16];
    TFloat z6 = nodesCoords[17];
    
    TFloat x7 = nodesCoords[18];
    TFloat y7 = nodesCoords[19];
    TFloat z7 = nodesCoords[20];
    
    TFloat x8 = nodesCoords[21];
    TFloat y8 = nodesCoords[22];
    TFloat z8 = nodesCoords[23];
    
    TFloat x9 = nodesCoords[24];
    TFloat y9 = nodesCoords[25];
    TFloat z9 = nodesCoords[26];
    
    TFloat x10 = nodesCoords[27];
    TFloat y10 = nodesCoords[28];
    TFloat z10 = nodesCoords[29];
    
    TFloat Jx1 = 4.0*(x1*(l1-0.25)+x5*l2+x7*l3+x8*l4);
    TFloat Jy1 = 4.0*(y1*(l1-0.25)+y5*l2+y7*l3+y8*l4);
    TFloat Jz1 = 4.0*(z1*(l1-0.25)+z5*l2+z7*l3+z8*l4);
    
    TFloat Jx2 = 4.0*(x5*l1+x2*(l2-0.25)+x6*l3+x9*l4);
    TFloat Jy2 = 4.0*(y5*l1+y2*(l2-0.25)+y6*l3+y9*l4);
    TFloat Jz2 = 4.0*(z5*l1+z2*(l2-0.25)+z6*l3+z9*l4);
    
    TFloat Jx3 = 4.0*(x7*l1+x6*l2+x3*(l3-0.25)+x10*l4);
    TFloat Jy3 = 4.0*(y7*l1+y6*l2+y3*(l3-0.25)+y10*l4);
    TFloat Jz3 = 4.0*(z7*l1+z6*l2+z3*(l3-0.25)+z10*l4);
    
    TFloat Jx4 = 4.0*(x8*l1+x9*l2+x10*l3+x4*(l4-0.25));
    TFloat Jy4 = 4.0*(y8*l1+y9*l2+y10*l3+y4*(l4-0.25));
    TFloat Jz4 = 4.0*(z8*l1+z9*l2+z10*l3+z4*(l4-0.25));
    
    TFloat Jx12 = Jx1-Jx2;
    TFloat Jx13 = Jx1-Jx3;
    TFloat Jx14 = Jx1-Jx4;
    TFloat Jx23 = Jx2-Jx3;
    TFloat Jx24 = Jx2-Jx4;
    TFloat Jx34 = Jx3-Jx4;
    
    TFloat Jy12 = Jy1-Jy2;
    TFloat Jy13 = Jy1-Jy3;
    TFloat Jy14 = Jy1-Jy4;
    TFloat Jy23 = Jy2-Jy3;
    TFloat Jy24 = Jy2-Jy4;
    TFloat Jy34 = Jy3-Jy4;
    
    TFloat Jz12 = Jz1-Jz2;
    TFloat Jz13 = Jz1-Jz3;
    TFloat Jz14 = Jz1-Jz4;
    TFloat Jz23 = Jz2-Jz3;
    TFloat Jz24 = Jz2-Jz4;
    TFloat Jz34 = Jz3-Jz4;
    
    TFloat Jx21 = -Jx12;
    TFloat Jx31 = -Jx13;
    
    //    TFloat Jx41 = -Jx14;
    TFloat Jx32 = -Jx23;
    TFloat Jx42 = -Jx24;
    TFloat Jx43 = -Jx34;
    TFloat Jy21 = -Jy12;
    TFloat Jy31 = -Jy13;
    
    //  TFloat Jy41 = -Jy14;
    TFloat Jy32 = -Jy23;
    TFloat Jy42 = -Jy24;
    TFloat Jy43 = -Jy34;
    TFloat Jz21 = -Jz12;
    TFloat Jz31 = -Jz13;
    
    //  TFloat Jz41 = -Jz14;
    TFloat Jz32 = -Jz23;
    TFloat Jz42 = -Jz24;
    TFloat Jz43 = -Jz34;
    detJ_ = Jx21*(Jy23*Jz34-Jy34*Jz23)+Jx32*(Jy34*Jz12-Jy12*Jz34)+Jx43*(Jy12*Jz23-Jy23*Jz12);
    
    
    TFloat a1 = Jy42*Jz32-Jy32*Jz42;
    TFloat a2 = Jy31*Jz43-Jy34*Jz13;
    TFloat a3 = Jy24*Jz14-Jy14*Jz24;
    TFloat a4 = Jy13*Jz21-Jy12*Jz31;
    
    TFloat b1 = Jx32*Jz42-Jx42*Jz32;
    TFloat b2 = Jx43*Jz31-Jx13*Jz34;
    TFloat b3 = Jx14*Jz24-Jx24*Jz14;
    TFloat b4 = Jx21*Jz13-Jx31*Jz12;
    
    TFloat c1 = Jx42*Jy32-Jx32*Jy42;
    TFloat c2 = Jx31*Jy43-Jx34*Jy13;
    TFloat c3 = Jx24*Jy14-Jx14*Jy24;
    TFloat c4 = Jx13*Jy21-Jx12*Jy31;
    
    TFloat factor = 4.0/detJ_;
    
    // Nfx
    dNdX[0]  = factor * (l1-0.25)*a1;
    dNdX[3]  = factor * (l2-0.25)*a2;
    dNdX[6]  = factor * (l3-0.25)*a3;
    dNdX[9]  = factor * (l4-0.25)*a4;
    dNdX[12] = factor * (l1*a2+l2*a1);
    dNdX[15] = factor * (l2*a3+l3*a2);
    dNdX[18] = factor * (l3*a1+l1*a3);
    dNdX[21] = factor * (l1*a4+l4*a1);
    dNdX[24] = factor * (l2*a4+l4*a2);
    dNdX[27] = factor * (l3*a4+l4*a3);
    
    // Nfy
    dNdX[0+1]  = factor * (l1-0.25)*b1;
    dNdX[3+1]  = factor * (l2-0.25)*b2;
    dNdX[6+1]  = factor * (l3-0.25)*b3;
    dNdX[9+1]  = factor * (l4-0.25)*b4;
    dNdX[12+1] = factor * (l1*b2+l2*b1);
    dNdX[15+1] = factor * (l2*b3+l3*b2);
    dNdX[18+1] = factor * (l3*b1+l1*b3);
    dNdX[21+1] = factor * (l1*b4+l4*b1);
    dNdX[24+1] = factor * (l2*b4+l4*b2);
    dNdX[27+1] = factor * (l3*b4+l4*b3);
    
    // Nfz
    dNdX[0+2]  = factor * (l1-0.25)*c1;
    dNdX[3+2]  = factor * (l2-0.25)*c2;
    dNdX[6+2]  = factor * (l3-0.25)*c3;
    dNdX[9+2]  = factor * (l4-0.25)*c4;
    dNdX[12+2] = factor * (l1*c2+l2*c1);
    dNdX[15+2] = factor * (l2*c3+l3*c2);
    dNdX[18+2] = factor * (l3*c1+l1*c3);
    dNdX[21+2] = factor * (l1*c4+l4*c1);
    dNdX[24+2] = factor * (l2*c4+l4*c2);
    dNdX[27+2] = factor * (l3*c4+l4*c3);
} // CBElementSolidT10::CalcShapeFunctionDerivatives

void CBElementSolidT10::CalcDeformationTensorsAtQuadraturePointsWithLocalBasis(const TFloat *nodesCoords, Matrix3<TFloat> *deformationTensors) {
    for (unsigned int n = 0; n < 4; n++) {
        TFloat *dNdX = &dNdXW_[30*n];         // Derivatives of shape functions at quadrature point n 0-4 = QuadP ; n=5/120 = center
        TFloat *f    = deformationTensors[n].GetArray();
        
        for (unsigned int i = 0; i < 3; i++) { // x/y/z
            for (int j = 0; j < 3; j++) { // dX/dY/dZ
                TFloat e = 0;
                
                for (int k = 0; k < 10; k++) {
                    e += dNdX[3*k+j] * nodesCoords[3*k+i];
                }
                f[3*i+j] = e;
            }
        }
    }
    
    
    for (int i = 0; i < 4; i++) {
        deformationTensors[i] =  GetBasisAtQuadraturePoint(i+1)->GetTranspose() * deformationTensors[i] *
        GetBasisAtQuadraturePoint(i+1)->GetInverse().GetTranspose();
    }
} // CBElementSolidT10::CalcDeformationTensorsAtQuadraturePointsWithLocalBasis

void CBElementSolidT10::CalcDeformationTensorsAtQuadraturePointsWithGlobalBasis(const TFloat *nodesCoords, Matrix3<TFloat> *deformationTensors) {
    for (unsigned int n = 0; n < 4; n++) {
        TFloat *dNdX = &dNdXW_[30*n];         // Derivatives of shape functions at quadrature point n 0-4 = QuadP ; n=5/120 = center
        TFloat *f    = deformationTensors[n].GetArray();
        
        for (unsigned int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++) {
                TFloat e = 0;
                
                for (int k = 0; k < 10; k++)
                    e += dNdX[3*k+j] * nodesCoords[3*k+i];
                
                f[3*i+j] = e;
            }
    }
}

void CBElementSolidT10::CalcDeformationTensorsAtCentroidWithLocalBasis(const TFloat *nodesCoords, Matrix3<TFloat> &deformationTensor) {
    TFloat *dNdX = &dNdXW_[120];
    TFloat *f    = deformationTensor.GetArray();
    
    for (unsigned int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) {
            TFloat e = 0;
            
            for (int k = 0; k < 10; k++) {
                e += dNdX[3*k+j] * nodesCoords[3*k+i];
            }
            f[3*i+j] = e;
        }
    
    deformationTensor = GetBasisAtQuadraturePoint(0)->GetTranspose() * deformationTensor * *GetBasisAtQuadraturePoint(0);
}

void CBElementSolidT10::CalcDeformationTensorsAtCentroidWithLocalBasisWithT4ShapeFunctions(const TFloat *nodesCoords, Matrix3<TFloat> &deformationTensor)
{
    TFloat *f = deformationTensor.GetArray();
    
    for (unsigned int i = 0; i < 3; i++) {
        f[i]   = dNdXt4_[i] * nodesCoords[0] + dNdXt4_[i+3] * nodesCoords[3] + dNdXt4_[i+6] * nodesCoords[6] +
        dNdXt4_[i+9] * nodesCoords[9];
        f[3+i] = dNdXt4_[i] * nodesCoords[1] + dNdXt4_[i+3] * nodesCoords[4] + dNdXt4_[i+6] * nodesCoords[7] +
        dNdXt4_[i+9] * nodesCoords[10];
        f[6+i] = dNdXt4_[i] * nodesCoords[2] + dNdXt4_[i+3] * nodesCoords[5] + dNdXt4_[i+6] * nodesCoords[8] +
        dNdXt4_[i+9] * nodesCoords[11];
    }
    
    deformationTensor = GetBasisAtQuadraturePoint(0)->GetTranspose() * deformationTensor * *GetBasisAtQuadraturePoint(0);
}

CBStatus CBElementSolidT10::CalcNodalForcesAndJacobian() {
    CBStatus rc;
    
    rc = CalcNodalForces();
    
    if (rc != CBStatus::SUCCESS)
        return rc;
    
    rc = CalcNodalForcesJacobian();
    return rc;
}

CBStatus CBElementSolidT10::CalcConsistentMassMatrix() {
    TFloat nodesCoords[30];
    bool   boundaryConditions[30];
    TInt   nodesCoordsIndices[30];
    
    for (unsigned int i = 0; i < 10; i++) {
        nodesCoordsIndices[3*i]   = 3*nodesIndices_[i];
        nodesCoordsIndices[3*i+1] = 3*nodesIndices_[i]+1;
        nodesCoordsIndices[3*i+2] = 3*nodesIndices_[i]+2;
    }
    
    Base::adapter_->GetNodesCoords(30, nodesCoordsIndices, nodesCoords);
    Base::adapter_->GetNodesComponentsBoundaryConditions(30, nodesCoordsIndices, boundaryConditions);
    
    TFloat c = (Base::material_->GetMassDensity() * detJ_) / 2520;
    
    TFloat massMatrixEntries[900];
    
    // adapted from http://arxiv.org/pdf/1411.1341.pdf
    TFloat r1[100] = {6,  1,  1,  1, -4, -6, -4, -4, -6, -6,
                      1,  6,  1,  1, -4, -4, -6, -6, -4, -6,
                      1,  1,  6,  1, -6, -4, -4, -6, -6, -4,
                      1,  1,  1,  6, -6, -6, -6, -4, -4, -4,
                      -4, -4, -6, -6, 32, 16, 16, 16, 16,  8,
                      -6, -4, -4, -6, 16, 32, 16,  8, 16, 16,
                      -4, -6, -4, -6, 16, 16, 32, 16,  8, 16,
                      -4, -6, -6, -4, 16,  8, 16, 32, 16, 16,
                      -6, -4, -6, -4, 16, 16,  8, 16, 32, 16,
                      -6, -6, -4, -4,  8, 16, 16, 16, 16, 32};
    
    for (int i = 0; i < 10; i++)
        for (int j = 0; j < 10; j++)
            for (int k = 0; k < 3; k++)
                for (int l = 0; l < 3; l++)
                    if ((k == l) && (boundaryConditions[3*i+k] == 0))
                        massMatrixEntries[30*(3*i+k)+(3*j+l)] = r1[10*i+j] * c;
                    else
                        massMatrixEntries[30*(3*i+k)+(3*j+l)] = 0;
    
    Base::adapter_->SetMassMatrixEntries(30, nodesCoordsIndices, 30, nodesCoordsIndices, massMatrixEntries);
    
    return CBStatus::SUCCESS;
} // CBElementSolidT10::CalcConsistentMassMatrix

CBStatus CBElementSolidT10::CalcLumpedMassMatrix() {
    TFloat nodesCoords[30];
    bool   boundaryConditions[30];
    TInt   nodesCoordsIndices[30];
    
    for (unsigned int i = 0; i < 10; i++) {
        nodesCoordsIndices[3*i]   = 3*nodesIndices_[i];
        nodesCoordsIndices[3*i+1] = 3*nodesIndices_[i]+1;
        nodesCoordsIndices[3*i+2] = 3*nodesIndices_[i]+2;
    }
    
    Base::adapter_->GetNodesCoords(30, nodesCoordsIndices, nodesCoords);
    Base::adapter_->GetNodesComponentsBoundaryConditions(30, nodesCoordsIndices, boundaryConditions);
    
    TFloat c = (Base::material_->GetMassDensity() * GetVolume()) / 10;
    TFloat massMatrixEntries[900];
    
    for (int i = 0; i < 30; i++)
        for (int j = 0; j < 30; j++)
            if ((i == j) && (boundaryConditions[i] == 0))
                massMatrixEntries[30*i+j] = c;
            else
                massMatrixEntries[30*i+j] = 0;
    
    Base::adapter_->SetMassMatrixEntries(30, nodesCoordsIndices, 30, nodesCoordsIndices, massMatrixEntries);
    return CBStatus::SUCCESS;
} // CBElementSolidT10::CalcLumpedMassMatrix

void CBElementSolidT10::GetNodesCoordsIndices(TInt *nodesCoordsIndices) {
    for (unsigned int i = 0; i < 10; i++) {
        nodesCoordsIndices[3*i]   = 3*nodesIndices_[i];
        nodesCoordsIndices[3*i+1] = 3*nodesIndices_[i]+1;
        nodesCoordsIndices[3*i+2] = 3*nodesIndices_[i]+2;
    }
}

CBStatus CBElementSolidT10::SetNodalForcesToZeroIfElementIsDefect() {
    if (isDefect_) {
        TInt     nodesCoordsIndices[30];
        GetNodesCoordsIndices(nodesCoordsIndices);
        TFloat forces[30] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        
        std::cout << "Deactivating element: " << GetIndex() << "\n";
        
        Base::adapter_->InsertNodalForcesComponents(30, nodesCoordsIndices, forces);
    }
    
    return CBStatus::SUCCESS;
}

CBStatus CBElementSolidT10::SetNodalForcesJacobianToZeroIfElementIsDefect() {
    if (isDefect_) {
        TInt     nodesCoordsIndices[30];
        GetNodesCoordsIndices(nodesCoordsIndices);
        TFloat jacobian[900];
        
        for (int i = 0; i < 30; i++)
            for (int j = 0; j < 30; j++)
                if (i == j)
                    jacobian[30*i+j] = 1;
                else
                    jacobian[30*i+j] = 0;
        std::cout << "Deactivating element: " << GetIndex() << "\n";
        Base::adapter_->InsertNodalForcesJacobianEntries(30, nodesCoordsIndices, 30, nodesCoordsIndices, jacobian);
    }
    
    return CBStatus::SUCCESS;
}

CBStatus CBElementSolidT10::CalcNodalForces() {
    CBStatus rc;
    TFloat   nodesCoords[30];
    bool     boundaryConditions[30];
    TInt     nodesCoordsIndices[30];
    
    GetNodesCoordsIndices(nodesCoordsIndices);
    
    Base::adapter_->GetNodesCoords(30, nodesCoordsIndices, nodesCoords);
    Base::adapter_->GetNodesComponentsBoundaryConditions(30, nodesCoordsIndices, boundaryConditions);
    
    TFloat forces[30];
    
    rc = CalcNodalForcesHelperFunction(nodesCoords, boundaryConditions, forces);
    
    if (rc != CBStatus::SUCCESS)
        return rc;
    
    Base::adapter_->AddNodalForcesComponents(30, nodesCoordsIndices, forces);
    
    return CBStatus::SUCCESS;
} // CBElementSolidT10::CalcNodalForces

CBStatus CBElementSolidT10::CalcNodalForcesWithoutActiveStress() {
    CBStatus rc;
    TFloat   nodesCoords[30];
    bool     boundaryConditions[30];
    TInt     nodesCoordsIndices[30];
    
    GetNodesCoordsIndices(nodesCoordsIndices);
    
    Base::adapter_->GetNodesCoords(30, nodesCoordsIndices, nodesCoords);
    Base::adapter_->GetNodesComponentsBoundaryConditions(30, nodesCoordsIndices, boundaryConditions);
    
    TFloat forces[30];
    
    rc = CalcNodalForcesHelperFunction(nodesCoords, boundaryConditions, forces);
    
    if (rc != CBStatus::SUCCESS)
        return rc;
    
    Base::adapter_->AddNodalForcesComponents(30, nodesCoordsIndices, forces);
    
    return CBStatus::SUCCESS;
}

CBStatus CBElementSolidT10::CalcNodalForcesJacobian() {
    CBStatus rc;
    TFloat   nodesCoords[30];
    TFloat   nodeCoord;
    bool     boundaryConditions[30];
    TInt     nodesCoordsIndices[30];
    TInt     indices[30];
    TFloat   f1[30];
    TFloat   f2[30];
    TFloat   forcesJacobian[30*30] = {0.0}; // may be initialized to zero, but not necessary, as entries for fixated nodes are not copied (negative indices)
    
    GetNodesCoordsIndices(nodesCoordsIndices);
    Base::adapter_->GetNodesComponentsBoundaryConditions(30, nodesCoordsIndices, boundaryConditions);
    memcpy(indices, nodesCoordsIndices, 30*sizeof(TInt));
    
    Base::adapter_->GetNodesCoords(30, nodesCoordsIndices, nodesCoords);
    
    TFloat epsilon = Base::adapter_->GetFiniteDifferencesEpsilon();
    TFloat epsilon2 = 2.0*epsilon;
    
    //    std::cout << epsilon << std::endl;
    // Calculate forces jacobian
    for (int i = 0; i < 30; i++) {
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
        
        for (int k = 0; k < 10; k++) {
            forcesJacobian[30*i + 3*k + 0] = (f1[3*k + 0] - f2[3*k + 0]) / epsilon2;
            forcesJacobian[30*i + 3*k + 1] = (f1[3*k + 1] - f2[3*k + 1]) / epsilon2;
            forcesJacobian[30*i + 3*k + 2] = (f1[3*k + 2] - f2[3*k + 2]) / epsilon2;
        }
    }
    
    Base::adapter_->AddNodalForcesJacobianEntries(30, indices, 30, indices, forcesJacobian);
    return CBStatus::SUCCESS;
} // CBElementSolidT10::CalcNodalForcesJacobian

TFloat CBElementSolidT10::GetVolume() {
    return detJ_/6;
}

void CBElementSolidT10::CheckNodeSorting() {
    if (Base::adapter_ != 0) {
        TFloat nodesCoords[30];
        TInt   nodesCoordsIndices[30];
        
        for (unsigned int i = 0; i < 10; i++) {
            nodesCoordsIndices[3*i]   = 3*nodesIndices_[i];
            nodesCoordsIndices[3*i+1] = 3*nodesIndices_[i]+1;
            nodesCoordsIndices[3*i+2] = 3*nodesIndices_[i]+2;
        }
        
        Base::adapter_->GetNodesCoords(30, nodesCoordsIndices, nodesCoords);
        
        TFloat z41 = (nodesCoords[11] - nodesCoords[2]);
        TFloat z31 = (nodesCoords[8] - nodesCoords[2]);
        TFloat z21 = (nodesCoords[5] - nodesCoords[2]);
        TFloat y41 = (nodesCoords[10] - nodesCoords[1]);
        TFloat y31 = (nodesCoords[7] - nodesCoords[1]);
        TFloat y21 = (nodesCoords[4] - nodesCoords[1]);
        TFloat x41 = (nodesCoords[9] - nodesCoords[0]);
        TFloat x31 = (nodesCoords[6] - nodesCoords[0]);
        TFloat x21 = (nodesCoords[3] - nodesCoords[0]);
        
        TFloat v = (x21 * (y31 * z41 - y41 * z31) + y21 * (x41 * z31 - x31 * z41) + z21 * (x31 * y41 - x41 * y31));
        
        if (v < 0) {
            TInt a = nodesIndices_[1];
            nodesIndices_[1] = nodesIndices_[2];
            nodesIndices_[2] = a;
            
            a = nodesIndices_[6];
            nodesIndices_[6] = nodesIndices_[4];
            nodesIndices_[4] = a;
            
            a = nodesIndices_[8];
            nodesIndices_[8] = nodesIndices_[9];
            nodesIndices_[9] = a;
            
            Matrix3<TFloat> b = *GetBasisAtQuadraturePoint(2);
            SetBasisAtQuadraturePoint(2, *GetBasisAtQuadraturePoint(3));
            SetBasisAtQuadraturePoint(3, b);
        }
    } else {
        throw std::runtime_error("Adapter has to be set before running void CBElementSolidT10::CheckNodeSorting()");
    }
} // CBElementSolidT10::CheckNodeSorting

CBStatus CBElementSolidT10::CalcNodalForcesHelperFunction(const TFloat *nodesCoords, const bool *boundaryConditions, TFloat *forces) {
    CBStatus rc;
    Matrix3<TFloat> activeStress    = {0, 0, 0, 0, 0, 0, 0, 0, 0};
    Matrix3<TFloat> deformationTensors[4];
    Matrix3<TFloat> stress[4]       = { {0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0}   };
    TFloat time = adapter_->GetSolver()->GetTiming().GetCurrentTime();
    
    CalcDeformationTensorsAtQuadraturePointsWithLocalBasis(nodesCoords, deformationTensors);
    
    for (int QPi = 0; QPi < 4; QPi++) {
        rc = Base::material_->GetConstitutiveModel()->CalcPK2Stress(deformationTensors[QPi], stress[QPi]);
        
        // print corrupt element messages from each processor
        if (rc == CBStatus::CORRUPT_ELEMENT) {
            std::cout << "SolidT10::CalcNodalForcesHelperFunction(): Element with index " << index_ <<
            " is corrupt. Det = " << deformationTensors[QPi].Det() << std::endl;
        }
        
        if (rc != CBStatus::SUCCESS)
            return rc;
        
        // until here, stress holds the PK2 stress
        activeStress = Base::tensionModel_->CalcActiveStress(deformationTensors[QPi], time);
        stress[QPi] += activeStress;
        
        // convert PK2 stress into nominal stress with respect to the local coordinate system aligned with the fibres
        stress[QPi] = (deformationTensors[QPi] * stress[QPi]).GetTranspose();
        
        // transform  stress back into the global coordinate system
        stress[QPi] = GetBasisAtQuadraturePoint(QPi+1)->GetInverse().GetTranspose() * stress[QPi] * GetBasisAtQuadraturePoint(QPi+1)->GetTranspose(); // +1 since centroid is at i = 0
    }
    
    // Calculate force contribution from all nodes with 2nd order quadrature rule: W_i = 1/4, tetCoords = a,b,b,b
    // f_i = V * 1/4 * [ f_j(a,b,b,b) + f_j(b,a,b,b) + f_j(b,b,a,b) + f_j(b,b,b,a) ]
    TFloat t[3];
    double V = detJ_/24.0; // technically here V = det(J)/6 * 1/4
    for (int i = 0; i < 10; i++) { // for each node
        t[0] = 0;
        t[1] = 0;
        t[2] = 0;
        
        for (int j = 0; j < 4; j++) { // for each QP
            t[0] += dNdXW_[30*j + 3*i+0] *
            stress[j].Get(0, 0)  + dNdXW_[30*j + 3*i+1] * stress[j].Get(1, 0) + dNdXW_[30*j + 3*i+2] * stress[j].Get(2, 0);
            t[1] += dNdXW_[30*j + 3*i+0] *
            stress[j].Get(0, 1)  + dNdXW_[30*j + 3*i+1] * stress[j].Get(1, 1) + dNdXW_[30*j + 3*i+2] * stress[j].Get(2, 1);
            t[2] += dNdXW_[30*j + 3*i+0] *
            stress[j].Get(0, 2)  + dNdXW_[30*j + 3*i+1] * stress[j].Get(1, 2) + dNdXW_[30*j + 3*i+2] * stress[j].Get(2, 2);
        }
        
        // check if node i has boundary Conditions in x direction
        if (boundaryConditions[3*i] == 0)
            forces[3*i] = V * t[0];
        else
            forces[3*i] = 0;
        
        // check if node has BC in y dir
        if (boundaryConditions[3*i+1] == 0)
            forces[3*i+1] = V * t[1];
        else
            forces[3*i+1] = 0;
        
        // check if node has BC in z dir
        if (boundaryConditions[3*i+2] == 0)
            forces[3*i+2] = V * t[2];
        else
            forces[3*i+2] = 0;
    }
    
    return rc;
} // CBElementSolidT10::CalcNodalForcesHelperFunction

void CBElementSolidT10::CalcT4ShapeFunctionsDerivatives() {
    TFloat nodesCoords[30];
    TInt   nodesCoordsIndices[30];
    
    GetNodesCoordsIndices(nodesCoordsIndices);
    Base::adapter_->GetNodesCoords(30, nodesCoordsIndices, nodesCoords);
    
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
    
    TFloat v = x21 * (y31 * z41 - y41 * z31) + y21 * (x41 * z31 - x31 * z41) + z21 * (x31 * y41 - x41 * y31);
    
    dNdXt4_[0] = 1.0 /v*(nodesCoords[4] * z43 - nodesCoords[7] * z42 + nodesCoords[10] * z32);
    dNdXt4_[3] = 1.0 /v*(-nodesCoords[1] * z43 + nodesCoords[7] * z41 - nodesCoords[10] * z31);
    dNdXt4_[6] = 1.0 /v*(nodesCoords[1] * z42 - nodesCoords[4] * z41 + nodesCoords[10] * z21);
    dNdXt4_[9] = 1.0 /v*(-nodesCoords[1] * z32 + nodesCoords[4] * z31 - nodesCoords[7] * z21);
    
    dNdXt4_[1] = 1.0 /v*(-nodesCoords[3] * z43 + nodesCoords[6] * z42 - nodesCoords[9] * z32);
    dNdXt4_[4] = 1.0 /v*(nodesCoords[0] * z43 - nodesCoords[6] * z41 + nodesCoords[9] * z31);
    dNdXt4_[7] = 1.0 /v*(-nodesCoords[0] * z42 + nodesCoords[3] * z41 - nodesCoords[9] * z21);
    dNdXt4_[10] = 1.0 /v*(nodesCoords[0] * z32 - nodesCoords[3] * z31 + nodesCoords[6] * z21);
    
    dNdXt4_[2] = 1.0 /v*(nodesCoords[3] * y43 - nodesCoords[6] * y42 + nodesCoords[9] * y32);
    dNdXt4_[5] = 1.0 /v*(-nodesCoords[0] * y43 + nodesCoords[6] * y41 - nodesCoords[9] * y31);
    dNdXt4_[8] = 1.0 /v*(nodesCoords[0] * y42 - nodesCoords[3] * y41 + nodesCoords[9] * y21);
    dNdXt4_[11] = 1.0 /v*(-nodesCoords[0] * y32 + nodesCoords[3] * y31 - nodesCoords[6] * y21);
} // CBElementSolidT10::CalcT4ShapeFunctionsDerivatives

TFloat *CBElementSolidT10::GetShapeFunctionsDerivatives() {
    return dNdXW_.data();
}

TFloat *CBElementSolidT10::GetT4ShapeFunctionsDerivatives() {
    return dNdXt4_.data();
}

bool CBElementSolidT10::IsElementInverted() {
    TFloat nodesCoords[30];
    TInt   nodesCoordsIndices[30];
    
    Matrix3<TFloat> deformationTensors[4];
    
    GetNodesCoordsIndices(nodesCoordsIndices);
    Base::adapter_->GetNodesCoords(30, nodesCoordsIndices, nodesCoords);
    CalcDeformationTensorsAtQuadraturePointsWithLocalBasis(nodesCoords, deformationTensors);
    
    for (int i = 0; i < 5; i++)
        if (deformationTensors[i].Det() < 0)
            return true;
    
    return false;
}

CBStatus CBElementSolidT10::GetDeformationTensor(Matrix3<TFloat> &f) {
    TFloat nodesCoords[30];
    TInt   nodesCoordsIndices[30];
    
    for (unsigned int i = 0; i < 10; i++) {
        nodesCoordsIndices[3*i]   = 3*nodesIndices_[i];
        nodesCoordsIndices[3*i+1] = 3*nodesIndices_[i]+1;
        nodesCoordsIndices[3*i+2] = 3*nodesIndices_[i]+2;
    }
    
    Base::adapter_->GetNodesCoords(30, nodesCoordsIndices, nodesCoords);
    CBElementSolidT10::CalcDeformationTensorsAtCentroidWithLocalBasis(nodesCoords, f);
    return CBStatus::SUCCESS;
}

CBStatus CBElementSolidT10::GetDeformationTensorAtQuadraturePoints(Matrix3<TFloat> *f) {
    TFloat nodesCoords[30];
    TInt   nodesCoordsIndices[30];
    
    for (unsigned int i = 0; i < 10; i++) {
        nodesCoordsIndices[3*i]   = 3*nodesIndices_[i];
        nodesCoordsIndices[3*i+1] = 3*nodesIndices_[i]+1;
        nodesCoordsIndices[3*i+2] = 3*nodesIndices_[i]+2;
    }
    
    Base::adapter_->GetNodesCoords(30, nodesCoordsIndices, nodesCoords);
    CBElementSolidT10::CalcDeformationTensorsAtQuadraturePointsWithLocalBasis(nodesCoords, f);
    return CBStatus::SUCCESS;
}

TFloat CBElementSolidT10::GetDeformationEnergy() {
    TFloat nodesCoords[30];
    TInt   nodesCoordsIndices[30];
    
    for (unsigned int i = 0; i < 10; i++) {
        nodesCoordsIndices[3*i]   = 3*nodesIndices_[i];
        nodesCoordsIndices[3*i+1] = 3*nodesIndices_[i]+1;
        nodesCoordsIndices[3*i+2] = 3*nodesIndices_[i]+2;
    }
    
    Base::adapter_->GetNodesCoords(30, nodesCoordsIndices, nodesCoords);
    Matrix3<TFloat> deformationTensors[4];
    CalcDeformationTensorsAtQuadraturePointsWithLocalBasis(nodesCoords, deformationTensors);
    TFloat e1, e2, e3, e4;
    
    Base::material_->GetConstitutiveModel()->CalcEnergy(deformationTensors[0], e1);
    Base::material_->GetConstitutiveModel()->CalcEnergy(deformationTensors[1], e2);
    Base::material_->GetConstitutiveModel()->CalcEnergy(deformationTensors[2], e3);
    Base::material_->GetConstitutiveModel()->CalcEnergy(deformationTensors[3], e4);
    
    return 0.25 *initialVolume_ * (e1+e2+e3+e4);
}

Matrix3<TFloat> CBElementSolidT10::GetPK2Stress() {
    TFloat nodesCoords[30];
    TInt   nodesCoordsIndices[30];
    
    for (unsigned int i = 0; i < 10; i++) {
        nodesCoordsIndices[3*i]   = 3*nodesIndices_[i];
        nodesCoordsIndices[3*i+1] = 3*nodesIndices_[i]+1;
        nodesCoordsIndices[3*i+2] = 3*nodesIndices_[i]+2;
    }
    
    Base::adapter_->GetNodesCoords(30, nodesCoordsIndices, nodesCoords);
    Matrix3<TFloat> deformationTensors[4];
    CalcDeformationTensorsAtQuadraturePointsWithLocalBasis(nodesCoords, deformationTensors);
    Matrix3<TFloat> e1, e2, e3, e4;
    
    Base::material_->GetConstitutiveModel()->CalcPK2Stress(deformationTensors[0], e1);
    Base::material_->GetConstitutiveModel()->CalcPK2Stress(deformationTensors[1], e2);
    Base::material_->GetConstitutiveModel()->CalcPK2Stress(deformationTensors[2], e3);
    Base::material_->GetConstitutiveModel()->CalcPK2Stress(deformationTensors[3], e4);
    
    return (e1+e2+e3+e4)/4;
}

CBStatus CBElementSolidT10::GetCauchyStress(Matrix3<TFloat> &cauchyStress) {
    Matrix3<TFloat> deformationTensor;
    
    GetDeformationTensor(deformationTensor);
    
    Matrix3<TFloat> pk2Stress = {0, 0, 0, 0, 0, 0, 0, 0, 0};
    CBStatus rc = Base::material_->GetConstitutiveModel()->CalcPK2Stress(deformationTensor, pk2Stress);
    if (rc != CBStatus::SUCCESS)
        return rc;
    
    Matrix3<TFloat> a;
    Base::adapter_->GetActiveStressTensor(localIndex_, a);
    pk2Stress += a;
    
    cauchyStress = 1.0 / deformationTensor.Det() *  deformationTensor * pk2Stress * deformationTensor.GetTranspose();
    
    return rc;
}

Matrix3<TFloat> CBElementSolidT10::GetDeformationTensorWithGlobalBasis(TFloat l1, TFloat l2, TFloat l3, TFloat l4) {
    TFloat dNdX[30];
    
    CalcShapeFunctionDerivatives(l1, l2, l3, l4, dNdX, true);
    
    Matrix3<TFloat> deformationTensor;
    TFloat *f = deformationTensor.GetArray();
    
    TFloat nodesCoords[30];
    TInt   nodesCoordsIndices[30];
    
    for (unsigned int i = 0; i < 10; i++) {
        nodesCoordsIndices[3*i]   = 3*nodesIndices_[i];
        nodesCoordsIndices[3*i+1] = 3*nodesIndices_[i]+1;
        nodesCoordsIndices[3*i+2] = 3*nodesIndices_[i]+2;
    }
    
    Base::adapter_->GetNodesCoords(30, nodesCoordsIndices, nodesCoords);
    
    for (unsigned int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            TFloat e = 0;
            for (int k = 0; k < 10; k++) {
                e += dNdX[3*k+j] * nodesCoords[3*k+i];
            }
            f[3*i+j] = e;
        }
    }
    
    return deformationTensor;
} // CBElementSolidT10::GetDeformationTensorWithGlobalBasis

Matrix3<TFloat> CBElementSolidT10::GetRightCauchyDeformationTensorWithGlobalBasis(TFloat l1, TFloat l2, TFloat l3,
                                                                                  TFloat l4) {
    Matrix3<TFloat> deformationTensor, C;
    
    deformationTensor = this->GetDeformationTensorWithGlobalBasis(l1, l2, l3, l4);
    C = deformationTensor.GetTranspose() * deformationTensor;
    return C;
}

CBStatus CBElementSolidT10::CalculateLaplacianT4() {
    TInt nodesCoordsIndices[30];
    
    GetNodesCoordsIndices(nodesCoordsIndices);
    TFloat nodesCoords[12];
    Base::adapter_->GetNodesCoords(12, nodesCoordsIndices, nodesCoords);
    TFloat laplacian[16];
    CalcShapeFunctionDerivativesAtQuadraturePoints();
    
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            laplacian[4*i+j] = initialVolume_ *
            (dNdXt4_[3*i]*dNdXt4_[3*j]  + dNdXt4_[3*i+1]*dNdXt4_[3*j+1] + dNdXt4_[3*i+2]*dNdXt4_[3*j+2]);
    Base::adapter_->AddLaplacianEntriesGlobal(4, nodesIndices_.data(), 4, nodesIndices_.data(), laplacian);
    return CBStatus::SUCCESS;
}

CBStatus CBElementSolidT10::CalculateLaplacian() {
    TInt nodesCoordsIndices[30];
    
    GetNodesCoordsIndices(nodesCoordsIndices);
    TFloat nodesCoords[30];
    Base::adapter_->GetNodesCoords(30, nodesCoordsIndices, nodesCoords);
    TFloat laplacian[100];
    CalcShapeFunctionDerivativesAtQuadraturePoints();
    for (int i = 0; i < 10; i++)
        for (int j = 0; j < 10; j++)
            laplacian[10*i+j] = 0;
    
    for (int i = 0; i < 10; i++)
        for (int j = 0; j < 10; j++)
            for (int k = 0; k < 4; k++)
                laplacian[10*i+j] += initialVolume_/4 *
                (dNdXW_[30*k+10*i+0]*dNdXW_[30*k+10*j+0]  + dNdXW_[30*k+10*i+1]*dNdXW_[30*k+10*j+1] + dNdXW_[30*k+10*i+2]*
                 dNdXW_[30*k+10*j+2]);
    
    Base::adapter_->AddLaplacianEntriesGlobal(10, nodesIndices_.data(), 10, nodesIndices_.data(), laplacian);
    return CBStatus::SUCCESS;
}

void CBElementSolidT10::RotateFiberBy(TFloat phi, TFloat theta) {
    Matrix3<TFloat> rotationY = GetRotationY(theta);
    Matrix3<TFloat> rotationZ = GetRotationZ(phi);
    
    Matrix3<TFloat> rotation = rotationY * rotationZ;
    
    for (int i = 0; i < GetNumberOfQuadraturePoints(); i++)
        basisAtQuadraturePoint_[i] = rotation * basisAtQuadraturePoint_[i];
}
