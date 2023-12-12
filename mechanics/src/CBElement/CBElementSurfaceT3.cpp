/*
 * File: CBElementSurfaceT3.cpp
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


#include "CBElementSurfaceT3.h"
#include "CBElementAdapter.h"
#include <functional>

CBElementSurfaceT3::CBElementSurfaceT3(CBElementSurfaceT3 &other) : CBElementCavity(other) {
    initialArea_ = other.initialArea_;
    nodesIndices_ = other.nodesIndices_;
    weight_ = other.weight_;
    basisAtQuadraturePoint_ = other.basisAtQuadraturePoint_;
}

CBElement *CBElementSurfaceT3::New() {
    return new CBElementSurfaceT3;
}

CBElement *CBElementSurfaceT3::Clone() {
    return new CBElementSurfaceT3(*this);
}

void CBElementSurfaceT3::SetWeight(unsigned int i, TFloat w) {
    if (i > 2)
        throw std::runtime_error("CBElementSurfaceT3<T>::SetWeight(unsigned int i,TFloat j) -> i = out of range");
    else
        weight_[i] = w;
}

TFloat CBElementSurfaceT3::GetWeight(unsigned int i) {
    if (i > 2)
        throw std::runtime_error("CBElementSurfaceT3<T>::GetWeight(unsigned int i) -> i = out of range");
    else
        return weight_[i];
}

void CBElementSurfaceT3::SetNodeIndex(unsigned int i, TInt j) {
    if (i > 2)
        throw std::runtime_error("CBElementSurfaceT3<T>::SetNode(unsigned int i,TInt j) -> i = out of range");
    else
        nodesIndices_[i] = j;
}

TInt CBElementSurfaceT3::GetNodeIndex(unsigned int i) {
    if (i > 2)
        throw std::runtime_error("CBElementSurfaceT3<T>::SetNode(unsigned int i,TInt j) -> i = out of range");
    else
        return nodesIndices_[i];
}

void CBElementSurfaceT3::SetBasisAtQuadraturePoint(int i, const Matrix3<TFloat> &basis) {
    if (i == 0) {
        basisAtQuadraturePoint_ = basis;
    } else {
        throw std::runtime_error(
                                 "void CBElementSurfaceT3::void SetBasisAtQuadraturePoint(int i, const Matrix3<TFloat>& basis): i out of range");
    }
}

Matrix3<TFloat> *CBElementSurfaceT3::GetBasisAtQuadraturePoint(int i) {
    {
        if (i == 0) {
            return &basisAtQuadraturePoint_;
        } else {
            throw std::runtime_error(
                                     "Matrix3<TFloat>* CBElementSurfaceT3::void GetBasisAtQuadraturePoint(int i, const Matrix3<TFloat>& basis): i out of range");
        }
    }
}

Vector3<TFloat> CBElementSurfaceT3::GetCentroid() {
    TInt   nodesCoordsIndices[9];
    TFloat nodesCoords[9];
    
    for (unsigned int i = 0; i < 3; i++) {
        nodesCoordsIndices[3 * i]     = 3 * nodesIndices_[i];
        nodesCoordsIndices[3 * i + 1] = 3 * nodesIndices_[i] + 1;
        nodesCoordsIndices[3 * i + 2] = 3 * nodesIndices_[i] + 2;
    }
    
    Base::adapter_->GetNodesCoords(9, nodesCoordsIndices, nodesCoords);
    
    return Vector3<TFloat>(1.0/3.0 * (nodesCoords[0] + nodesCoords[3] + nodesCoords[6]),
                           1.0/3.0 * (nodesCoords[1] + nodesCoords[4] + nodesCoords[7]),
                           1.0/3.0 * (nodesCoords[2] + nodesCoords[5] + nodesCoords[8]));
}

Vector3<TFloat> CBElementSurfaceT3::GetNormalVector() {
    TInt    nodesCoordsIndices[9];
    TFloat  nodesCoords[9];
    TFloat *p = nodesCoords;
    
    for (unsigned int i = 0; i < 3; i++) {
        nodesCoordsIndices[3 * i]     = 3 * nodesIndices_[i];
        nodesCoordsIndices[3 * i + 1] = 3 * nodesIndices_[i] + 1;
        nodesCoordsIndices[3 * i + 2] = 3 * nodesIndices_[i] + 2;
    }
    
    Base::adapter_->GetNodesCoords(9, nodesCoordsIndices, nodesCoords);
    
    Base::adapter_->GetNodesCoords(9, nodesCoordsIndices, nodesCoords);
    Vector3<TFloat> a(p);
    
    p += 3;
    Vector3<TFloat> b(p);
    
    p += 3;
    Vector3<TFloat> c(p);
    
    Vector3<TFloat> n = CrossProduct((a-b), (a-c));
    if (n.Norm() != 0)
        n.Normalize();
    return n;
} // CBElementSurfaceT3::GetNormalVector

TFloat CBElementSurfaceT3::GetArea() {
    TInt    nodesCoordsIndices[9];
    TFloat  nodesCoords[9];
    TFloat *p = nodesCoords;
    
    for (unsigned int i = 0; i < 3; i++) {
        nodesCoordsIndices[3 * i]     = 3 * nodesIndices_[i];
        nodesCoordsIndices[3 * i + 1] = 3 * nodesIndices_[i] + 1;
        nodesCoordsIndices[3 * i + 2] = 3 * nodesIndices_[i] + 2;
    }
    
    Base::adapter_->GetNodesCoords(9, nodesCoordsIndices, nodesCoords);
    Vector3<TFloat> a(p);
    
    p += 3;
    Vector3<TFloat> b(p);
    
    p += 3;
    Vector3<TFloat> c(p);
    
    Vector3<TFloat> n = CrossProduct((a-b), (a-c));
    return 0.5 * n.Norm();
}

TFloat CBElementSurfaceT3::CalcContributionToVolume(const TFloat *referenceCoords) {
    TFloat nodesCoords[9];
    TInt   nodesCoordsIndices[9];
    
    for (unsigned int i = 0; i < 3; i++) {
        nodesCoordsIndices[3 * i]     = 3 * nodesIndices_[i];
        nodesCoordsIndices[3 * i + 1] = 3 * nodesIndices_[i] + 1;
        nodesCoordsIndices[3 * i + 2] = 3 * nodesIndices_[i] + 2;
    }
    
    Base::adapter_->GetNodesCoords(9, nodesCoordsIndices, nodesCoords);
    
    TFloat a[3];
    a[0] = nodesCoords[3] - nodesCoords[0];
    a[1] = nodesCoords[4] - nodesCoords[1];
    a[2] = nodesCoords[5] - nodesCoords[2];
    
    TFloat b[3];
    b[0] = nodesCoords[6] - nodesCoords[0];
    b[1] = nodesCoords[7] - nodesCoords[1];
    b[2] = nodesCoords[8] - nodesCoords[2];
    
    TFloat c[3];
    c[0] = a[1] * b[2] - a[2] * b[1];
    c[1] = a[2] * b[0] - a[0] * b[2];
    c[2] = a[0] * b[1] - a[1] * b[0];
    
    TFloat d[3];
    d[0] = nodesCoords[0] - referenceCoords[0];
    d[1] = nodesCoords[1] - referenceCoords[1];
    d[2] = nodesCoords[2] - referenceCoords[2];
    
    return (d[0] * c[0] + d[1] * c[1] + d[2] * c[2]) / 6.0;
    
    //// Original code using the centroid
    //// However it should be more accurate to directly use the common starting point of vectors a and b --> volume of tetrahedron = 1/6 * volume of parallelepiped
    // a[0] = (nodesCoords[0] + nodesCoords[3] + nodesCoords[6]) / 3.0;
    // a[1] = (nodesCoords[1] + nodesCoords[4] + nodesCoords[7]) / 3.0;
    // a[2] = (nodesCoords[2] + nodesCoords[5] + nodesCoords[8]) / 3.0;
    // return((a[0] * c[0] + a[1] * c[1] + a[2] * c[2]) / 6.0);
} // CBElementSurfaceT3::CalcContributionToVolume

void CBElementSurfaceT3::CalcContributionToVolumeJacobian(TFloat *volumeJacobianEntries,
                                                          TInt *volumeJacobianEntriesIndices) {
    TFloat nodesCoords[9];
    TInt   nodesCoordsIndices[9];
    
    for (unsigned int i = 0; i < 3; i++) {
        nodesCoordsIndices[3 * i]     = 3 * nodesIndices_[i];
        nodesCoordsIndices[3 * i + 1] = 3 * nodesIndices_[i] + 1;
        nodesCoordsIndices[3 * i + 2] = 3 * nodesIndices_[i] + 2;
    }
    
    Base::adapter_->GetNodesCoords(9, nodesCoordsIndices, nodesCoords);
    
    TFloat a[3];
    a[0] = nodesCoords[3] - nodesCoords[0];
    a[1] = nodesCoords[4] - nodesCoords[1];
    a[2] = nodesCoords[5] - nodesCoords[2];
    
    TFloat b[3];
    b[0] = nodesCoords[6] - nodesCoords[0];
    b[1] = nodesCoords[7] - nodesCoords[1];
    b[2] = nodesCoords[8] - nodesCoords[2];
    
    TFloat c[3];
    c[0] = a[1] * b[2] - a[2] * b[1];
    c[1] = a[2] * b[0] - a[0] * b[2];
    c[2] = a[0] * b[1] - a[1] * b[0];
    
    a[0] = (nodesCoords[0] + nodesCoords[3] + nodesCoords[6]) / 3.0;
    a[1] = (nodesCoords[1] + nodesCoords[4] + nodesCoords[7]) / 3.0;
    a[2] = (nodesCoords[2] + nodesCoords[5] + nodesCoords[8]) / 3.0;
    
    double vol = (a[0] * c[0] + a[1] * c[1] + a[2] * c[2]) / 6.0;
    
    TFloat epsilon = Base::adapter_->GetFiniteDifferencesEpsilon();
    
    for (unsigned int i = 0; i < 3; i++)
        for (unsigned int j = 0; j < 3; j++) {
            double tmp = nodesCoords[3 * i + j];
            nodesCoords[3 * i + j] += epsilon;
            
            a[0] = nodesCoords[3] - nodesCoords[0];
            a[1] = nodesCoords[4] - nodesCoords[1];
            a[2] = nodesCoords[5] - nodesCoords[2];
            
            b[0] = nodesCoords[6] - nodesCoords[0];
            b[1] = nodesCoords[7] - nodesCoords[1];
            b[2] = nodesCoords[8] - nodesCoords[2];
            
            c[0] = a[1] * b[2] - a[2] * b[1];
            c[1] = a[2] * b[0] - a[0] * b[2];
            c[2] = a[0] * b[1] - a[1] * b[0];
            
            a[0] = (nodesCoords[0] + nodesCoords[3] + nodesCoords[6]) / 3;
            a[1] = (nodesCoords[1] + nodesCoords[4] + nodesCoords[7]) / 3;
            a[2] = (nodesCoords[2] + nodesCoords[5] + nodesCoords[8]) / 3;
            
            TFloat dv = (((a[0] * c[0] + a[1] * c[1] + a[2] * c[2]) / 6) - vol) / epsilon;
            volumeJacobianEntriesIndices[3 * i + j] = nodesCoordsIndices[3 * i + j];
            volumeJacobianEntries[3 * i + j]        = dv;
            nodesCoords[3 * i + j] = tmp;
        }
} // CBElementSurfaceT3::CalcContributionToVolumeJacobian

void CBElementSurfaceT3::CalcForcesDueToPressure(TFloat pressure, const TInt *nodesCoordsIndices,
                                                 const TFloat *nodesCoords, TFloat *forces) {
    // we use a one point quadrature rule for the integration on the linear triangle element e
    // therefore, the force due to a pressure p at node I is given with
    // f_i = - A_e * sum_i^n[ W * N_i(l1, l2, l3) * p * normalVec ]
    // n = 1
    // W = 1
    // l1 = l2 = l3 = 1/3
    bool bc[9];
    Base::adapter_->GetNodesComponentsBoundaryConditions(9, nodesCoordsIndices, bc);
    
    TFloat area = GetArea();
    Vector3<TFloat> normalVector = GetNormalVector();
    TFloat W = 1.0;
    
    // Shape functions T3 element
    std::function<double(double, double, double)> Ni[3];
    
    Ni[0] = [](double l1, double l2, double l3) {return l1;};
    Ni[1] = [](double l1, double l2, double l3) {return l2;};
    Ni[2] = [](double l1, double l2, double l3) {return l3;};
    
    for (unsigned int i = 0; i < 3; i++) { // iterate over nodes
        Vector3<TFloat> f = -area * pressure * W * Ni[i](1./3., 1./3., 1./3.)*normalVector;
        
        if (bc[3*i+0] != 0)
            forces[3 * i + 0] = 0;
        else
            forces[3 * i + 0] = f.X();
        
        if (bc[3*i+1] != 0)
            forces[3 * i + 1] = 0;
        else
            forces[3 * i + 1] = f.Y();
        
        if (bc[3*i+2] != 0)
            forces[3 * i + 2] = 0;
        else
            forces[3 * i + 2] = f.Z();
    }
} // CBElementSurfaceT3::CalcForcesDueToPressure

void CBElementSurfaceT3::ApplyPressure(TFloat pressure) {
    TFloat nodesCoords[9];
    TInt   nodesCoordsIndices[9];
    
    for (unsigned int i = 0; i < 3; i++) {
        nodesCoordsIndices[3 * i]     = 3 * nodesIndices_[i];
        nodesCoordsIndices[3 * i + 1] = 3 * nodesIndices_[i] + 1;
        nodesCoordsIndices[3 * i + 2] = 3 * nodesIndices_[i] + 2;
    }
    
    Base::adapter_->GetNodesCoords(9, nodesCoordsIndices, nodesCoords);
    
    TFloat forces[9];
    CalcForcesDueToPressure(pressure, nodesCoordsIndices, nodesCoords, forces);
    
    Base::adapter_->AddNodalForcesComponents(9, nodesCoordsIndices, forces);
}

Triangle<TFloat> CBElementSurfaceT3::GetTriangle() {
    TInt   nodesCoordsIndices[9];
    TFloat nodesCoords[9];
    
    for (unsigned int i = 0; i < 3; i++) {
        nodesCoordsIndices[3 * i]     = 3 * nodesIndices_[i];
        nodesCoordsIndices[3 * i + 1] = 3 * nodesIndices_[i] + 1;
        nodesCoordsIndices[3 * i + 2] = 3 * nodesIndices_[i] + 2;
    }
    
    Base::adapter_->GetNodesCoords(9, nodesCoordsIndices, nodesCoords);
    return Triangle<TFloat>(nodesCoords);
}

void CBElementSurfaceT3::CalcPressureJacobian(TFloat pressure) {
    TFloat nodesCoords[9];
    TInt   nodesCoordsIndices[9];
    
    for (unsigned int i = 0; i < 3; i++) {
        nodesCoordsIndices[3 * i]     = 3 * nodesIndices_[i];
        nodesCoordsIndices[3 * i + 1] = 3 * nodesIndices_[i] + 1;
        nodesCoordsIndices[3 * i + 2] = 3 * nodesIndices_[i] + 2;
    }
    
    Base::adapter_->GetNodesCoords(9, nodesCoordsIndices, nodesCoords);
    
    TFloat epsilon = Base::adapter_->GetFiniteDifferencesEpsilon();
    
    TFloat dfdx[9*9];
    
    bool bc[9];
    Base::adapter_->GetNodesComponentsBoundaryConditions(9, nodesCoordsIndices, bc);
    
    for (unsigned int i = 0; i < 3; i++) // iter over nodes
        for (unsigned int j = 0; j < 3; j++) { // iter over coordinates x,y,z
            TFloat f1[9];
            TFloat f2[9];
            
            double tmp = nodesCoords[3 * i + j];
            nodesCoords[3 * i + j] += epsilon;
            CalcForcesDueToPressure(pressure, nodesCoordsIndices, nodesCoords, f1);
            nodesCoords[3 * i + j] -= 2.0*epsilon;
            CalcForcesDueToPressure(pressure, nodesCoordsIndices, nodesCoords, f2);
            
            for (unsigned int k = 0; k < 3; k++)
                for (unsigned int l = 0; l < 3; l++) {
                    if ((bc[3*k+l] == 0) && (bc[3*i+j] == 0))
                        dfdx[9 * (3 * k + l) + (3 * i + j)] = (f1[3*k+l] - f2[3*k+l])/(2*epsilon);
                    else
                        dfdx[9 * (3 * k + l) + (3 * i + j)] = 0;
                }
            nodesCoords[3 * i + j] = tmp;
        }
    
    Base::adapter_->AddNodalForcesJacobianEntries(9, nodesCoordsIndices, 9, nodesCoordsIndices, dfdx);
} // CBElementSurfaceT3::CalcPressureJacobian

void CBElementSurfaceT3::GetAreaNormalAndIndices(TFloat *c, TInt *nodesCoordsIndices) {
    TFloat nodesCoords[9];
    
    for (unsigned int i = 0; i < 3; i++) {
        nodesCoordsIndices[3 * i]     = 3 * nodesIndices_[i];
        nodesCoordsIndices[3 * i + 1] = 3 * nodesIndices_[i] + 1;
        nodesCoordsIndices[3 * i + 2] = 3 * nodesIndices_[i] + 2;
    }
    
    Base::adapter_->GetNodesCoords(9, nodesCoordsIndices, nodesCoords);
    
    TFloat a[3];
    a[0] = nodesCoords[3] - nodesCoords[0];
    a[1] = nodesCoords[4] - nodesCoords[1];
    a[2] = nodesCoords[5] - nodesCoords[2];
    
    TFloat b[3];
    b[0] = nodesCoords[6] - nodesCoords[0];
    b[1] = nodesCoords[7] - nodesCoords[1];
    b[2] = nodesCoords[8] - nodesCoords[2];
    
    c[0] = (a[1] * b[2] - a[2] * b[1]) / 2.0;
    c[1] = (a[2] * b[0] - a[0] * b[2]) / 2.0;
    c[2] = (a[0] * b[1] - a[1] * b[0]) / 2.0;
}
