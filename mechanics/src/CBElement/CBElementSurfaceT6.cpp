/*
 *  CBElementSurfaceT6.cpp
 *  CardioMechanics
 *
 *  Created by Thomas Fritz on 13.09.11.
 *  Copyright 2011 Institute for Biomedical Engineering, Karlsruhe Institute of Technology (KIT). All rights reserved.
 *
 */

#include "CBElementSurfaceT6.h"
#include "CBElementAdapter.h"
#include <functional>


CBElementSurfaceT6::CBElementSurfaceT6(CBElementSurfaceT6& other) : CBElementCavity(other) {
    nodesIndices_ = other.nodesIndices_;
    basisAtQuadraturePoint_ = other.basisAtQuadraturePoint_;
}

CBElement* CBElementSurfaceT6::New()
{
    return(new CBElementSurfaceT6);
}

CBElement* CBElementSurfaceT6::Clone() {
    return new CBElementSurfaceT6(*this);
}


void CBElementSurfaceT6::SetNodeIndex(unsigned int i, TInt j)
{
    if(i > 5)
        throw std::runtime_error("CBElementSurfaceT6<T>::SetNode(unsigned int i,TInt j) -> i = out of range");
    else
        nodesIndices_[i] = j;
}


TInt CBElementSurfaceT6::GetNodeIndex(unsigned int i)
{
    if(i > 5)
        throw std::runtime_error("CBElementSurfaceT6<T>::SetNode(unsigned int i,TInt j) -> i = out of range");
    else
        return(nodesIndices_[i]);
}

void CBElementSurfaceT6::SetBasisAtQuadraturePoint(int i, const Matrix3<TFloat>& basis)
{
    if(i == 0)
        basisAtQuadraturePoint_ = basis;
    else
        throw std::runtime_error("void CBElementSurfaceT6::void SetBasisAtQuadraturePoint(int i, const Matrix3<TFloat>& basis): i out of range");
}


Matrix3<TFloat>* CBElementSurfaceT6::GetBasisAtQuadraturePoint(int i)
{
    {
        if(i == 0)
            return(&basisAtQuadraturePoint_);
        else
            throw std::runtime_error("Matrix3<TFloat>* CBElementSurfaceT6::void GetBasisAtQuadraturePoint(int i, const Matrix3<TFloat>& basis): i out of range");
    }
}


Vector3<TFloat> CBElementSurfaceT6::GetCentroid()
{
    TInt   nodesCoordsIndices[9];
    TFloat nodesCoords[9];
    
    for(unsigned int i = 0; i < 3; i++)
    {
        nodesCoordsIndices[3 * i]     = 3 * nodesIndices_[i];
        nodesCoordsIndices[3 * i + 1] = 3 * nodesIndices_[i] + 1;
        nodesCoordsIndices[3 * i + 2] = 3 * nodesIndices_[i] + 2;
    }
    
    Base::adapter_->GetNodesCoords(9, nodesCoordsIndices, nodesCoords);
    
    return(Vector3<TFloat>(1.0/3.0 * (nodesCoords[0] + nodesCoords[3] + nodesCoords[6]),
                           1.0/3.0 * (nodesCoords[1] + nodesCoords[4] + nodesCoords[7]),
                           1.0/3.0 * (nodesCoords[2] + nodesCoords[5] + nodesCoords[8])));
}


Matrix3<TFloat> CBElementSurfaceT6::GetNormalVectorAtQuadraturePoints(const TFloat* nodesCoords)
{
    const TFloat* p = nodesCoords;
    
    Vector3<TFloat> n1(&p[0]);
    Vector3<TFloat> n2(&p[3]);
    Vector3<TFloat> n3(&p[6]);
    Vector3<TFloat> n4(&p[9]);
    Vector3<TFloat> n5(&p[12]);
    Vector3<TFloat> n6(&p[15]);
    Vector3<TFloat> n7(&p[18]);
    
    Triangle<TFloat> t1(n1,n4,n6);
    Triangle<TFloat> t2(n4,n5,n6);
    Triangle<TFloat> t3(n4,n2,n5);
    Triangle<TFloat> t4(n6,n5,n3);
    
    Matrix3<TFloat> normalVectors;
    
    Vector3<TFloat> nv1 = t1.GetNormalVector()+t2.GetNormalVector();
    nv1.Normalize();
    Vector3<TFloat> nv2 = t3.GetNormalVector()+t2.GetNormalVector();
    nv2.Normalize();
    Vector3<TFloat> nv3 = t4.GetNormalVector()+t2.GetNormalVector();
    nv3.Normalize();
    
    normalVectors.SetCol(0, nv1);
    normalVectors.SetCol(1, nv2);
    normalVectors.SetCol(2, nv3);
    return normalVectors;
}

TFloat CBElementSurfaceT6::GetArea()
{
    TInt    nodesCoordsIndices[18];
    TFloat  nodesCoords[18];
    
    for(unsigned int i = 0; i < 6; i++)
    {
        nodesCoordsIndices[3 * i]     = 3 * nodesIndices_[i];
        nodesCoordsIndices[3 * i + 1] = 3 * nodesIndices_[i] + 1;
        nodesCoordsIndices[3 * i + 2] = 3 * nodesIndices_[i] + 2;
    }
    
    Base::adapter_->GetNodesCoords(18,nodesCoordsIndices, nodesCoords);
    
    return GetArea(nodesCoords);
}

TFloat CBElementSurfaceT6::GetArea(const TFloat* nodesCoords)
{
    const TFloat* p = nodesCoords;
    
    Vector3<TFloat> n1(&p[0]);
    Vector3<TFloat> n2(&p[3]);
    Vector3<TFloat> n3(&p[6]);
    Vector3<TFloat> n4(&p[9]);
    Vector3<TFloat> n5(&p[12]);
    Vector3<TFloat> n6(&p[15]);
    
    Triangle<TFloat> t1(n1,n4,n6);
    Triangle<TFloat> t2(n4,n5,n6);
    Triangle<TFloat> t3(n4,n2,n5);
    Triangle<TFloat> t4(n6,n5,n3);
    
    return t1.GetArea()+t2.GetArea()+t3.GetArea()+t4.GetArea();
}

//TFloat CBElementSurfaceT6::CalcContributionToVolume()
//{
//    TInt    nodesCoordsIndices[18];
//    TFloat  nodesCoords[18];
//    
//    for(unsigned int i = 0; i < 6; i++)
//    {
//        nodesCoordsIndices[3 * i]     = 3 * nodesIndices_[i];
//        nodesCoordsIndices[3 * i + 1] = 3 * nodesIndices_[i] + 1;
//        nodesCoordsIndices[3 * i + 2] = 3 * nodesIndices_[i] + 2;
//    }
//    
//    Base::adapter_->GetNodesCoords(18,nodesCoordsIndices, nodesCoords);
//    const TFloat* p = nodesCoords;
//    
//    Vector3<TFloat> n1(&p[0]);
//    Vector3<TFloat> n2(&p[3]);
//    Vector3<TFloat> n3(&p[6]);
//    Vector3<TFloat> n4(&p[9]);
//    Vector3<TFloat> n5(&p[12]);
//    Vector3<TFloat> n6(&p[15]);
//    
//    Triangle<TFloat> t1(n1,n4,n6);
//    Triangle<TFloat> t2(n4,n5,n6);
//    Triangle<TFloat> t3(n4,n2,n5);
//    Triangle<TFloat> t4(n6,n5,n3);
//    
//    return 1.0/3.0*(t1.GetArea()*t1.GetNormalVector())*t1.GetCentroid()+(t2.GetArea()*t2.GetNormalVector())*t2.GetCentroid()+(t3.GetArea()*t3.GetNormalVector())*t3.GetCentroid()+(t4.GetArea()*t4.GetNormalVector())*t4.GetCentroid();
//
//}

TFloat CBElementSurfaceT6::CalcContributionToVolume(const TFloat* referenceCoords)
{
    TInt    nodesCoordsIndices[18];
    for(unsigned int i = 0; i < 6; i++)
    {
        nodesCoordsIndices[3 * i]     = 3 * nodesIndices_[i];
        nodesCoordsIndices[3 * i + 1] = 3 * nodesIndices_[i] + 1;
        nodesCoordsIndices[3 * i + 2] = 3 * nodesIndices_[i] + 2;
    }
    
    TFloat  nodesCoords[18];
    Base::adapter_->GetNodesCoords(18, nodesCoordsIndices, nodesCoords);
    
    TFloat a1[3] = {nodesCoords[ 9]-nodesCoords[0], nodesCoords[10]-nodesCoords[1], nodesCoords[11]-nodesCoords[2]};
    TFloat a2[3] = {nodesCoords[ 3]-nodesCoords[0], nodesCoords[ 4]-nodesCoords[1], nodesCoords[ 5]-nodesCoords[2]};
    TFloat a3[3] = {nodesCoords[12]-nodesCoords[0], nodesCoords[13]-nodesCoords[1], nodesCoords[14]-nodesCoords[2]};
    TFloat a4[3] = {nodesCoords[ 6]-nodesCoords[0], nodesCoords[ 7]-nodesCoords[1], nodesCoords[ 8]-nodesCoords[2]};
    TFloat a5[3] = {nodesCoords[15]-nodesCoords[0], nodesCoords[16]-nodesCoords[1], nodesCoords[17]-nodesCoords[2]};
    
    TFloat c1[3];
    c1[0] = a1[1] * a2[2] - a1[2] * a2[1];
    c1[1] = a1[2] * a2[0] - a1[0] * a2[2];
    c1[2] = a1[0] * a2[1] - a1[1] * a2[0];
    
    TFloat c2[3];
    c2[0] = a2[1] * a3[2] - a2[2] * a3[1];
    c2[1] = a2[2] * a3[0] - a2[0] * a3[2];
    c2[2] = a2[0] * a3[1] - a2[1] * a3[0];
    
    TFloat c3[3];
    c3[0] = a3[1] * a4[2] - a3[2] * a4[1];
    c3[1] = a3[2] * a4[0] - a3[0] * a4[2];
    c3[2] = a3[0] * a4[1] - a3[1] * a4[0];
    
    TFloat c4[3];
    c4[0] = a4[1] * a5[2] - a4[2] * a5[1];
    c4[1] = a4[2] * a5[0] - a4[0] * a5[2];
    c4[2] = a4[0] * a5[1] - a4[1] * a5[0];
    
    TFloat d[3];
    d[0] = nodesCoords[0] - referenceCoords[0];
    d[1] = nodesCoords[1] - referenceCoords[1];
    d[2] = nodesCoords[2] - referenceCoords[2];
    
    TFloat v1 = d[0] * c1[0] + d[1] * c1[1] + d[2] * c1[2];
    TFloat v2 = d[0] * c2[0] + d[1] * c2[1] + d[2] * c2[2];
    TFloat v3 = d[0] * c3[0] + d[1] * c3[1] + d[2] * c3[2];
    TFloat v4 = d[0] * c4[0] + d[1] * c4[1] + d[2] * c4[2];
    
    return (v1+v2+v3+v4)/6.0;
}

void CBElementSurfaceT6::CalcContributionToVolumeJacobian(TFloat* volumeJacobianEntries, TInt* volumeJacobianEntriesIndices)
{
    throw std::runtime_error("CBElementSurfaceT6::CalcContributionToVolumeJacobian is not implemented! Either try to implement or use T3 instead.");
}

void CBElementSurfaceT6::CalcForcesDueToPressure(TFloat pressure, const TInt* nodesCoordsIndices, const TFloat* nodesCoords, TFloat* forces)
{
    
    bool bc[18];
    Base::adapter_->GetNodesComponentsBoundaryConditions(18, nodesCoordsIndices, bc);
    
    TFloat area = GetArea(nodesCoords);
    Matrix3<TFloat> normalVectors = GetNormalVectorAtQuadraturePoints(nodesCoords);
    
    std::function<double (double,double,double)> Ni[6];
    
    Ni[0] = [](double l1,double l2, double l3){return l1*(2*l1 - 1);};
    Ni[1] = [](double l1,double l2, double l3){return l2*(2*l2 - 1);};
    Ni[2] = [](double l1,double l2, double l3){return l3*(2*l3 - 1);};
    Ni[3] = [](double l1,double l2, double l3){return 4*l1*l2;};
    Ni[4] = [](double l1,double l2, double l3){return 4*l1*l3;};
    Ni[5] = [](double l1,double l2, double l3){return 4*l2*l3;};
    
    for(unsigned int i = 0; i < 6; i++)
    {
        Vector3<TFloat> f = -area * pressure * 1.0/3.0 * (Ni[i](2.0/3.0,1.0/6.0,1.0/6.0)*normalVectors.GetCol(0) + Ni[i](1.0/6.0,2.0/3.0,1.0/6.0)*normalVectors.GetCol(1) + Ni[i](1.0/6.0,1.0/6.0,2.0/3.0)*normalVectors.GetCol(2));
        
        if (bc[3*i+0] !=0)
            forces[3 * i + 0] = 0;
        else
            forces[3 * i + 0] = f.X();
        
        if (bc[3*i+1] !=0)
            forces[3 * i + 1] = 0;
        else
            forces[3 * i + 1] = f.Y();
        
        if (bc[3*i+2] !=0)
            forces[3 * i + 2] = 0;
        else
            forces[3 * i + 2] = f.Z();
    }
}

void CBElementSurfaceT6::ApplyPressure(TFloat pressure)
{
    TInt    nodesCoordsIndices[18];
    TFloat  nodesCoords[18];
    
    for(unsigned int i = 0; i < 6; i++)
    {
        nodesCoordsIndices[3 * i]     = 3 * nodesIndices_[i];
        nodesCoordsIndices[3 * i + 1] = 3 * nodesIndices_[i] + 1;
        nodesCoordsIndices[3 * i + 2] = 3 * nodesIndices_[i] + 2;
    }
    
    Base::adapter_->GetNodesCoords(18,nodesCoordsIndices, nodesCoords);
    TFloat forces[18];
    CalcForcesDueToPressure(pressure, nodesCoordsIndices, nodesCoords, forces);
    Base::adapter_->AddNodalForcesComponents(18, nodesCoordsIndices, forces);
}

Triangle<TFloat> CBElementSurfaceT6::GetTriangle()
{
    TInt   nodesCoordsIndices[9];
    TFloat nodesCoords[9];
    
    for(unsigned int i = 0; i < 3; i++)
    {
        nodesCoordsIndices[3 * i]     = 3 * nodesIndices_[i];
        nodesCoordsIndices[3 * i + 1] = 3 * nodesIndices_[i] + 1;
        nodesCoordsIndices[3 * i + 2] = 3 * nodesIndices_[i] + 2;
    }
    
    Base::adapter_->GetNodesCoords(9, nodesCoordsIndices, nodesCoords);
    return(Triangle<TFloat>(nodesCoords));
}

void CBElementSurfaceT6::CalcPressureJacobian(TFloat pressure)
{
    TInt    nodesCoordsIndices[18];
    TFloat  nodesCoords[18];
    
    for(unsigned int i = 0; i < 6; i++)
    {
        nodesCoordsIndices[3 * i + 0] = 3 * nodesIndices_[i];
        nodesCoordsIndices[3 * i + 1] = 3 * nodesIndices_[i] + 1;
        nodesCoordsIndices[3 * i + 2] = 3 * nodesIndices_[i] + 2;
    }
    
    Base::adapter_->GetNodesCoords(18,nodesCoordsIndices, nodesCoords);
    
    TFloat epsilon = Base::adapter_->GetFiniteDifferencesEpsilon();
    
    TFloat dfdx[18*18];
    
    bool bc[18];
    Base::adapter_->GetNodesComponentsBoundaryConditions(18, nodesCoordsIndices, bc);
    
    for(unsigned int i = 0; i < 6; i++)
        for(unsigned int j = 0; j < 3; j++)
        {
            TFloat f1[18];
            TFloat f2[18];
            
            double tmp = nodesCoords[3 * i + j];
            nodesCoords[3 * i + j] += epsilon;
            CalcForcesDueToPressure(pressure, nodesCoordsIndices, nodesCoords, f1);
            nodesCoords[3 * i + j] -= 2.0*epsilon;
            CalcForcesDueToPressure(pressure, nodesCoordsIndices, nodesCoords, f2);
            
            for(unsigned int k = 0; k < 6; k++)
                for(unsigned int l = 0; l < 3; l++)
                {
                    if(bc[3*k+l] == 0 && bc[3*i+j]==0)
                        dfdx[18 * (3 * k + l) + (3 * i + j)] = (f1[3*k+l] - f2[3*k+l])/(2*epsilon);
                    else
                        dfdx[18 * (3 * k + l) + (3 * i + j)] = 0;
                }
            nodesCoords[3 * i + j] = tmp;
        }
    
    Base::adapter_->AddNodalForcesJacobianEntries(18, nodesCoordsIndices, 18, nodesCoordsIndices, dfdx);
}
