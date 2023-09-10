/*
 *  CBElementSolid.h
 *  CardioMechanics
 *
 *  Created by Thomas Fritz on 31.03.11.
 *  Copyright 2011 Institute for Biomedical Engineering, Karlsruhe Institute of Technology (KIT). All rights reserved.
 *
 */

#ifndef CB_ELEMENT_SOLID_H
#define CB_ELEMENT_SOLID_H

#include "CBElement.h"
#include <typeinfo>

using namespace math_pack;

class CBElementSolid : public CBElement {
public:
    CBElementSolid() : CBElement() {
        SetProperty(ElementProperty::Solid);
    }
    
    CBElementSolid(CBElementSolid &other);
    
    virtual ~CBElementSolid() {}
    
    // ------- Pure virtual functions -------
    virtual CBStatus CalcNodalForcesJacobian() = 0;
    
    virtual CBStatus CalcNodalForces()         = 0;
    
    virtual CBStatus SetNodalForcesToZeroIfElementIsDefect() {
        throw std::runtime_error(
                                 "CBElementSolid::CalculateLaplacian(): This function is not implemented for the requested child class:");
    }
    
    virtual CBStatus SetNodalForcesJacobianToZeroIfElementIsDefect() {
        throw std::runtime_error(
                                 "CBElementSolid::CalculateLaplacian(): This function is not implemented for the requested child class:");
    }
    
    virtual TFloat GetDeformationEnergy() = 0;
    virtual CBStatus GetCauchyStress(Matrix3<TFloat> &cauchyStress) = 0;
    virtual CBStatus GetDeformationTensor(Matrix3<TFloat> &f) = 0;
    
    virtual Matrix3<TFloat> GetPK2Stress() {
        throw std::runtime_error(
                                 "CBElementSolid::GetPK2Stress(): This function is not implemented for the requested child class:");
    }
    
    virtual CBStatus GetDeformationTensorAtQuadraturePoints(Matrix3<TFloat> *f) {return CBStatus::SUCCESS;}
    
    virtual CBStatus CalcNodalForcesAndJacobian() = 0;
    virtual CBStatus CalcConsistentMassMatrix()   = 0;
    virtual CBStatus CalcLumpedMassMatrix()       = 0;
    
    virtual bool IsElementInverted() {return false;}
    
    virtual void CheckNodeSorting()    = 0;
    virtual void UpdateShapeFunctions() = 0;
    virtual TFloat GetVolume()          = 0;
    
    bool GetStatus() {return isDefect_;}
    
    virtual TFloat *GetShapeFunctionsDerivatives() = 0;
    
    virtual void CalcShapeFunctionDerivatives(TFloat l1, TFloat l2, TFloat l3, TFloat l4, TFloat *dNdX,
                                              bool useReferenceNodes = 0) {}
    
    virtual void RotateFiberBy(TFloat phi, TFloat theta) {}
    
    // --------------------------------------
    virtual bool CheckFiberLength() override;
    
    virtual CBStatus CalculateLaplacian() {
        throw std::runtime_error(
                                 "CBElementSolid::CalculateLaplacian(): This function is not implemented for the requested child class:");
    }
    
    virtual CBStatus CalcNodalForcesDueActiveStress() {
        throw std::runtime_error(
                                 "CBElementSolid::CalcNodalForcesDueActiveStress(): This function is not implemented for the requested child class:");
    }
    
    virtual Matrix3<TFloat> GetRightCauchyDeformationTensorWithGlobalBasis(TFloat l1, TFloat l2, TFloat l3, TFloat l4) {
        throw std::runtime_error(
                                 "CBElementSolid::GetRightCauchyDeformationTensorWithGlobalBasis(TFloat l1, TFloat l2, TFloat l3, TFloat l4): This function is not implemented for the requested child class:");
    }
    
    TFloat GetInitialVolume() {return initialVolume_;}
    
    static void RepairDeformationTensorIfInverted(Matrix3<TFloat> &deformationTensor);
    
protected:
    TFloat initialVolume_;
    bool isDefect_ = false;
    Matrix3<TFloat> targetStrain = Matrix3<TFloat>::Identity();
    
private:
    void operator=(const CBElementSolid &) = delete;
};
#endif // ifndef CB_ELEMENT_SOLID_H
