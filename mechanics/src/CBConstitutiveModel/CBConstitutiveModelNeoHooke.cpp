/* -------------------------------------------------------
 
 CBConstitutiveModelNeoHooke.cpp
 
 Ver. 1.0.0
 
 Created:       Tobias Gerach      (13.07.2021)
 Last modified: Tobias Gerach      (13.07.2021)
 
 Institute of Biomedical Engineering
 Karlsruhe Institute of Technology (KIT)
 
 http://www.ibt.kit.edu
 
 Copyright 2000-2009 - All rights reserved.
 
 ------------------------------------------------------ */

/// Taken from Bonet and Wood: NONLINEAR CONTINUUM MECHANICS FOR FINITE ELEMENT ANALYSIS Chapter 5, page 14

#include "CBConstitutiveModelNeoHooke.h"
#include "DCCtrlPETSc.h"

void CBConstitutiveModelNeoHooke::Init(ParameterMap *parameters, TInt materialIndex) {
    Base::Init(parameters, materialIndex);
    
    std::stringstream mi;
    mi << materialIndex;
    
    if ((parameters->IsAvailable("Materials.Mat_" + mi.str() + ".NeoHooke.a") == false) &&
        (parameters->IsAvailable("Materials.Mat_Default.NeoHooke.a") == true) )
        a_ =  parameters->Get<double>("Materials.Mat_Default.NeoHooke.a");
    else
        a_ =  parameters->Get<double>("Materials.Mat_" + mi.str() + ".NeoHooke.a");
    
    if ((parameters->IsAvailable("Materials.Mat_" + mi.str() + ".NeoHooke.k") == false) &&
        (parameters->IsAvailable("Materials.Mat_Default.NeoHooke.k") == true) )
        k_ =  parameters->Get<double>("Materials.Mat_Default.NeoHooke.k");
    else
        k_ =  parameters->Get<double>("Materials.Mat_" + mi.str() + ".NeoHooke.k");
}  // CBConstitutiveModelNeoHooke::Init

CBStatus CBConstitutiveModelNeoHooke::CalcEnergy(const Matrix3<TFloat> &deformationTensor, TFloat &energy) {
    if (!Base::ignoreCorruptElements_) {
        if (deformationTensor.Det() <= 0)
            return CBStatus::CORRUPT_ELEMENT;
    }
    
    TFloat J              = deformationTensor.Det();
    Matrix3<TFloat> C     = deformationTensor.GetTranspose() * deformationTensor;
    TFloat I3_C           = C.Invariant3();
    Matrix3<TFloat> C_dist = pow(I3_C, -1/3) * C;
    
    /// distortional component
    TFloat dist = (a_/2) * (C_dist.Invariant1() - 3);
    
    /// volumetric component
    TFloat vol = (k_/2) * (J - 1) * (J - 1);
    
    energy = dist + vol;
    
    if (std::isinf(energy))
        return CBStatus::INFINITIVE;
    else if (std::isnan(energy))
        return CBStatus::NOT_A_NUMBER;
    else
        return CBStatus::SUCCESS;
}  // CBConstitutiveModelNeoHooke::CalcEnergy

CBStatus CBConstitutiveModelNeoHooke::CalcPK2Stress(const Matrix3<TFloat> &deformationTensor,
                                                    Matrix3<TFloat> &pk2Stress) {
    if (!Base::ignoreCorruptElements_) {
        if (deformationTensor.Det() <= 0)
            return CBStatus::CORRUPT_ELEMENT;
    }
    
    TFloat J              = deformationTensor.Det();
    Matrix3<TFloat> C     = deformationTensor.GetTranspose() * deformationTensor;
    Matrix3<TFloat> C_inv = C.GetInverse();
    TFloat I1_C           = C.Invariant1();
    TFloat I3_C           = C.Invariant3();
    TFloat p              = k_*(J - 1);
    
    pk2Stress = a_ * pow(I3_C, -1/3) * (identity_ - (I1_C * C_inv)/3) + p * J * C_inv;
    
    return CBStatus::SUCCESS;
}  // CBConstitutiveModelNeoHooke::CalcPK2Stress
