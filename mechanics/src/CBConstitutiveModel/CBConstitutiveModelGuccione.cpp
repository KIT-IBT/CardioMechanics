/*
 *  CBConstitutiveModelGuccione.cpp
 *  CardioMechanics
 *
 *  Created by Thomas Fritz on 14.09.11.
 *  Copyright 2011 Institute for Biomedical Engineering, Karlsruhe Institute of Technology (KIT). All rights reserved.
 *
 */


#include "CBConstitutiveModelGuccione.h"
#include "DCCtrlPETSc.h"
#include <algorithm>

void CBConstitutiveModelGuccione::Init(ParameterMap *parameters, TInt materialIndex) {
    Base::Init(parameters, materialIndex);
    
    std::stringstream mi;
    mi << materialIndex;
    if ((parameters->IsAvailable("Materials.Mat_" + mi.str() + ".Guccione.C") == false) &&
        (parameters->IsAvailable("Materials.Mat_Default.Guccione.C") == true) )
        C_ =  parameters->Get<double>("Materials.Mat_Default.Guccione.C");
    else
        C_ =  parameters->Get<double>("Materials.Mat_" + mi.str() + ".Guccione.C");
    
    if ((parameters->IsAvailable("Materials.Mat_" + mi.str() + ".Guccione.bf") == false) &&
        (parameters->IsAvailable("Materials.Mat_Default.Guccione.bf") == true) )
        bf_ =  parameters->Get<double>("Materials.Mat_Default.Guccione.bf");
    else
        bf_ =  parameters->Get<double>("Materials.Mat_" + mi.str() + ".Guccione.bf");
    
    if ((parameters->IsAvailable("Materials.Mat_" + mi.str() + ".Guccione.bt") == false) &&
        (parameters->IsAvailable("Materials.Mat_Default.Guccione.bt") == true) )
        bt_ =  parameters->Get<double>("Materials.Mat_Default.Guccione.bt");
    else
        bt_ =  parameters->Get<double>("Materials.Mat_" + mi.str() + ".Guccione.bt");
    
    if ((parameters->IsAvailable("Materials.Mat_" + mi.str() + ".Guccione.bfs") == false) &&
        (parameters->IsAvailable("Materials.Mat_Default.Guccione.bfs") == true) )
        bfs_ =  parameters->Get<double>("Materials.Mat_Default.Guccione.bfs");
    else
        bfs_ =  parameters->Get<double>("Materials.Mat_" + mi.str() + ".Guccione.bfs");
    
    if ((parameters->IsAvailable("Materials.Mat_" + mi.str() + ".Guccione.K") == false) &&
        (parameters->IsAvailable("Materials.Mat_Default.Guccione.K") == true) )
        k_ =  parameters->Get<double>("Materials.Mat_Default.Guccione.K");
    else
        k_ =  parameters->Get<double>("Materials.Mat_" + mi.str() + ".Guccione.K");
    
    if ((parameters->IsAvailable("Materials.Mat_" + mi.str() + ".Guccione.aScale") == false) &&
        (parameters->IsAvailable("Materials.Mat_Default.Guccione.aScale") == true) )
        aScale_ =  parameters->Get<double>("Materials.Mat_Default.Guccione.aScale");
    else
        aScale_ =  parameters->Get<double>("Materials.Mat_" + mi.str() + ".Guccione.aScale", 1);
    
    if ((parameters->IsAvailable("Materials.Mat_" + mi.str() + ".Guccione.bScale") == false) &&
        (parameters->IsAvailable("Materials.Mat_Default.Guccione.bScale") == true) )
        bScale_ =  parameters->Get<double>("Materials.Mat_Default.Guccione.bScale");
    else
        bScale_ =  parameters->Get<double>("Materials.Mat_" + mi.str() + ".Guccione.bScale", 1);
    
    C_  *= aScale_;
    bf_ *= bScale_;
    bt_ *= bScale_;
    bfs_ *= bScale_;
} // CBConstitutiveModelGuccione::Init

CBStatus CBConstitutiveModelGuccione::CalcEnergy(const Matrix3<TFloat> &deformationTensor, TFloat &energy) {
    if (!Base::ignoreCorruptElements_) {
        if (deformationTensor.Det() <= 0)
            return CBStatus::CORRUPT_ELEMENT;
    }
    
    /// multiplicative decomposition of F into volume-changing (vol) and volume-preserving (iso) parts.
    TFloat J              = deformationTensor.Det();
    Matrix3<TFloat> C     = deformationTensor.GetTranspose() * deformationTensor;
    Matrix3<TFloat> E     = 0.5 * (C - identity_);
    TFloat *e             = E.GetArray();
    
    /// strain energy function W
    TFloat Q = bf_ * e[0] * e[0] + bt_ * (e[4] * e[4] + e[8] * e[8] + e[5] * e[5] + e[7] * e[7]) + bfs_ *
    (e[1] * e[1] + e[3] * e[3] + e[2] * e[2] + e[6] * e[6]);
    
    /// iso
    energy = C_/2.  * (exp(Q) - 1.0);
    
    /// vol
    energy += k_/2. * (J - 1) * (J - 1);
    
    if (std::isinf(energy))
        return CBStatus::INFINITIVE;
    else if (std::isnan(energy))
        return CBStatus::NOT_A_NUMBER;
    else
        return CBStatus::SUCCESS;
} // CBConstitutiveModelGuccione::CalcEnergy

CBStatus CBConstitutiveModelGuccione::CalcPK2Stress(const Matrix3<TFloat> &deformationTensor,
                                                    Matrix3<TFloat> &pk2Stress) {
    if (!Base::ignoreCorruptElements_) {
        if (deformationTensor.Det() <= 0)
            return CBStatus::CORRUPT_ELEMENT;
    }
    
    /// multiplicative decomposition of F into volume-changing (vol) and volume-preserving (iso) parts.
    TFloat J              = deformationTensor.Det();
    Matrix3<TFloat> C     = deformationTensor.GetTranspose() * deformationTensor;
    Matrix3<TFloat> E     = 0.5 * (C - identity_);
    TFloat *e             = E.GetArray();
    TFloat *p             = pk2Stress.GetArray();
    Matrix3<TFloat> C_inv = C.GetInverse();
    
    /// strain energy function W
    TFloat Q = bf_ * e[0] * e[0] + bt_ * (e[4] * e[4] + e[8] * e[8] + e[5] * e[5] + e[7] * e[7]) + bfs_ *
    (e[1] * e[1] + e[3] * e[3] + e[2] * e[2] + e[6] * e[6]);
    TFloat expQ = exp(Q);
    
    /// iso
    p[0] = C_ * bf_ * e[0] * expQ;
    p[1] = C_ * bfs_ * e[1] * expQ;
    p[2] = C_ * bfs_ * e[2] * expQ;
    p[3] = C_ * bfs_ * e[3] * expQ;
    p[4] = C_ * bt_ * e[4] * expQ;
    p[5] = C_ * bt_ * e[5] * expQ;
    p[6] = C_ * bfs_ * e[6] * expQ;
    p[7] = C_ * bt_ * e[7] * expQ;
    p[8] = C_ * bt_ * e[8] * expQ;
    
    /// vol
    pk2Stress += k_ * (J - 1) * J * C_inv;
    
    return CBStatus::SUCCESS;
} // CBConstitutiveModelGuccione::CalcPK2Stress
