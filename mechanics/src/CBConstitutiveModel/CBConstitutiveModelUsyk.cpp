/*
 * File: CBConstitutiveModelUsyk.cpp
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


/* -------------------------------------------------------
 Usyk et al. 2000: Effect of Laminar Orthotropic Myofiber Architecture on Regional Stress and Strain in the Canine
 Left Ventricle (https://intern.ibt.kit.edu/ilibrarian/stable.php?id=14569)
 a   = 880 Pa
 bff = 5
 bss = 6
 bnn = 3
 bfs = 10
 bfn = 2
 bns = 2
 k  >> 0
 ------------------------------------------------------*/

#include "CBConstitutiveModelUsyk.h"
#include "DCCtrlPETSc.h"

void CBConstitutiveModelUsyk::Init(ParameterMap *parameters, TInt materialIndex) {
    Base::Init(parameters, materialIndex);
    
    std::stringstream mi;
    mi << materialIndex;
    
    if ((parameters->IsAvailable("Materials.Mat_" + mi.str() + ".Usyk.a") == false) &&
        (parameters->IsAvailable("Materials.Mat_Default.Usyk.a") == true) )
        a_ =  parameters->Get<double>("Materials.Mat_Default.Usyk.a");
    else
        a_ =  parameters->Get<double>("Materials.Mat_" + mi.str() + ".Usyk.a");
    
    if ((parameters->IsAvailable("Materials.Mat_" + mi.str() + ".Usyk.bff") == false) &&
        (parameters->IsAvailable("Materials.Mat_Default.Usyk.bff") == true) )
        bff_ =  parameters->Get<double>("Materials.Mat_Default.Usyk.bff");
    else
        bff_ =  parameters->Get<double>("Materials.Mat_" + mi.str() + ".Usyk.bff");
    
    if ((parameters->IsAvailable("Materials.Mat_" + mi.str() + ".Usyk.bss") == false) &&
        (parameters->IsAvailable("Materials.Mat_Default.Usyk.bss") == true) )
        bss_ =  parameters->Get<double>("Materials.Mat_Default.Usyk.bss");
    else
        bss_ =  parameters->Get<double>("Materials.Mat_" + mi.str() + ".Usyk.bss");
    
    if ((parameters->IsAvailable("Materials.Mat_" + mi.str() + ".Usyk.bnn") == false) &&
        (parameters->IsAvailable("Materials.Mat_Default.Usyk.bnn") == true) )
        bnn_ =  parameters->Get<double>("Materials.Mat_Default.Usyk.bnn");
    else
        bnn_ =  parameters->Get<double>("Materials.Mat_" + mi.str() + ".Usyk.bnn");
    
    if ((parameters->IsAvailable("Materials.Mat_" + mi.str() + ".Usyk.bfs") == false) &&
        (parameters->IsAvailable("Materials.Mat_Default.Usyk.bfs") == true) )
        bfs_ =  parameters->Get<double>("Materials.Mat_Default.Usyk.bfs");
    else
        bfs_ =  parameters->Get<double>("Materials.Mat_" + mi.str() + ".Usyk.bfs");
    
    if ((parameters->IsAvailable("Materials.Mat_" + mi.str() + ".Usyk.bfn") == false) &&
        (parameters->IsAvailable("Materials.Mat_Default.Usyk.bfn") == true) )
        bfn_ =  parameters->Get<double>("Materials.Mat_Default.Usyk.bfn");
    else
        bfn_ =  parameters->Get<double>("Materials.Mat_" + mi.str() + ".Usyk.bfn");
    
    if ((parameters->IsAvailable("Materials.Mat_" + mi.str() + ".Usyk.bns") == false) &&
        (parameters->IsAvailable("Materials.Mat_Default.Usyk.bns") == true) )
        bns_ =  parameters->Get<double>("Materials.Mat_Default.Usyk.bns");
    else
        bns_ =  parameters->Get<double>("Materials.Mat_" + mi.str() + ".Usyk.bns");
    
    if ((parameters->IsAvailable("Materials.Mat_" + mi.str() + ".Usyk.k") == false) &&
        (parameters->IsAvailable("Materials.Mat_Default.Usyk.k") == true) )
        k_ =  parameters->Get<double>("Materials.Mat_Default.Usyk.k");
    else
        k_ =  parameters->Get<double>("Materials.Mat_" + mi.str() + ".Usyk.k");
    
    if ((parameters->IsAvailable("Materials.Mat_" + mi.str() + ".Usyk.aScale") == false) &&
        (parameters->IsAvailable("Materials.Mat_Default.Usyk.aScale") == true) )
        aScale_ =  parameters->Get<double>("Materials.Mat_Default.Usyk.aScale");
    else
        aScale_ =  parameters->Get<double>("Materials.Mat_" + mi.str() + ".Usyk.aScale", 1);
    
    if ((parameters->IsAvailable("Materials.Mat_" + mi.str() + ".Usyk.bScale") == false) &&
        (parameters->IsAvailable("Materials.Mat_Default.Usyk.bScale") == true) )
        bScale_ =  parameters->Get<double>("Materials.Mat_Default.Usyk.bScale");
    else
        bScale_ =  parameters->Get<double>("Materials.Mat_" + mi.str() + ".Usyk.bScale", 1);
    
    a_   *= aScale_;
    bff_ *= bScale_;
    bss_ *= bScale_;
    bnn_ *= bScale_;
    bfs_ *= bScale_;
    bfn_ *= bScale_;
    bns_ *= bScale_;
}  // CBConstitutiveModelUsyk::Init

CBStatus CBConstitutiveModelUsyk::CalcEnergy(const Matrix3<TFloat> &deformationTensor, TFloat &energy) {
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
    TFloat Q = bff_ * e[0] * e[0] + bss_ * e[4] * e[4] + bnn_ * e[8] * e[8] + bfs_ * (e[1] * e[1] + e[3] * e[3]) + bfn_ *
    (e[2] * e[2] + e[6] * e[6]) + bns_ * (e[5] * e[5] + e[7] * e[7]);
    
    /// iso
    energy = a_/2. * (exp(Q) - 1);
    
    /// vol
    energy += k_/2. * log(J) * log(J);
    
    if (std::isinf(energy))
        return CBStatus::INFINITIVE;
    else if (std::isnan(energy))
        return CBStatus::NOT_A_NUMBER;
    else
        return CBStatus::SUCCESS;
}  // CBConstitutiveModelUsyk::CalcEnergy

CBStatus CBConstitutiveModelUsyk::CalcPK2Stress(const Matrix3<TFloat> &deformationTensor, Matrix3<TFloat> &pk2Stress) {
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
    TFloat *c             = C.GetInverse().GetArray();
    
    /// second Piola Kirchoff stress tensor
    TFloat Q = bff_ * e[0] * e[0] + bss_ * e[4] * e[4] + bnn_ * e[8] * e[8] + bfs_ * (e[1] * e[1] + e[3] * e[3]) + bfn_ *
    (e[2] * e[2] + e[6] * e[6]) + bns_ * (e[5] * e[5] + e[7] * e[7]);
    TFloat expQ = exp(Q);
    
    /// iso
    p[0] = a_ * bff_ * e[0] * expQ;
    p[1] = a_ * bfs_ * e[1] * expQ;
    p[2] = a_ * bfn_ * e[2] * expQ;
    p[3] = a_ * bfs_ * e[3] * expQ;
    p[4] = a_ * bss_ * e[4] * expQ;
    p[5] = a_ * bns_ * e[5] * expQ;
    p[6] = a_ * bfn_ * e[6] * expQ;
    p[7] = a_ * bns_ * e[7] * expQ;
    p[8] = a_ * bnn_ * e[8] * expQ;
    
    /// vol
    p[0] += k_ * log(J) * c[0];
    p[1] += k_ * log(J) * c[1];
    p[2] += k_ * log(J) * c[2];
    p[3] += k_ * log(J) * c[3];
    p[4] += k_ * log(J) * c[4];
    p[5] += k_ * log(J) * c[5];
    p[6] += k_ * log(J) * c[6];
    p[7] += k_ * log(J) * c[7];
    p[8] += k_ * log(J) * c[8];
    
    return CBStatus::SUCCESS;
}  // CBConstitutiveModelUsyk::CalcPK2Stress
