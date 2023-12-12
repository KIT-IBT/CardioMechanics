/*
 * File: CBConstitutiveModelHolzapfel.cpp
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


#include "CBConstitutiveModelHolzapfel.h"
#include "DCCtrlPETSc.h"

void CBConstitutiveModelHolzapfel::Init(ParameterMap *parameters, TInt materialIndex) {
    Base::Init(parameters, materialIndex);
    
    std::stringstream mi;
    mi << materialIndex;
    if ((parameters->IsAvailable("Materials.Mat_" + mi.str() + ".Holzapfel.a") == false) &&
        (parameters->IsAvailable("Materials.Mat_Default.Holzapfel.a") == true) )
        a_ =  parameters->Get<double>("Materials.Mat_Default.Holzapfel.a");
    else
        a_ =  parameters->Get<double>("Materials.Mat_" + mi.str() + ".Holzapfel.a");
    
    if ((parameters->IsAvailable("Materials.Mat_" + mi.str() + ".Holzapfel.b") == false) &&
        (parameters->IsAvailable("Materials.Mat_Default.Holzapfel.b") == true) )
        b_ =  parameters->Get<double>("Materials.Mat_Default.Holzapfel.b");
    else
        b_ =  parameters->Get<double>("Materials.Mat_" + mi.str() + ".Holzapfel.b");
    
    if ((parameters->IsAvailable("Materials.Mat_" + mi.str() + ".Holzapfel.af") == false) &&
        (parameters->IsAvailable("Materials.Mat_Default.Holzapfel.af") == true) )
        af_ =  parameters->Get<double>("Materials.Mat_Default.Holzapfel.af");
    else
        af_ =  parameters->Get<double>("Materials.Mat_" + mi.str() + ".Holzapfel.af");
    
    if ((parameters->IsAvailable("Materials.Mat_" + mi.str() + ".Holzapfel.bf") == false) &&
        (parameters->IsAvailable("Materials.Mat_Default.Holzapfel.bf") == true) )
        bf_ =  parameters->Get<double>("Materials.Mat_Default.Holzapfel.bf");
    else
        bf_ =  parameters->Get<double>("Materials.Mat_" + mi.str() + ".Holzapfel.bf");
    
    if ((parameters->IsAvailable("Materials.Mat_" + mi.str() + ".Holzapfel.as") == false) &&
        (parameters->IsAvailable("Materials.Mat_Default.Holzapfel.as") == true) )
        as_ =  parameters->Get<double>("Materials.Mat_Default.Holzapfel.as");
    else
        as_ =  parameters->Get<double>("Materials.Mat_" + mi.str() + ".Holzapfel.as");
    
    if ((parameters->IsAvailable("Materials.Mat_" + mi.str() + ".Holzapfel.bs") == false) &&
        (parameters->IsAvailable("Materials.Mat_Default.Holzapfel.bs") == true) )
        bs_ =  parameters->Get<double>("Materials.Mat_Default.Holzapfel.bs");
    else
        bs_ =  parameters->Get<double>("Materials.Mat_" + mi.str() + ".Holzapfel.bs");
    
    if ((parameters->IsAvailable("Materials.Mat_" + mi.str() + ".Holzapfel.afs") == false) &&
        (parameters->IsAvailable("Materials.Mat_Default.Holzapfel.afs") == true) )
        afs_ =  parameters->Get<double>("Materials.Mat_Default.Holzapfel.afs");
    else
        afs_ =  parameters->Get<double>("Materials.Mat_" + mi.str() + ".Holzapfel.afs");
    
    if ((parameters->IsAvailable("Materials.Mat_" + mi.str() + ".Holzapfel.bfs") == false) &&
        (parameters->IsAvailable("Materials.Mat_Default.Holzapfel.bfs") == true) )
        bfs_ =  parameters->Get<double>("Materials.Mat_Default.Holzapfel.bfs");
    else
        bfs_ =  parameters->Get<double>("Materials.Mat_" + mi.str() + ".Holzapfel.bfs");
    
    if ((parameters->IsAvailable("Materials.Mat_" + mi.str() + ".Holzapfel.kappa") == false) &&
        (parameters->IsAvailable("Materials.Mat_Default.Holzapfel.kappa") == true) )
        kappa_ =  parameters->Get<double>("Materials.Mat_Default.Holzapfel.kappa");
    else
        kappa_ =  parameters->Get<double>("Materials.Mat_" + mi.str() + ".Holzapfel.kappa");
    
    /// This parameter is optional. k_ determines the slope in the approximation of the Heavyside function
    /// The Heavyside function H(x) was mainly implemented for the 2021 Benchmark problem
    /// If k_ is set to 0, we explicitly take H(x) = 1 if x >= 1, else 0
    if ((parameters->IsAvailable("Materials.Mat_" + mi.str() + ".Holzapfel.k") == false) &&
        (parameters->IsAvailable("Materials.Mat_Default.Holzapfel.k") == true) )
        k_ =  parameters->Get<double>("Materials.Mat_Default.Holzapfel.k", 100);
    else
        k_ =  parameters->Get<double>("Materials.Mat_" + mi.str() + ".Holzapfel.k", 100);
} // CBConstitutiveModelHolzapfel::Init

CBStatus CBConstitutiveModelHolzapfel::CalcEnergy(const Matrix3<TFloat> &deformationTensor, TFloat &energy) {
    TFloat  J = deformationTensor.Det();
    TFloat Jm = pow(J, -2/3);
    
    if (!Base::ignoreCorruptElements_  && (J <= 0))
        return CBStatus::CORRUPT_ELEMENT;
    
    Matrix3<TFloat> rightCauchyGreenTensor = deformationTensor.GetTranspose() * deformationTensor;
    TFloat *C = rightCauchyGreenTensor.GetArray();
    
    TFloat I1   = Jm * rightCauchyGreenTensor.Invariant1();
    TFloat I4f  = C[0]; // f * Cf
    TFloat I4s  = C[4]; // s * Cs
    TFloat I8fs = 0.5 * (C[1] + C[3]); // f * Cs
    
    TFloat vol      = (kappa_ / 4) * (pow(J, 2) - 1.0 - 2 * log(J) );
    TFloat iso      = a_ / (2.0 * b_) * (exp(b_ * (I1 - 3.0)) - 1.0);
    TFloat aniso_f  = Heavyside(I4f) * (af_ / (2.0 * bf_)  * (exp(bf_ * (I4f - 1.0) * (I4f - 1.0)) - 1.0) );
    TFloat aniso_s  = Heavyside(I4s) * (as_ / (2.0 * bs_)  * (exp(bs_ * (I4s - 1.0) * (I4s - 1.0)) - 1.0) );
    TFloat aniso_fs = afs_ / (2.0 * bfs_) * (exp(bfs_* I8fs * I8fs) - 1.0);
    
    energy = vol + iso + aniso_f + aniso_s + aniso_fs;
    
    if (std::isinf(energy))
        return CBStatus::INFINITIVE;
    else if (std::isnan(energy))
        return CBStatus::NOT_A_NUMBER;
    else
        return CBStatus::SUCCESS;
} // CBConstitutiveModelHolzapfel::CalcEnergy

CBStatus CBConstitutiveModelHolzapfel::CalcPK2Stress(const Matrix3<TFloat> &deformationTensor,
                                                     Matrix3<TFloat> &pk2Stress) {
    TFloat J  = deformationTensor.Det();
    TFloat Jm = pow(J, -2/3);
    
    if (!Base::ignoreCorruptElements_  && (J <= 0))
        return CBStatus::CORRUPT_ELEMENT;
    
    Matrix3<TFloat> rightCauchyGreenTensor                  = deformationTensor.GetTranspose() * deformationTensor;
    Matrix3<TFloat> rightCauchyGreenTensor_inv              = rightCauchyGreenTensor.GetInverse();
    TFloat *C = rightCauchyGreenTensor.GetArray();
    
    TFloat I1   = rightCauchyGreenTensor.Invariant1();
    TFloat I4f  = C[0]; // f * Cf
    TFloat I4s  = C[4]; // s * Cs
    TFloat I8fs = 0.5 * (C[1] + C[3]); // f * Cs
    
    // passive contributions to PK2Stress
    Matrix3<TFloat> pk2Vol      = (kappa_ / 2.0) * (J - 1.0/J) * J * rightCauchyGreenTensor_inv;
    Matrix3<TFloat> pk2Iso      = Jm * a_ * exp(b_ * (Jm * I1 - 3.0)) *
    (identity_ - 1.0/3.0 * I1 * rightCauchyGreenTensor_inv);
    Matrix3<TFloat> pk2Aniso_f  = 2.0 * af_ * Heavyside(I4f) * (I4f - 1.0) * exp(bf_ * (I4f - 1.0) * (I4f - 1.0)) * fxf_;
    Matrix3<TFloat> pk2Aniso_s  = 2.0 * as_ * Heavyside(I4s) * (I4s - 1.0) * exp(bs_ * (I4s - 1.0) * (I4s - 1.0)) * sxs_;
    Matrix3<TFloat> pk2Aniso_fs = afs_ * I8fs * exp(bfs_ * I8fs * I8fs) * (fxs_ + sxf_);
    
    pk2Stress = pk2Vol + pk2Iso + pk2Aniso_f + pk2Aniso_s + pk2Aniso_fs;
    
    return CBStatus::SUCCESS;
} // CBConstitutiveModelHolzapfel::CalcPK2Stress

TFloat CBConstitutiveModelHolzapfel::Heavyside(TFloat I4) {
    if (k_ == 0) {
        return (I4 >= 1.0) ? 1.0 : 0.0;
    } else {
        return 1.0 / (1.0 + exp(-k_ * (I4 - 1.0)));
    }
}
