/*
 *  CBConstitutiveModelMooneyRivlin.h.cpp
 *  CardioMechanics
 *
 *  Created by Thomas Fritz on 14.09.11.
 *  Copyright 2011 Institute for Biomedical Engineering, Karlsruhe Institute of Technology (KIT). All rights reserved.
 *
 */

#include "CBConstitutiveModelMooneyRivlin.h"



void CBConstitutiveModelMooneyRivlin::Init(ParameterMap* parameters, TInt materialIndex)
{
    Base::Init(parameters,materialIndex);
    
    std::stringstream mi;
    mi << materialIndex;
    
    if(parameters->IsAvailable("Materials.Mat_" + mi.str() + ".MooneyRivlin.C10"))
        c10_ =  parameters->Get<double>("Materials.Mat_" + mi.str() + ".MooneyRivlin.C10");
    else if(parameters->IsAvailable("Materials.Mat_" + mi.str() + ".MooneyRivlin.C1"))
        c10_ =  parameters->Get<double>("Materials.Mat_" + mi.str() + ".MooneyRivlin.C1");
    else if(parameters->IsAvailable("Materials.Mat_Default.MooneyRivlin.C10"))
        c10_ =  parameters->Get<double>("Materials.Mat_Default.MooneyRivlin.C10");
    else if(parameters->IsAvailable("Materials.Mat_Default.MooneyRivlin.C1"))
        c10_ =  parameters->Get<double>("Materials.Mat_Default.MooneyRivlin.C1");
    else
        throw( std::runtime_error("CBConstitutiveModelMooneyRivlin::Init(ParameterMap* parameters, TInt materialIndex): No Mooney-Rivlin parameter C10 or C1 found for material " + mi.str() + ". This parameter has to be defined !!!"));
    
    
    if(parameters->IsAvailable("Materials.Mat_" + mi.str() + ".MooneyRivlin.C01"))
        c01_ =  parameters->Get<double>("Materials.Mat_" + mi.str() + ".MooneyRivlin.C01");
    else if(parameters->IsAvailable("Materials.Mat_" + mi.str() + ".MooneyRivlin.C2"))
        c01_ =  parameters->Get<double>("Materials.Mat_" + mi.str() + ".MooneyRivlin.C2");
    else if(parameters->IsAvailable("Materials.Mat_Default.MooneyRivlin.C01"))
        c01_ =  parameters->Get<double>("Materials.Mat_Default.MooneyRivlin.C01");
    else if(parameters->IsAvailable("Materials.Mat_Default.MooneyRivlin.C2"))
        c01_ =  parameters->Get<double>("Materials.Mat_Default.MooneyRivlin.C2");
    else
        c01_ = 0;
    
    
    if(!parameters->IsAvailable("Materials.Mat_" + mi.str() + ".MooneyRivlin.C20") && parameters->IsAvailable("Materials.Mat_Default.MooneyRivlin.C20"))
        c20_ =  parameters->Get<double>("Materials.Mat_Default.MooneyRivlin.C20");
    else
        c20_ =  parameters->Get<double>("Materials.Mat_" + mi.str() + ".MooneyRivlin.C20",0);
    
    if(!parameters->IsAvailable("Materials.Mat_" + mi.str() + ".MooneyRivlin.C02") && parameters->IsAvailable("Materials.Mat_Default.MooneyRivlin.C02"))
        c02_ =  parameters->Get<double>("Materials.Mat_Default.MooneyRivlin.C02");
    else
        c02_ =  parameters->Get<double>("Materials.Mat_" + mi.str() + ".MooneyRivlin.C02",0);
    
    if(!parameters->IsAvailable("Materials.Mat_" + mi.str() + ".MooneyRivlin.C11") && parameters->IsAvailable("Materials.Mat_Default.MooneyRivlin.C11"))
        c11_ =  parameters->Get<double>("Materials.Mat_Default.MooneyRivlin.C11");
    else
        c11_ =  parameters->Get<double>("Materials.Mat_" + mi.str() + ".MooneyRivlin.C11",0);
    
    
    if(!parameters->IsAvailable("Materials.Mat_" + mi.str() + ".MooneyRivlin.B") && parameters->IsAvailable("Materials.Mat_Default.MooneyRivlin.B"))
        b_ =  parameters->Get<double>("Materials.Mat_Default.MooneyRivlin.B");
    else
        b_ =  parameters->Get<double>("Materials.Mat_" + mi.str() + ".MooneyRivlin.B",1e5);
    
}




CBStatus CBConstitutiveModelMooneyRivlin::CalcEnergy(const Matrix3<TFloat>& deformationTensor, TFloat& energy)
{
    
    if (!Base::ignoreCorruptElements_)
        if(deformationTensor.Det() <= 0)
            return(CBStatus::CORRUPT_ELEMENT);
    
    Matrix3<TFloat> rightGreenStrain = deformationTensor.GetTranspose() * deformationTensor;
    
    TFloat invariant1   = rightGreenStrain.Invariant1();
    TFloat invariant2   = rightGreenStrain.Invariant2();
    TFloat invariant3   = rightGreenStrain.Invariant3();
    TFloat lnInvariant3 = log(invariant3);
    
    energy =  (c10_ * (invariant1 - 3))  +  (c01_ * (invariant2 - 3)) + (c20_ * (invariant1 - 3) * (invariant1 - 3)) + (c02_ * (invariant1 - 3) * (invariant2 -3)) + (c11_ * (invariant2 - 3) * (invariant2 - 3)) + (-(c10_ + 2 * c01_) * lnInvariant3) + (0.5 * b_ * lnInvariant3 * lnInvariant3);
    if(std::isnan(energy))
        return(CBStatus::NOT_A_NUMBER);
    if(std::isinf(energy))
        return(CBStatus::INFINITIVE);
    else
        
        return(CBStatus::SUCCESS);
}


CBStatus CBConstitutiveModelMooneyRivlin::CalcPK2Stress(const Matrix3<TFloat>& deformationTensor, Matrix3<TFloat>& pk2Stress)
{
    if(deformationTensor.Det() <= 0)
        if (!Base::ignoreCorruptElements_)
            return(CBStatus::CORRUPT_ELEMENT);
    
    Matrix3<TFloat> rightGreenStrain        =  deformationTensor.GetTranspose() * deformationTensor;
    Matrix3<TFloat> rightGreenStrainInverse = rightGreenStrain.GetInverse();
    
    TFloat invariant1 = rightGreenStrain.Invariant1();
    TFloat invariant2 = rightGreenStrain.Invariant2();
    
    TFloat invariant3 = rightGreenStrain.Invariant3();
    if(deformationTensor.Det() <= 0)
        invariant3 = 1;
    
    
    TFloat dWdI1 = c10_ + 2.0 * c20_ * (invariant1	- 3.0) + c11_ * (invariant2 - 3.0);
    TFloat dWdI2 = c01_ + 2.0 * c02_ * (invariant2  - 3.0) + c11_ * (invariant1 - 3.0);
    
    
    // passive deformation
    pk2Stress = 2.0  * (dWdI1 + dWdI2 * invariant1) * identity_ - 2.0 * dWdI2 * rightGreenStrain  - 2.0 * (c10_ + 2.0 * c01_) * rightGreenStrainInverse;
    
    // volume penalty
    double d = 0.8;
    Matrix3<TFloat> h;
    if(invariant3 > d)
        h = (b_ * log(invariant3)) * rightGreenStrainInverse;
    else
        h = ( b_ * log(d) +  b_/d * (invariant3-d))*rightGreenStrainInverse;
    pk2Stress += h;
    
    for(int i = 0; i < 9; i++)
    {
        if(std::isnan(pk2Stress(i)))
            return(CBStatus::NOT_A_NUMBER);
        if(std::isinf(pk2Stress(i)))
            return(CBStatus::INFINITIVE);
    }
    return(CBStatus::SUCCESS);
}

