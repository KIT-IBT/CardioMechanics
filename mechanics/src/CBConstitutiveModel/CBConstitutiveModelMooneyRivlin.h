
/*
 *  CBConstitutiveModelMooneyRivlin.h
 *  CardioMechanics
 *
 *  Created by Thomas Fritz on 29.09.10.
 *  Copyright 2010 Institute for Biomedical Engineering, Karlsruhe Institute of Technology (KIT). All rights reserved.
 *
 */


#ifndef CB_CONSTITUTIVE_MODEL_MOONEY_RIVLIN
#define CB_CONSTITUTIVE_MODEL_MOONEY_RIVLIN

#include <string>
#include <sstream>
#include <math.h>

#include "Matrix3.h"
#include "ParameterMap.h"
#include "CBConstitutiveModel.h"



class CBConstitutiveModelMooneyRivlin : public CBConstitutiveModel
{
public:
    CBConstitutiveModelMooneyRivlin()
    {
        identity_.SetToIdentityMatrix();
    }
    
    void Init(ParameterMap* parameters, TInt materialIndex);
    
    CBStatus CalcEnergy(const Matrix3<TFloat>& deformationTensor, TFloat& energy);
    CBStatus CalcPK2Stress(const Matrix3<TFloat>& deformationTensor, Matrix3<TFloat>& pk2Stress);
protected:
private:
    typedef CBConstitutiveModel Base;
    
    TFloat          c10_ = 0;
    TFloat          c20_ = 0;
    TFloat          c01_ = 0;
    TFloat          c02_ = 0;
    TFloat          c11_ = 0;
    TFloat            b_ = 0;
    
    Matrix3<TFloat> identity_;
};

#endif
