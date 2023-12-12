/*
 * File: CBConstitutiveModelMooneyRivlin.h
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
