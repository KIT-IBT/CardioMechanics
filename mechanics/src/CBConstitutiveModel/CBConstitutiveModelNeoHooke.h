/* -------------------------------------------------------
 
 CBConstitutiveModelNeoHooke.h
 
 Ver. 1.0.0
 
 Created:       Tobias Gerach      (09.02.2021)
 Last modified: Tobias Gerach      (09.02.2021)
 
 Institute of Biomedical Engineering
 Karlsruhe Institute of Technology (KIT)
 
 http://www.ibt.kit.edu
 
 Copyright 2000-2009 - All rights reserved.
 
 ------------------------------------------------------ */

#ifndef CB_CONSTITUTIVE_MODEL_NEOHOOKE
#define CB_CONSTITUTIVE_MODEL_NEOHOOKE

#include <string>
#include <sstream>

#include "Matrix3.h"

#include "ParameterMap.h"

#include "CBConstitutiveModel.h"


class CBConstitutiveModelNeoHooke : public CBConstitutiveModel {
public:
    CBConstitutiveModelNeoHooke() : a_(0), k_(0) {
        identity_.SetToIdentityMatrix();
    }
    
    CBStatus CalcEnergy(const Matrix3<TFloat> &deformationTensor, TFloat &energy);
    CBStatus CalcPK2Stress(const Matrix3<TFloat> &deformationTensor, Matrix3<TFloat> &pk2Stress);
    void     Init(ParameterMap *parameters, TInt materialIndex);
    
protected:
private:
    typedef CBConstitutiveModel Base;
    TFloat a_, k_;
    Matrix3<TFloat> identity_;
};

#endif  // ifndef CB_CONSTITUTIVE_MODEL_NEOHOOKE
