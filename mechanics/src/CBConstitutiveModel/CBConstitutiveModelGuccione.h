/*
 *  CBConstitutiveModelGuccione.h
 *  CardioMechanics
 *
 *  Created by Thomas Fritz on 08.06.10.
 *  Copyright 2010 Institute for Biomedical Engineering, Karlsruhe Institute of Technology (KIT). All rights reserved.
 *
 */

#ifndef CB_CONSTITUTIVE_MODEL_GUCCIONE
#define CB_CONSTITUTIVE_MODEL_GUCCIONE

#include <string>
#include <sstream>

#include "Matrix3.h"

#include "ParameterMap.h"

#include "CBConstitutiveModel.h"


class CBConstitutiveModelGuccione : public CBConstitutiveModel {
public:
    CBConstitutiveModelGuccione() : C_(0), bf_(0), bt_(0), bfs_(0), k_(0) {identity_.SetToIdentityMatrix(); }
    
    CBStatus CalcEnergy(const Matrix3<TFloat> &deformationTensor, TFloat &energy);
    CBStatus CalcPK2Stress(const Matrix3<TFloat> &deformationTensor, Matrix3<TFloat> &pk2Stress);
    void Init(ParameterMap *parameters, TInt materialIndex);
    
protected:
    typedef CBConstitutiveModel Base;
    TFloat          C_, bf_, bt_, bfs_, k_, aScale_, bScale_;
    Matrix3<TFloat> identity_;
    
private:
};

#endif // ifndef CB_CONSTITUTIVE_MODEL_GUCCIONE
