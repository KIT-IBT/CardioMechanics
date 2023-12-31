/*
 * File: CBConstitutiveModelUsyk.h
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


#ifndef CB_CONSTITUTIVE_MODEL_USYK
#define CB_CONSTITUTIVE_MODEL_USYK

#include <string>
#include <sstream>

#include "Matrix3.h"

#include "ParameterMap.h"

#include "CBConstitutiveModel.h"


class CBConstitutiveModelUsyk : public CBConstitutiveModel {
public:
    CBConstitutiveModelUsyk() : a_(0), bff_(0), bss_(0), bnn_(0), bfs_(0), bfn_(0), bns_(0), k_(0) {
        identity_.SetToIdentityMatrix();
    }
    
    CBStatus CalcEnergy(const Matrix3<TFloat> &deformationTensor, TFloat &energy);
    CBStatus CalcPK2Stress(const Matrix3<TFloat> &deformationTensor, Matrix3<TFloat> &pk2Stress);
    void     Init(ParameterMap *parameters, TInt materialIndex);
    
protected:
private:
    typedef CBConstitutiveModel Base;
    TFloat a_, bff_, bss_, bnn_, bfs_, bfn_, bns_, k_, aScale_, bScale_;
    Matrix3<TFloat> identity_;
};

#endif  // ifndef CB_CONSTITUTIVE_MODEL_USYK
