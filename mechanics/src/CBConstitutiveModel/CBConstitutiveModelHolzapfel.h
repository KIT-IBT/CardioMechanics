/*
 *  CBConstitutiveModelHolzapfel.h
 *  CardioMechanics
 *
 *  Created by Robin Moss, 10.18.
 *
 *
 */

#ifndef CB_CONSTITUTIVE_MODEL_CBConstitutiveModelHolzapfel
#define CB_CONSTITUTIVE_MODEL_CBConstitutiveModelHolzapfel

#include <string>
#include <sstream>

#include "Matrix3.h"

#include "ParameterMap.h"

#include "CBConstitutiveModelGuccione.h"


class CBConstitutiveModelHolzapfel : public CBConstitutiveModel {
public:
    CBConstitutiveModelHolzapfel() {}
    
    void Init(ParameterMap *parameters, TInt materialIndex);
    CBStatus CalcEnergy(const Matrix3<TFloat> &deformationTensor, TFloat &energy);
    CBStatus CalcPK2Stress(const Matrix3<TFloat> &deformationTensor, Matrix3<TFloat> &pk2Stress);
    TFloat Heavyside(TFloat I4);
    
protected:
    typedef CBConstitutiveModel Base;
    TFloat          a_, b_, af_, bf_, as_, bs_, afs_, bfs_, k_, kappa_;
    Matrix3<TFloat> identity_ = Matrix3<TFloat>::Identity();
    Matrix3<TFloat> fxf_ =      Matrix3<TFloat>(1, 0, 0, 0, 0, 0, 0, 0, 0); // 11 ff
    Matrix3<TFloat> sxs_ =      Matrix3<TFloat>(0, 0, 0, 0, 1, 0, 0, 0, 0); // 22 ss
    Matrix3<TFloat> fxs_ =      Matrix3<TFloat>(0, 1, 0, 0, 0, 0, 0, 0, 0); // 12 fs
    Matrix3<TFloat> sxf_ =      Matrix3<TFloat>(0, 0, 0, 1, 0, 0, 0, 0, 0); // 21 sf
    
private:
    void ThrowError() {
        std::runtime_error("CBConstitutiveModelHolzapfel:: functions needs ResidualDeformation");
    }
};

#endif // ifndef CB_CONSTITUTIVE_MODEL_CBConstitutiveModelHolzapfel
