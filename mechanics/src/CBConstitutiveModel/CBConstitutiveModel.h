/*
 * File: CBConstitutiveModel.h
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


#ifndef CB_CONSTITUTIVE_MODEL
#define CB_CONSTITUTIVE_MODEL

#include "Matrix3.h"
#include "ParameterMap.h"
#include "CBStatus.h"
#include "DCType.h"

using namespace math_pack;
double relaxedExp(double a, double criticalValue = 1e20);
class CBConstitutiveModel {
public:
    CBConstitutiveModel() : ignoreCorruptElements_(false) {}
    
    virtual ~CBConstitutiveModel() {}
    
    virtual void Init(ParameterMap *parameters, TInt materialIndex)                                            = 0;
    virtual CBStatus CalcEnergy(const Matrix3<TFloat> &deformationTensor, TFloat &energy)                      = 0;
    
    virtual CBStatus CalcEnergy(const Matrix3<TFloat> &deformationTensor, const Matrix3<TFloat> &ResidualDeformation,
                                TFloat &energy) {
        std::runtime_error("This material law does not provide a function for residual deformation");
        return CBStatus::INTERRUPTED;
    }
    
    virtual CBStatus CalcPK2Stress(const Matrix3<TFloat> &deformationTensor, const Matrix3<TFloat> &ResidualDeformation,
                                   Matrix3<TFloat> &pk2Stress) {
        std::runtime_error("This material law does not provide a function for residual deformation");
        return CBStatus::INTERRUPTED;
    }
    
    virtual CBStatus CalcPK2Stress(const Matrix3<TFloat> &deformationTensor, Matrix3<TFloat> &pk2Stress)       = 0;
    
    TInt GetModelIndex() {return modelIndex_; }
    
    void SetModelIndex(TInt modelIndex) {modelIndex_ = modelIndex; }
    
    bool ShouldIgnoreCorruptElements() {return ignoreCorruptElements_;}
    
protected:
    TInt modelIndex_;
    bool ignoreCorruptElements_;
    TFloat criticalVolumeChange_;
    
    TFloat Tmax_;
    
private:
};


#endif // ifndef CB_CONSTITUTIVE_MODEL
