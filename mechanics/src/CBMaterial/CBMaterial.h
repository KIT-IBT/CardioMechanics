/*
 *  CBMaterial.h
 *  CardioMechanics_local
 *
 *  Created by Thomas Fritz on 06.10.10.
 *  Copyright 2010 Institute for Biomedical Engineering, Karlsruhe Institute of Technology (KIT). All rights reserved.
 *
 */

#ifndef CB_MATERIAL_H
#define CB_MATERIAL_H

#include "Matrix3.h"
#include "ParameterMap.h"

#include "CBStatus.h"
//#include "CBStressModel.h"
#include "CBConstitutiveModel.h"
#include "DCType.h"

using namespace math_pack;
class CBMaterial
{
private:
    class Properties
    {
    public:
        //        CBStressModel* stressModel_=0;
        CBConstitutiveModel* constitutiveModel_=0;
        TInt initAsMaterial_=-1;
        TFloat massDensity_=0;
        TFloat rayleighAlpha_=0;
        TFloat rayleighBeta_=0;
        TFloat tensionMax_=0;
        std::string tensionName_ = "not defined yet";
        bool useSimpleElectroMechanicalFeedback_ = false;
        std::string simpleElectroMechanicalFeedbackType_;
        Matrix3<TFloat> conductivity_;
    };
public:
    ~CBMaterial(){delete properties_;}
    void Init(ParameterMap* parameters, TInt materialIndex);
    TInt GetMaterialIndex(){return(materialIndex_); }
    TInt GetInitAsMaterial(){return properties_->initAsMaterial_;}
    void SetMaterialIndex(TInt materialIndex){materialIndex_ =  materialIndex; }
    //    void SetStressModel(CBStressModel* stressModel){properties_->stressModel_ = stressModel; }
    //    inline CBStressModel* GetStressModel(){return(properties_->stressModel_); }
    void SetTensionName(std::string tensionName) { properties_->tensionName_ = tensionName; }
    std::string GetTensionName() { return properties_->tensionName_; }
    void SetConstitutiveModel(CBConstitutiveModel* constitutiveModel){properties_->constitutiveModel_ = constitutiveModel; }
    inline CBConstitutiveModel* GetConstitutiveModel(){return(properties_->constitutiveModel_); }
    void SetMassDensity(TFloat massDensity){properties_->massDensity_ = massDensity; }
    void SetRayleighAlpha(TFloat rayleighAlpha){properties_->rayleighAlpha_ = rayleighAlpha; }
    void SetRayleighBeta(TFloat rayleighBeta){properties_->rayleighBeta_ = rayleighBeta; }
    
    void BehaveLikeMaterial(CBMaterial* mat){originalProperties_ = properties_; properties_ = mat->GetProperties();}
    void UseOriginalProperites(){properties_ = originalProperties_;}
    
    Properties* GetProperties(){return properties_;}
    TFloat GetMassDensity(){return(properties_->massDensity_); }
    TFloat GetRayleighAlpha(){return(properties_->rayleighAlpha_); }
    TFloat GetRayleighBeta(){return(properties_->rayleighBeta_); }
    TFloat GetTensionMax(){return(properties_->tensionMax_); }
    Matrix3<TFloat> GetConductivity(){return(properties_->conductivity_);}
protected:
    TInt                               materialIndex_;
    Properties* properties_;
    Properties* originalProperties_;
};

#endif
