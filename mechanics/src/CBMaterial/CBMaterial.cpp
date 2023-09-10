/*
 *  CBMaterial.cpp
 *  CardioMechanics
 *
 *  Created by Thomas Fritz on 14.09.11.
 *  Copyright 2011 Institute for Biomedical Engineering, Karlsruhe Institute of Technology (KIT). All rights reserved.
 *
 */

#include "CBMaterial.h"



void CBMaterial::Init(ParameterMap* parameters, TInt materialIndex)
{
    materialIndex_ = materialIndex;
    properties_ = new CBMaterial::Properties; 
    
    if(parameters->IsAvailable("Materials.Mat_" + std::to_string(materialIndex_) + ".Density") == false && parameters->IsAvailable("Materials.Mat_Default.Density") == true)
        properties_->massDensity_ = parameters->Get<double>("Materials.Mat_Default.Density");
    else
        properties_->massDensity_ = parameters->Get<double>("Materials.Mat_" + std::to_string(materialIndex_) + ".Density", 0);
    
    if(parameters->IsAvailable("Materials.Mat_" + std::to_string(materialIndex_) + ".Rayleigh.Alpha") == false && parameters->IsAvailable("Materials.Mat_Default.Rayleigh.Alpha") == true)
        properties_->rayleighAlpha_ = parameters->Get<double>("Materials.Mat_Default.Rayleigh.Alpha");
    else
        properties_->rayleighAlpha_ = parameters->Get<double>("Materials.Mat_" + std::to_string(materialIndex_) + ".Rayleigh.Alpha", 0);
    
    if(parameters->IsAvailable("Materials.Mat_" + std::to_string(materialIndex_) + ".Rayleigh.Beta") == false && parameters->IsAvailable("Materials.Mat_Default.Rayleigh.Beta") == true)
        properties_->rayleighBeta_ = parameters->Get<double>("Materials.Mat_Default.Rayleigh.Beta");
    else
        properties_->rayleighBeta_ = parameters->Get<double>("Materials.Mat_" + std::to_string(materialIndex_) + ".Rayleigh.Beta", 0);
    
    if(parameters->IsAvailable("Materials.Mat_" + std::to_string(materialIndex_) + ".TensionMax") == false && parameters->IsAvailable("Materials.Mat_Default.TensionMax") == true)
        properties_->tensionMax_ = parameters->Get<double>("Materials.Mat_Default.TensionMax");
    else
        properties_->tensionMax_ = parameters->Get<double>("Materials.Mat_" + std::to_string(materialIndex_) + ".TensionMax", 0);
    
    if(parameters->IsAvailable("Materials.Mat_" + std::to_string(materialIndex_) + ".TensionModel") == false && parameters->IsAvailable("Materials.Mat_Default.TensionModel") == true)
        properties_->tensionName_ = parameters->Get<std::string>("Materials.Mat_Default.TensionModel");
    else
        properties_->tensionName_ = parameters->Get<std::string>("Materials.Mat_" + std::to_string(materialIndex_) + ".TensionModel", "None");
    
    if(parameters->IsAvailable("Materials.Mat_" + std::to_string(materialIndex_) + ".SimpleElectroMechanicalFeedback") == false && parameters->IsAvailable("Materials.Mat_Default.SimpleElectroMechanicalFeedback") == true)
        properties_->useSimpleElectroMechanicalFeedback_ = parameters->Get<bool>("Materials.Mat_Default.SimpleElectroMechanicalFeedback");
    else
        properties_->useSimpleElectroMechanicalFeedback_ = parameters->Get<bool>("Materials.Mat_" + std::to_string(materialIndex_) + ".SimpleElectroMechanicalFeedback", 0);
    
    if(parameters->IsAvailable("Materials.Mat_" + std::to_string(materialIndex_) + ".SimpleElectroMechanicalFeedbackType") == false && parameters->IsAvailable("Materials.Mat_Default.SimpleElectroMechanicalFeedbackType") == true)
        properties_->simpleElectroMechanicalFeedbackType_ = parameters->Get<std::string>("Materials.Mat_Default.SimpleElectroMechanicalFeedbackType", "Linear");
    else
        properties_->simpleElectroMechanicalFeedbackType_ = parameters->Get<std::string>("Materials.Mat_" + std::to_string(materialIndex_) + ".SimpleElectroMechanicalFeedbackType", "Linear");
    
    if(parameters->IsAvailable("Materials.Mat_" + std::to_string(materialIndex_) + ".Conductivity") == false && parameters->IsAvailable("Materials.Mat_Default.Conductivity") == true)
    {
        std::vector<TFloat> c = parameters->GetArray<double>("Materials.Mat_Default.Conductivity");
        
        if(c.size() != 3)
            throw std::runtime_error("void CBMaterial::Init(ParameterMap* parameters, TInt materialIndex): Materials.Mat_.Conductivity");
        
        properties_->conductivity_.SetToZero();
        properties_->conductivity_(0,0) = c[0];
        properties_->conductivity_(1,1) = c[1];
        properties_->conductivity_(2,2) = c[2];
        
    }
    else
        properties_->conductivity_.SetToIdentityMatrix();
}
