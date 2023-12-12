/*
 * File: CBConstitutiveModelFactory.cpp
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



#include <string>

#include "CBConstitutiveModelFactory.h"

CBConstitutiveModel *CBConstitutiveModelFactory::New(ParameterMap *parameters, TInt materialIndex) {
    std::string modelType;
    
    if ((parameters->IsAvailable("Materials.Mat_Default.Type") == true) &&
        (parameters->IsAvailable("Materials.Mat_" + std::to_string(materialIndex) + ".Type") == false) )
        modelType = parameters->Get<std::string>("Materials.Mat_Default.Type");
    else
        modelType = parameters->Get<std::string>("Materials.Mat_" + std::to_string(materialIndex) + ".Type");
    
    CBConstitutiveModel *model;
    
    if (modelType == std::string("Guccione"))
        model = new CBConstitutiveModelGuccione;
    else if (modelType == std::string("MooneyRivlin"))
        model = new CBConstitutiveModelMooneyRivlin;
    else if (modelType == std::string("Holzapfel"))
        model = new CBConstitutiveModelHolzapfel;
    else if (modelType == std::string("Usyk"))
        model = new CBConstitutiveModelUsyk;
    else if (modelType == std::string("NeoHooke"))
        model = new CBConstitutiveModelNeoHooke;
    
    if (!model) {
        throw std::runtime_error(
                                 "Unkown constitutive model type: " + modelType + " You might have to extend CBConstitutiveModelFactory\n");
    } else {
        model->SetModelIndex(constitutiveModels_.size());
        constitutiveModels_.insert(model);
        return model;
    }
}  // CBConstitutiveModelFactory::New
