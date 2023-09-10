/*
 *  CBTensionFactory.h
 *  CardioMechanics
 *
 *  Created by Lukas Baron on Thu Apr 27 2017.
 *  Copyright 2017 Institute for Biomedical Engineering, Karlsruhe Institute of Technology (KIT). All rights reserved.
 *
 */

#include "CBTensionFactory.h"

#include "CBTensionModelLumens.h"
#include "CBTensionFromFunction.h"
#include "CBTensionModelLand17.h"
#include "CBTensionModelTanH.h"
#include "CBTensionModelBestel.h"

#include "CBTiming.h"
#include "CBFileManager.h"
#include "DCCtrl.h"

CBTensionFactory::CBTensionFactory() {}

void CBTensionFactory::Init(ParameterMap *parameters, CBTiming *timing, CBFileManager *fileManager) {
    // tensionFactory has a filemanager object to track the files depending on
    // the elements and materials
    fileManager_   = fileManager;
    parameters_    = parameters;
    timing_        = timing;
    isInitialized_ = true;
    
    producers_["None"]          = [](CBElementSolid *ele) {return new CBNoTension(); };
    producers_["File"]          = [this](CBElementSolid *ele) {return new CBFileTension(fileManager_, ele); };
    producers_["Lumens"]        = [this](CBElementSolid *ele) {return new CBLumensTension(ele, parameters_); };
    producers_["FromFunction"]  = [this](CBElementSolid *ele) {return new CBTensionFromFunction(parameters_, ele); };
    producers_["Land17"]        = [this](CBElementSolid *ele) {return new CBTensionModelLand17(ele, parameters_); };
    producers_["TanH"]          = [this](CBElementSolid *ele) {return new CBTensionModelTanH(ele, parameters_); };
    producers_["Bestel"]        = [this](CBElementSolid *ele) {return new CBTensionModelBestel(ele, parameters_); };
}

/// returns a new tensionModel, that fits to the tensionName found in the passed elements material properties
CBTensionModel *CBTensionFactory::New(CBElementSolid *ele) {
    CBTensionModel *tensionModel = 0;
    
    std::string elementType = ele->GetMaterial()->GetTensionName();
    auto producer = producers_.find(elementType);
    
    if (producer == producers_.end()) {
        DCCtrl::print << "Material " << ele->GetMaterial()->GetMaterialIndex() << " has strange tension type." << std::endl;
        throw std::runtime_error(std::string("Unkown tension type: [" + elementType +
                                             "] You might have to extend CBTensionFactory \nAvailable tension models are: \n")
                                 + "\t None\t: No tension (default)\n"
                                 + "\t File\t: read tensions from input files listed in a tensionlist\n"
                                 + "\t Lumens\t: ode-based tension model\n"
                                 + "\t TanH\t: length-dependent analytic function\n"
                                 + "\t Bestel\t, Bestel2001\t: ode-based, time-dependent stress function\n"
                                 + "\t Land17\t: ode-based tension model designed for humans\n"
                                 + "\t FromFunction\t: calculate tension through functions \n"
                                 = ")");
    } else {
        ProducerFunction pf = producer->second;
        tensionModel = pf(ele);
    }
    assert(tensionModel != 0);
    return tensionModel;
} // CBTensionFactory::New
