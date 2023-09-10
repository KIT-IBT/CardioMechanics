/* -------------------------------------------------------
 
 CBConstitutiveModelFactory.cpp
 
 Ver. 1.1.0
 
 Created:       Thomas Fritz       (14.09.2011)
 Last modified: Tobias Gerach      (09.02.2021)
 
 Institute of Biomedical Engineering
 Karlsruhe Institute of Technology (KIT)
 
 http://www.ibt.kit.edu
 
 Copyright 2000-2009 - All rights reserved.
 
 ------------------------------------------------------ */

#ifndef CB_CONSTITUTIVE_MODEL_FACTORY
#define CB_CONSTITUTIVE_MODEL_FACTORY

#include <iostream>
#include <sstream>
#include <vector>

#include "CBConstitutiveModelGuccione.h"
#include "CBConstitutiveModelMooneyRivlin.h"
#include "CBConstitutiveModelHolzapfel.h"
#include "CBConstitutiveModelUsyk.h"
#include "CBConstitutiveModelNeoHooke.h"

class CBConstitutiveModelFactory {
public:
    ~CBConstitutiveModelFactory() {
        for (auto i : constitutiveModels_)
            delete i;
    }
    
    CBConstitutiveModel *New(ParameterMap *parameters, TInt materialIndex);
    
    bool Delete(CBConstitutiveModel *model) {
        std::set<CBConstitutiveModel *>::iterator it = constitutiveModels_.find(model);
        
        if (it != constitutiveModels_.end()) {
            delete *it;
            constitutiveModels_.erase(it);
            return true;
        } else {
            return false;
        }
    }
    
protected:
private:
    std::set<CBConstitutiveModel *> constitutiveModels_;
};

#endif  // ifndef CB_CONSTITUTIVE_MODEL_FACTORY
