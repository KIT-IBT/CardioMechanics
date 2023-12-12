/*
 * File: CBConstitutiveModelFactory.h
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
