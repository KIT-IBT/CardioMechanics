/*
 * File: CBMaterialFactory.h
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


#ifndef CB_MATERIAL_FACTORY_H
#define CB_MATERIAL_FACTORY_H

#include "CBMaterial.h"
#include "CBConstitutiveModelFactory.h"
#include "ParameterMap.h"
#include <vector>



class CBMaterialFactory
{
public:
    ~CBMaterialFactory()
    {
        for(auto i:materials_)
            delete i;
    }
    
    CBMaterial* New(ParameterMap* parameters, TInt materialIndex);
    bool Delete(CBMaterial* material)
    {
        std::set<CBMaterial* >::iterator it = materials_.find(material);
        if(it != materials_.end())
        {
            delete *it;
            materials_.erase(it);
            return(true);
        }
        else
        {
            return(false);
        }
    }
protected:
private:
    std::set<CBMaterial* >     materials_;
    CBConstitutiveModelFactory constitutiveModelFactory_;
};

#endif
