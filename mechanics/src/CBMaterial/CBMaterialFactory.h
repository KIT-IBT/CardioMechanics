/*
 *  CBMaterialFactory.h
 *  CardioMechanics_local
 *
 *  Created by Thomas Fritz on 06.10.10.
 *  Copyright 2010 Institute for Biomedical Engineering, Karlsruhe Institute of Technology (KIT). All rights reserved.
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
