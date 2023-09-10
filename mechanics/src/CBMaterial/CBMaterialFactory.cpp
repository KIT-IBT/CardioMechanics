/*
 *  CBMaterialFactory.cpp
 *  CardioMechanics
 *
 *  Created by Thomas Fritz on 14.09.11.
 *  Copyright 2011 Institute for Biomedical Engineering, Karlsruhe Institute of Technology (KIT). All rights reserved.
 *
 */

#include "CBMaterialFactory.h"



CBMaterial* CBMaterialFactory::New(ParameterMap* parameters, TInt materialIndex)
{
    CBMaterial*          material = new CBMaterial;
    material->Init(parameters, materialIndex);
    CBConstitutiveModel* model = constitutiveModelFactory_.New(parameters, materialIndex);
    model->Init(parameters, materialIndex);
    material->SetConstitutiveModel(model);
    materials_.insert(material);
    return(material);
}

