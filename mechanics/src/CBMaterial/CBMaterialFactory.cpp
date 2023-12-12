/*
 * File: CBMaterialFactory.cpp
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

