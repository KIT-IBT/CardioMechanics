/*
 * File: CBDataPerMaterial.cpp
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


#include "CBDataPerMaterial.h"

CBDataPerMaterial::CBDataPerMaterial(ParameterMap* parameters)
{
    std::string type = parameters->Get<std::string>("ActiveStress.Type", "notset");
    if(type == "notset")
    {
        std::vector<std::string> driverKeys = parameters->GetChildNodes("ActiveStress");
        for(auto &driverKey : driverKeys)
        {
            type = parameters->Get<std::string>(driverKey + ".Type", "None");
            if(type == "None")
                continue;
            else if(type == "FromFunction")
                driver_ = std::make_shared<CBDataFromFunction>(parameters, driverKey + ".FromFunction");
            else if(type == "FromFile")
                driver_ = std::make_shared<CBDataFromFile>(parameters, driverKey + ".FromFile");
            else
                throw std::runtime_error("CBDataPerMaterial::CBDataPerMaterial: ActiveStress Type " + type + " is unkown.");
            
            std::string materialsString = parameters->Get<std::string>(driverKey + ".Materials", "All");
            if(materialsString == "All" || materialsString == "all")
            {
                allMaterials = true;
                break;
            }
            else
            {
                std::vector<TInt> materials = parameters->GetArray<TInt>(driverKey + ".Materials");
                for(auto material : materials)
                    driverMap_.insert(std::pair<TInt, std::shared_ptr<CBData>>(material, driver_));
            }
        }
    }
    else
    {
        allMaterials = true;
        
        if(type == "None")
            ;
        else if(type == "FromFunction")
            driver_ = std::make_shared<CBDataFromFunction>(parameters, "ActiveStress.FromFunction");
        else if(type == "FromFile")
            driver_ = std::make_shared<CBDataFromFile>(parameters, "ActiveStress.FromFile");
        else
            throw std::runtime_error("CBDataPerMaterial::CBDataPerMaterial: ActiveStress Type " + type + " is unkown.");
    }
}

TFloat CBDataPerMaterial::Get(TFloat time, TInt index, TInt material)
{
    if(allMaterials)
    {
        if(driver_ == 0)
            return 0.0;
        else
            return driver_->Get(time, index);
    }
    else
    {
        std::map<TInt, std::shared_ptr<CBData>>::iterator it = driverMap_.find(material);
        if(it == driverMap_.end())
            return 0.0;
        else
            return it->second->Get(time, index);
    }
}
