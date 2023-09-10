/*
 *  CBDataPerMaterial.h
 *  CardioMechanics
 *
 *  Created by Lukas Baron and Steffen Schuler on 02.02.2016.
 *  Copyright 2016 Institute for Biomedical Engineering, Karlsruhe Institute of Technology (KIT). All rights reserved.
 *
 */

#ifndef CB_DATA_PER_MATERIAL_H
#define CB_DATA_PER_MATERIAL_H

#include <stdint.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <cmath>
#include <memory>

#include "ParameterMap.h"

#include "CBData.h"
#include "CBDataFromFunction.h"
#include "CBDataFromFile.h"
#include "DCType.h"

class CBDataPerMaterial : public CBData
{
public:
    CBDataPerMaterial(ParameterMap* parameters);
    ~CBDataPerMaterial() {}
    
    TFloat Get(TFloat time, TInt index) { return 0.0; }
    TFloat Get(TFloat time, TInt index, TInt material);
    
protected:
    
private:
    bool allMaterials = false;
    std::shared_ptr<CBData> driver_ = 0;
    std::map<TInt, std::shared_ptr<CBData> > driverMap_;
};
#endif
