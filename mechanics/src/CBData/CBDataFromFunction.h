/*
 *  CBDataFromFunction.h
 *  CardioMechanics
 *
 *  Created by Thomas Fritz on 01.11.10.
 *  Copyright 2010 Institute for Biomedical Engineering, Karlsruhe Institute of Technology (KIT). All rights reserved.
 *
 */

#ifndef CB_DATA_FROM_FUNCTION_H
#define CB_DATA_FROM_FUNCTION_H

#include "ParameterMap.h"

#include "CBData.h"
#include "DCType.h"

class CBData;

class CBAnalyticFunction;

class CBDataFromFunction : public CBData
{
public:
    CBDataFromFunction(ParameterMap* parameters, std::string parameterKey);
    ~CBDataFromFunction();
    virtual TFloat Get(TFloat time, TInt index);
    
protected:
    
private:
    TFloat offset_ = 0;
    CBAnalyticFunction* getImpl_ = 0;
};
#endif
