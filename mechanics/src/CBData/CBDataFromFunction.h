/*
 * File: CBDataFromFunction.h
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
