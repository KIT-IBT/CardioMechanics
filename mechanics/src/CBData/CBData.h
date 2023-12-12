/*
 * File: CBData.h
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


#ifndef CB_DATA_H
#define CB_DATA_H

#include <iostream>
#include "ParameterMap.h"

#include "DCType.h"

class CBData
{
public:
    CBData() : time_(0) {}
    virtual ~CBData(){}
    virtual void Init(){}
    virtual void Init(TInt size){}
    virtual void Init(std::string filename){}
    virtual TFloat Get(TFloat time, TInt index) = 0;
    virtual TFloat Get(TFloat time, TInt index, TInt material) { std::cout << "warning, this get function should not be called at all! It returns a value of zero, which might be unexpected." << std::endl; return 0.0; }
    virtual void Set(TFloat time, TInt index, TFloat val) {}
    void SetTime(TFloat time) { time_ = time; }
    
    virtual TFloat Get(TInt index, TInt material) { return Get(time_, index, material); }
    
protected:
private:
    TFloat time_;
};

#endif
