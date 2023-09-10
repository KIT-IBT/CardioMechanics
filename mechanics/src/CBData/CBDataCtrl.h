/*
 *  CBDataCtrl.h
 *  CardioMechanics
 *
 *  Created by Thomas Fritz on 01.11.10.
 *  Copyright 2010 Institute for Biomedical Engineering, Karlsruhe Institute of Technology (KIT). All rights reserved.
 *
 */

#ifndef CB_DATA_CTRL_H
#define CB_DATA_CTRL_H

#include <stdint.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <cmath>

#include "ParameterMap.h"

#include "CBData.h"
#include "DCType.h"

class CBData;



class CBDataCtrl : public CBData
{
public:
    CBDataCtrl() {}
    void Init(TInt size);
    virtual TFloat Get(TFloat time, TInt index){return val_[index];}
    virtual void Set(TFloat time, TInt index,TFloat val){val_[index] = val;}
    virtual void Backup(){lastVal_ = val_;}
    virtual void LoadFromBackup(){val_ = lastVal_;}
protected:
private:
    std::vector<TFloat> val_;
    std::vector<TFloat> lastVal_;
};

#endif
