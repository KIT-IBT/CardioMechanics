/*
 * File: CBDataCtrl.h
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
