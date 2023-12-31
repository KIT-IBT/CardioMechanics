/*
 * File: CBDataFromFile.h
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


#ifndef CB_DATA_FROM_FILE_H
#define CB_DATA_FROM_FILE_H

#include <stdint.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>

#include "ParameterMap.h"

#include "CBData.h"
#include "DCType.h"

class CBData;


class CBDataFromFile : public CBData
{
public:
    CBDataFromFile();
    CBDataFromFile(ParameterMap* parameters, std::string parameterKey);
    void Init(std::string filename);
    virtual TFloat Get(TFloat time, TInt index);
    std::vector<TFloat> GetTimePoints();
    
protected:
    
private:
    virtual void LoadFileList(std::string filename);
    bool LoadDataSet(TFloat time);
    
    TFloat startTime_;
    TFloat period_;
    
    bool isDataSetLoaded_;
    TFloat                                      timeDataBegin_ = 0;
    TFloat                                      timeDataEnd_ = 0;
    std::vector<TFloat>                         dataBegin_;
    std::vector<TFloat>                         dataEnd_;
    std::vector<std::pair<TFloat, std::string>> fileList_;
};


#endif
