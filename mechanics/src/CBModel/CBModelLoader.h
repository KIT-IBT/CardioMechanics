/*
 * File: CBModelLoader.h
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

#ifndef CB_MODEL_LOADER
#define CB_MODEL_LOADER

#include<string>
#include "CBModel.h"



class CBModelLoader
{
public:
    CBModelLoader(ParameterMap* parameters){InitParameters(parameters); model_ = new CBModel;}
    virtual ~CBModelLoader(){}
    virtual CBModel* Load(std::string Tag = "Mesh") = 0;
    
    void InitParameters(ParameterMap* parameters){parameters_ = parameters;}
    void TransformT4toT10(){}
protected:
    void InitModelOutputFile();
    CBModel* model_;
    ParameterMap*          parameters_;
private:
    CBModelLoader();
};


#endif
