/*
 *  CBModelLoader.h
 *  CardioMechanics
 *
 *  Created by Thomas Fritz on 13.04.10.
 *  Copyright 2010 Institute for Biomedical Engineering, Karlsruhe Institute of Technology (KIT). All rights reserved.
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
