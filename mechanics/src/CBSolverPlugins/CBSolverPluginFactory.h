/*
 *  CBSolverPluginFactory.h
 *  CardioMechanics
 *
 *  Created by Thomas Fritz on 19.04.11.
 *  Copyright 2011 Institute for Biomedical Engineering, Karlsruhe Institute of Technology (KIT). All rights reserved.
 *
 */

#ifndef CB_SOLVER_PLUGIN_FACTORY
#define CB_SOLVER_PLUGIN_FACTORY

#include <vector>
#include <set>
#include "ParameterMap.h"

class CBSolverPlugin;


class CBSolverPluginFactory
{
public:
    
    CBSolverPluginFactory(){}
    
    ~CBSolverPluginFactory()
    {
        for( std::vector<CBSolverPlugin* >::iterator it = plugins_.begin(); it != plugins_.end(); it++)
            delete(*it);
    }
    CBSolverPlugin* New(std::string pluginName);
    
    void LoadAllPlugins(std::vector<CBSolverPlugin* >& plugins, ParameterMap* parameters);
    
protected:
private:
    std::vector<CBSolverPlugin* > plugins_;
};

#endif
