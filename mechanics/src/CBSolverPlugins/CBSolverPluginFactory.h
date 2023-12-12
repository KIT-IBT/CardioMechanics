/*
 * File: CBSolverPluginFactory.h
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
