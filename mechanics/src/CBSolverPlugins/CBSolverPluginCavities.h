/*
 * File: CBSolverPluginCavities.h
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


#pragma once

#include <map>
#include <memory>

#include "CBSolverPlugin.h"
#include "CBCirculationCavity.h"

/// Virtual Base class providing cavities for each surface number. Cavities can get pressure and support volume computation.
class CBSolverPluginCavities : public CBSolverPlugin
{
public:
    CBSolverPluginCavities() : CBSolverPlugin() {}
    virtual ~CBSolverPluginCavities(){}
    
    virtual void Init() = 0;
    virtual void Apply(TFloat time)  = 0;
    virtual std::string GetName() override { return("SolverPluginCavities"); };
    
    virtual void InitCavities(std::string pluginRootKey);
    virtual void CreateCavities();
    virtual void SyncCavities();
    
    // std::map<TInt, CBCirculationCavity*>& GetCavities();
    // CBCirculationCavity* GetCavity(int cavityIndex);
    
protected:
    std::map<TInt, CBCirculationCavity*> cavities_; // maps surface indices with CBCirculationCavities
private:
    std::vector<TInt> ignoredSurfaceIndices_;
    std::vector<TInt> ignoredSurfaceElementsIndices_;
    std::vector<TInt> pokedNodeIndices_;
    
    typedef CBSolverPlugin Base;
};
