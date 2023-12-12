/*
 * File: CBSolverPlugin.h
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


#ifndef CB_SOLVER_PLUGIN
#define CB_SOLVER_PLUGIN

#include "Matrix3.h"


#include "ParameterMap.h"

#include "CBElementAdapter.h"
#include "CBStatus.h"


using namespace math_pack;

class CBElementAdapter;
class ParameterMap;
class CBSolver;

class CBSolverPlugin
{ 
public:
    CBSolverPlugin(){}
    virtual ~CBSolverPlugin(){}
    // ---- Pure Virtual Functions ----
    virtual void Apply(TFloat time) = 0;
    virtual std::string GetName()=0;
    virtual void Init()=0;
    // --------------------------------
    
    virtual void AnalyzeResults(){}
    virtual bool WantsToAnalyzeResults(){return false;}
    virtual void StepBack(){}
    virtual void ApplyToNodalForces(){}
    virtual void ApplyToNodalForcesJacobian(){}
    virtual bool ExitCheck(TFloat time) { return false; }
    virtual void WriteToFile(TFloat time){}
    virtual void Export(TFloat time){}
    virtual void ExportState(ParameterMap* PluginParameters){}
    virtual void Prepare(){status_ = CBStatus::DACCORD;}
    virtual double GetPreparationProgress(){return -1.0;}
    virtual void Reset(){}
    virtual void SetAdapter(CBElementAdapter* adapter){adapter_ = adapter; }
    virtual CBStatus GetStatus(){return status_;}
    void SetParameters(ParameterMap* parameters){parameters_ = parameters; }
    void SwitchOff(){isActive_ = false;}
    void SwitchOn(){isActive_ = true;}
    bool IsActive(){return isActive_;}
    ParameterMap* GetParameters(){return(parameters_); }
    CBElementAdapter* GetAdapter(){return(adapter_); }
    CBSolver* GetSolver(){return adapter_->GetSolver();}
    
    virtual void LoadState(ParameterMap* pluginParameters){}
    
protected:
    
    CBElementAdapter* adapter_;
    ParameterMap* parameters_;
    CBStatus status_ = CBStatus::WAITING;
    
private:
    bool isActive_ = true;
};

#endif
