/*
 * File: CBLoadUnloadedState.h
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

#ifndef CB_LOAD_UNLOADED_STATE
# define CB_LOAD_UNLOADED_STATE

# include <vector>

# include "CBSolverPlugin.h"

using namespace math_pack;

class CBLoadUnloadedState : public CBSolverPlugin {
public:
    CBLoadUnloadedState() {}
    
    virtual ~CBLoadUnloadedState() {}
    
    void Init() override;
    void Apply(TFloat time) override;
    
    CBStatus GetStatus() override { return status_; }
    
    std::string GetName() override { return "LoadUnloadedState"; }
    
private:
    void InitPetscVectors();
    void DeInitPetscVectors();
    void LoadNodes(std::string filename);
    
    bool isFirstStep_ = true;
    bool settleDown_;
    
    Vec originalCoords_;
    Vec coordsSeq_;
    Vec unloadedCoords_;
    Vec inflatedCoords_;
    VecScatter coordsScatter_;
};

#endif // ifndef CB_LOAD_UNLOADED_STATE
