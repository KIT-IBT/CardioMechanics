/* -------------------------------------------------------
 
 CBLoadUnloadedState.h
 
 Ver. 1.1.0
 
 Created:       Steffen Schuler (07.02.2017)
 Last modified: Tobias Gerach  (13.02.2023)
 
 Institute of Biomedical Engineering
 Karlsruhe Institute of Technology (KIT)
 
 http://www.ibt.kit.edu
 
 Copyright 2000-2009 - All rights reserved.
 
 ------------------------------------------------------ */

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
