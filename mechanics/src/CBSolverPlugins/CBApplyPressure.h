/*
 * File: CBApplyPressure.h
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

#include "CBSolverPlugin.h"
#include "CBElementCavity.h"

#include "Matrix3.h"
#include <petscksp.h>

#include <vector>

using namespace math_pack;

class CBApplyPressure: public CBSolverPlugin
{
public:
    CBApplyPressure(){}
    virtual ~CBApplyPressure(){}
    
    void Init() override;
    void Apply(TFloat time) override;
    void ApplyToNodalForces() override;
    void ApplyToNodalForcesJacobian() override;
    TFloat CalcVolume();
    void StepBack() override;
    void WriteToFile(TFloat time) override;
    std::string GetName() override {return "Apply Pressure";}
protected:
private:
    std::vector<CBElementCavity*> targetElements_;
    TInt surface_;
    TFloat currentPressure_=0;
    TFloat startTime_;
    TFloat stopTime_;
    TFloat currentVolume_;
    TFloat maxPressure_;
    TFloat volumeWork_=0;
    TFloat lastVolumeWork_=0;
    bool keepMaxPressure_ = true;
    typedef CBSolverPlugin   Base;
};
