/*
 * File: CBTiming.cpp
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


#include "CBTiming.h"

#include "DCCtrl.h"



CBTiming::CBTiming() {}

CBTiming::~CBTiming() {}

void CBTiming::Init(ParameterMap &parameter) {
    startTime_   = parameter.Get<TFloat>("Solver.StartTime", 0);
    stopTime_    = parameter.Get<TFloat>("Solver.StopTime", 1);
    
    firstTimeStep_ = parameter.Get<TFloat>("Solver.TimeStep", 1e-4);
    minTimeStep_ = parameter.Get<TFloat>("Solver.MinTimeStep", 1e-9);
    maxSimTimeStep_  = parameter.Get<TFloat>("Solver.MaxTimeStep", firstTimeStep_);
    maxInitTimeStep_ = parameter.Get<TFloat>("Solver.MaxInitTimeStep", maxSimTimeStep_);
    
    fastRelaxation_ = parameter.Get<bool>("Solver.FastRelaxation", false);
    minNumCnt_      = parameter.Get<int>("Solver.MinSteps", 4);
    
    // reset state variables
    // same as ResetForSim() if no special MaxInitTimeStep was defined
    ResetForInit();
}

TFloat CBTiming::GetTimeStep() {
    return timeStep_;
}

TFloat CBTiming::GetStartTime() {
    return startTime_;
}

TFloat CBTiming::GetStopTime() {
    return stopTime_;
}

/// returns the largest time step size possible during simulation (could be used for export)
TFloat CBTiming::GetMaxSimTimeStep() {
    return maxSimTimeStep_;
}

/// returns the value previously set using SetCurrentTime()
TFloat CBTiming::GetCurrentTime() {
    return currentTime_;
}

/// stores the time t in the variable currentTime
void CBTiming::SetCurrentTime(TFloat t) {
    currentTime_ = t;
}

/// updates the internal solver status, should be set before calling ComputeNextTimeStep()
void CBTiming::SetSolverStatus(CBStatus solverStatus) {
    solverStatus_ = solverStatus;
}

/// updates the internal variables for timeStep_ and counter depending on the
/// solvers return code. Returns FAILED if minimum time step was reached,
/// SUCCESS in all other cases.
CBStatus CBTiming::ComputeNextTimeStep() {
    bool shallReduce   = false;
    bool shallIncrease = false;
    bool shallRepeat   = false;
    
    if (solverStatus_ == CBStatus::FAILED)
        shallReduce = true;
    if (solverStatus_ == CBStatus::SUCCESS)
        shallIncrease = true;
    if (solverStatus_ == CBStatus::REPEAT)
        shallRepeat = true;
    
    // reduce time step size
    if (shallReduce) {
        lowestLevelCounter_ = 0;
        counter_ *= 2;
        timeStep_ /= 2;
        if (timeStep_ < minTimeStep_) {
            DCCtrl::print << "\n\nMinimum time step size has been reached with no success: Simulation has failed !!!\n";
            return CBStatus::FAILED;
        }
        DCCtrl::verbose << "\n\nRepeating with reduced time step size: " << GetTimeStep() << "\n";
        return CBStatus::SUCCESS;
    }
    
    // advance to next iteration, increase timestep if possible
    else if (shallIncrease) {
        counter_++;
        lowestLevelCounter_++;
        bool didIncrease = false;
        
        // the modulo division ensures that increase happens at time steps that fit the 'packets of four'
        // lowestLevelCounter_ is needed because of bug 2)
        while (counter_%2 == 0 && counter_%minNumCnt_ == 0 && lowestLevelCounter_ >= minNumCnt_) {
            counter_ = counter_/2;
            if (timeStep_ < maxInitTimeStep_) {
                timeStep_ *= 2;
                
                // fast relaxation: Going just once through the while loop prevents the
                // timestep to increase by multiple levels at the same time. Up to
                // maximal time step length at every export step, which can be a
                // unwanted effect if larger steps are very likely to fail.
                if (!fastRelaxation_)
                    break;
            }
        }
        if (didIncrease)
            DCCtrl::verbose << "\n\nIncreasing time step size: " << GetTimeStep() << "s \n";
        return CBStatus::SUCCESS;
    }
    
    // repeat with same time step -> nothing to do to the internal variables
    else if (shallRepeat) {
        DCCtrl::verbose << "\n\nRepeating with same time step: " << GetTimeStep() << "s \n";
        return CBStatus::SUCCESS;
    } else {
        throw std::runtime_error("CBTiming::ComputeNextTimeStep() Solver returned an unknown return code!");
    }
    return CBStatus::SUCCESS;
} // CBTiming::ComputeNextTimeStep

/// reset state variables for the initialization phase, different from ResetForSim() is the MaxTimeStep
void CBTiming::ResetForInit() {
    maxTimeStep_ = maxInitTimeStep_;
    timeStep_ = std::min(firstTimeStep_, maxTimeStep_);
    counter_ = 0;
    lowestLevelCounter_ = 0;
}

/// reset state variables for the simulation phase, different from ResetForInit() is the MaxTimeStep
void CBTiming::ResetForSim() {
    maxTimeStep_ = maxSimTimeStep_;
    timeStep_ = std::min(firstTimeStep_, maxTimeStep_);
    counter_ = 0;
    lowestLevelCounter_ = 0;
}

/// resets the internal counter variables, but not the timeStep size. Can be useful for ignoring certain simulation steps.
void CBTiming::ResetCounter() {
    counter_ = 0;
    lowestLevelCounter_ = 0;
}

bool CBTiming::IsValid(TFloat time) {
    return time < stopTime_ && stopTime_-time > 1e-14;
}
