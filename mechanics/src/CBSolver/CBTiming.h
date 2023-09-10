//
//  CBTiming.hpp
//  CardioMechanics
//
//  Created by Emanuel Poremba on 16.01.17.
//
//

#ifndef CBTiming_hpp
#define CBTiming_hpp

#include "DCType.h"
#include "ParameterMap.h"
#include "CBStatus.h"

#include <vector>

class CBTiming {
protected:
    TFloat startTime_ = 0;  //< global time when the simulation start
    TFloat stopTime_  = 0;  //< global time at which the simulation finished successfully
    TFloat currentTime_ = 0;  //< current time variable, gets only updated on SetCurrentTime()
    
    TFloat firstTimeStep_   = 0;  //< time step size used for the first initialization/simulation steps
    TFloat minTimeStep_     = 0;  //< reducing below this value will make the simulation fail
    TFloat maxInitTimeStep_ = 0;  //< max time step size for the INITIALIZATION phase
    TFloat maxSimTimeStep_  = 0;  //< max time step size for the SIMULATION phase
    TInt   minNumCnt_ = 0;        //< number of steps that needs to be done before allowing time steps to increase
    bool   fastRelaxation_ = false; //< allow time steps to increase by more than one level
    
    // state variables
    TFloat timeStep_    = 0;  //< variable countaining the current time step size
    TFloat maxTimeStep_ = 0;  //< max time step size, changes between initialization and simulation using ResetForSim()
    unsigned long long counter_            = 0;  //< global counter tracking reduction level and number of steps at the same time
    unsigned long long lowestLevelCounter_ = 0;  //< specific counter ensuring a certain number of smaller steps is done before trying a larger one
    
    CBStatus solverStatus_;
    
public:
    CBTiming();
    ~CBTiming();
    
    void Init(ParameterMap& parameter);
    
    TFloat GetTimeStep();
    TFloat GetStartTime();
    TFloat GetStopTime();
    TFloat GetMaxSimTimeStep();
    
    void SetCurrentTime(TFloat t);
    TFloat GetCurrentTime();
    
    void SetSolverStatus(CBStatus solverStatus);
    CBStatus ComputeNextTimeStep();
    
    void ResetForInit();
    void ResetForSim();
    void ResetCounter();
    
    bool IsValid(TFloat time);
    void CheckTime(TFloat time, TFloat timestep);
};

#endif /* CBTiming_hpp */
