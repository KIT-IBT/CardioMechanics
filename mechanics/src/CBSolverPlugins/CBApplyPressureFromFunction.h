//
//  CBApplyPressureFromFunction.h
//  CardioMechanics
//
//  Created by Steffen Schuler on 01.03.16.
//
//

#ifndef CB_APPLY_PRESSURE_FROM_FUNCTION
#define CB_APPLY_PRESSURE_FROM_FUNCTION

#include <map>
#include <memory>

#include "CBSolverPlugin.h"
#include "CBDataFromFunction.h"
#include "CBCirculationCavity.h"
#include "CBSolverPluginCavities.h"

class CBApplyPressureFromFunction : public CBSolverPluginCavities
{
public:
    CBApplyPressureFromFunction(){}
    virtual ~CBApplyPressureFromFunction(){}
    
    virtual void Init() override;
    virtual void Apply(TFloat time) override;
    virtual void ApplyToNodalForces() override;
    virtual void ApplyToNodalForcesJacobian() override;
    virtual void AnalyzeResults() override;
    virtual void WriteHeaderToFile();
    virtual void WriteToFile(TFloat time) override;
    bool WantsToAnalyzeResults() { return true; }
    CBStatus GetStatus() override {return status_;}
    virtual std::string GetName() override { return("ApplyPressureFromFunction"); };
    
private:
    std::ofstream file_;
    std::string filename_;
    bool headerWritten_ = false;
    
    struct IntervalStruct
    {
        std::shared_ptr<CBDataFromFunction> function;
        bool relaxElementsAtStart;
        std::vector<TInt> materialsToRelax;
        bool invert;
        TFloat startTime;
        TFloat stopTime;
        TFloat offset;
        TFloat amplitude;
    };
    
    struct ValueStruct
    {
        TFloat pressure;
        TFloat volume;
    };
    
    std::map<TInt, std::vector<IntervalStruct>> groups_;
    std::map<TInt, std::vector<IntervalStruct>>::iterator groupsIt_;
    std::map<TInt, ValueStruct> valuesStructs_;
    std::map<TInt, ValueStruct>::iterator valuesIt_;
    
    std::map<TInt, TFloat> materialsRelaxed_;
    std::map<TInt, TFloat>::iterator materialsRelaxedIt_;
    
    typedef CBSolverPlugin Base;
};

#endif
