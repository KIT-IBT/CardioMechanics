/*
 * File: CBApplyPressureFromFunctionNodeExport.h
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


#ifndef CB_APPLY_PRESSURE_FROM_FUNCTION_NODE_EXPORT
#define CB_APPLY_PRESSURE_FROM_FUNCTION_NODE_EXPORT

#include <map>
#include <memory>

#include "CBSolverPlugin.h"
#include "CBDataFromFunction.h"
#include "CBCirculationCavity.h"
#include "CBSolverPluginCavities.h"

class CBApplyPressureFromFunctionNodeExport : public CBSolverPluginCavities
{
public:
    CBApplyPressureFromFunctionNodeExport(){}
    virtual ~CBApplyPressureFromFunctionNodeExport(){}
    
    virtual void Init() override;
    virtual void Apply(TFloat time) override;
    virtual void ApplyToNodalForces() override;
    virtual void ApplyToNodalForcesJacobian() override;
    virtual void AnalyzeResults() override;
    virtual void WriteHeaderToFile();
    virtual void WriteToFile(TFloat time) override;
    bool WantsToAnalyzeResults() { return true; }
    bool ExitCheck(TFloat time) { return shallExit_; };
    CBStatus GetStatus() override {return status_;}
    virtual std::string GetName() override { return("ApplyPressureFromFunctionNodeExport"); };
    
private:
    TInt stopSurface_;
    TFloat stopVolume_;
    TFloat nodeExportTime_;
    bool shallExit_ = false;
    std::string nodeFilename_;
    std::ofstream nodeFile_;
    Vec coords_;
    Vec coordsSeq_;
    VecScatter coordsScatter_;
    void InitPetscVectors();
    void ExportNodeFile();
    TFloat currentTime_;
    
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

