//
//  CBApplyPressureFromFunction.cpp
//  CardioMechanics
//
//  Created by Steffen Schuler on 01.03.16.
//
//

#include "CBApplyPressureFromFunction.h"
#include "CBSolver.h"
#include "CBElementSurfaceT3.h"
#include "CBElementSurfaceT6.h"

void CBApplyPressureFromFunction::Init()
{
    // init filenames
    std::string defaultFilename = adapter_->GetSolver()->GetModel()->GetExporter()->GetExportDir() + "/ApplyPressureFromFunction.dat";
    filename_ = parameters_->Get<std::string>("Plugins.ApplyPressureFromFunction.ExportFile", defaultFilename);
    
    // init cavities map using function from base class
    CBSolverPluginCavities::InitCavities("Plugins.ApplyPressureFromFunction");
    
    // init groups and intervals
    std::vector<std::string> groupKeys = parameters_->GetChildNodes("Plugins.ApplyPressureFromFunction.Groups");
    for (auto &groupKey : groupKeys)
    {
        std::vector<TInt> surfaces = parameters_->GetArray<TInt>(groupKey + ".Surfaces");
        
        std::vector<std::string> intervalKeys = parameters_->GetChildNodes(groupKey + ".Intervals");
        if (intervalKeys.size() == 0)
            throw std::runtime_error("void CBApplyPressureFromFunction::Init(): No Intervals specified in " + groupKey);
        
        TFloat prevStopTime = 0.0;
        std::vector<IntervalStruct> intervals;
        for (auto &intervalKey : intervalKeys)
        {
            IntervalStruct interval;
            
            interval.function  = std::make_shared<CBDataFromFunction>(parameters_, intervalKey);
            interval.relaxElementsAtStart = parameters_->Get<bool>(intervalKey + ".RelaxElementsAtStart", false);
            if (interval.relaxElementsAtStart)
                interval.materialsToRelax = parameters_->GetArray<TInt>(intervalKey + ".MaterialsToRelax");
            interval.invert    = parameters_->Get<bool>(intervalKey + ".Invert", false);
            interval.startTime = parameters_->Get<TFloat>(intervalKey + ".StartTime", 0.0);
            if (interval.startTime < prevStopTime)
                throw std::runtime_error("void CBApplyPressureFromFunction::Init(): StartTime of " + intervalKey + " is smaller than StopTime of previous Interval");
            interval.stopTime  = parameters_->Get<TFloat>(intervalKey + ".StopTime", 0.0);
            prevStopTime = interval.stopTime;
            interval.offset    = parameters_->Get<TFloat>(intervalKey + ".Offset", 0.0);
            interval.amplitude = parameters_->Get<TFloat>(intervalKey + ".Amplitude", 0.0);
            
            intervals.push_back(interval);
        }
        
        for (auto &surface : surfaces)
        {
            if (cavities_.find(surface) == cavities_.end())
                throw std::runtime_error("void CBApplyPressureFromFunction::Init(): No cavity with surface index " + std::to_string(surface) + " found.");
            
            groups_.insert(std::pair<TInt, std::vector<IntervalStruct>>(surface, intervals));
            ValueStruct vs;
            vs.pressure = 0.0;
            vs.volume = 0.0;
            valuesStructs_.insert(std::pair<TInt, ValueStruct>(surface, vs));
        }
    }
}


void CBApplyPressureFromFunction::Apply(TFloat time)
{
    auto& timing = adapter_->GetSolver()->GetTiming();
    for (groupsIt_ = groups_.begin(); groupsIt_ != groups_.end(); groupsIt_++)
    {
        for (auto &interval : groupsIt_->second)
        {
            TFloat dt = timing.GetTimeStep();
            if (time >= interval.startTime && time < interval.stopTime+0.5*dt)
            {
                if (interval.relaxElementsAtStart)
                {
                    for (auto &mat : interval.materialsToRelax)
                    {
                        materialsRelaxedIt_ = materialsRelaxed_.find(mat);
                        if (materialsRelaxedIt_ == materialsRelaxed_.end() || materialsRelaxedIt_->second != interval.startTime)
                        {
                            DCCtrl::print << "Relaxing elements with material " << mat << std::endl;
                            Base::adapter_->GetSolver()->RelaxElementsAndBasesByMaterialIndex(mat);
                            materialsRelaxed_.insert(std::pair<TInt, TFloat>(mat, interval.startTime));
                        }
                    }
                    interval.relaxElementsAtStart = false;
                }
                
                TFloat f = interval.function->Get(time, 0);
                if (interval.invert)
                    valuesStructs_[groupsIt_->first].pressure = interval.offset + interval.amplitude * (1.0-f);
                else
                    valuesStructs_[groupsIt_->first].pressure = interval.offset + interval.amplitude * f;
                
                break;
            }
        }
    }
}


void CBApplyPressureFromFunction::ApplyToNodalForces()
{
    for (valuesIt_ = valuesStructs_.begin(); valuesIt_ != valuesStructs_.end(); valuesIt_++)
        cavities_[valuesIt_->first]->ApplyPressure(valuesIt_->second.pressure);
}


void CBApplyPressureFromFunction::ApplyToNodalForcesJacobian()
{
    for (valuesIt_ = valuesStructs_.begin(); valuesIt_ != valuesStructs_.end(); valuesIt_++)
        cavities_[valuesIt_->first]->ApplyPressureToJacobian(valuesIt_->second.pressure);
}


void CBApplyPressureFromFunction::AnalyzeResults()
{
    for (valuesIt_ = valuesStructs_.begin(); valuesIt_ != valuesStructs_.end(); valuesIt_++)
        valuesIt_->second.volume = cavities_[valuesIt_->first]->CalcVolume();
}


void CBApplyPressureFromFunction::WriteHeaderToFile()
{
    file_.open(filename_.c_str());
    if (!file_.good())
        throw std::runtime_error("void CBApplyPressureFromFunction::WriteToFile(TFloat time): Couldn't create " + filename_);
    file_ << std::setw(16) << "time";
    for (valuesIt_ = valuesStructs_.begin(); valuesIt_ != valuesStructs_.end(); valuesIt_++)
    {
        std::ostringstream p, v;
        p << "pressure" << valuesIt_->first;
        v << "volume" << valuesIt_->first;
        file_ << std::setw(16) << p.str() << std::setw(16) << v.str();
    }
    file_ << std::endl;
    file_.close();
}


void CBApplyPressureFromFunction::WriteToFile(TFloat time)
{
    // write header on first call
    if (!headerWritten_)
        WriteHeaderToFile();
    headerWritten_ = true;
    
    // write plugin data
    file_.open(filename_.c_str(), std::ios::app);
    if (!file_.good())
        throw std::runtime_error("void CBApplyPressureFromFunction::WriteToFile(TFloat time): Couldn't create " + filename_);
    file_ << std::setprecision(7) << std::fixed << std::setw(16) << time;
    for (valuesIt_ = valuesStructs_.begin(); valuesIt_ != valuesStructs_.end(); valuesIt_++)
    {
        file_ << std::setw(16) << valuesIt_->second.pressure/133.322 << std::setw(16) << (valuesIt_->second).volume*1e6;
        std::cout << std::setw(16) << valuesIt_->second.pressure/133.322 << std::setw(16) << (valuesIt_->second).volume*1e6;
    }
    std::cout << std::endl;
    file_ << std::endl;
    file_.close();
}
