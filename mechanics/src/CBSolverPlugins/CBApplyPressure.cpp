//
//  CBApplyPressure.cpp
//  CardioMechanics
//
//  Created by Thomas Fritz on 15.10.13.
//
//
#include "DCCtrl.h"
#include "CBApplyPressure.h"

#include "CBSolver.h"
#include "CBElementSurfaceT3.h"
#include "CBElementSurfaceT6.h"
#include "CBElementSolidT10.h"

void CBApplyPressure::Init()
{
    auto& timing = adapter_->GetSolver()->GetTiming();
    TFloat startTimeDefault = timing.GetStartTime();
    TFloat stopTimeDefault = timing.GetStopTime();
    maxPressure_ = parameters_->Get<TFloat>("Plugins.ApplyPressure.MaxPressure",100);
    surface_ = parameters_->Get<TFloat>("Plugins.ApplyPressure.SurfaceIndex");
    startTime_ = parameters_->Get<TFloat>("Plugins.ApplyPressure.StartTime",startTimeDefault);
    stopTime_ =parameters_->Get<TFloat>("Plugins.ApplyPressure.StopTime",stopTimeDefault);
    keepMaxPressure_ = parameters_->Get<bool>("Plugins.ApplyPressure.KeepMaxPressure",true);
    currentVolume_ = CalcVolume();
    
    // extract cavity elements
    for(auto& it : adapter_->GetSolver()->GetElementVector()) {
        CBElementCavity* s = dynamic_cast<CBElementCavity*>(it);
        if (s!=0)
            if (s->GetSurfaceIndex() == surface_ )
                targetElements_.push_back(s);
    }
    
    if (targetElements_.size() == 0)
        DCCtrl::print << "CBApplyPressure WARNING: The number of target elements is zero. This might not be the intention." << std::endl;
    DCCtrl::print << "CBApplyPressure target elements size: " << targetElements_.size() << std::endl;
    
}

void CBApplyPressure::StepBack()
{
    volumeWork_ = lastVolumeWork_;
}

void CBApplyPressure::Apply(TFloat time)
{
    TFloat lastVolume = currentVolume_;
    currentVolume_ = CalcVolume();
    
    currentPressure_ = 0;
    if(keepMaxPressure_)
        currentPressure_ = maxPressure_;
    if(time < stopTime_)
        currentPressure_ = maxPressure_  *  time / ( stopTime_ - startTime_ );
    
    lastVolumeWork_ = volumeWork_;
    volumeWork_ += (currentVolume_ - lastVolume)*(currentPressure_ );
    
    Petsc::debug << "ApplyPressure(p/V/W): " << currentPressure_ << " " << currentVolume_ << " " << volumeWork_ << "\n";
}

void CBApplyPressure::ApplyToNodalForces()
{
    for (auto s : targetElements_)
        s->ApplyPressure(currentPressure_);
}

TFloat CBApplyPressure::CalcVolume()
{
    double volume      = 0;
    double localVolume = 0;
    
    for (auto s : targetElements_)
        localVolume+=s->CalcContributionToVolume();
    
    MPI_Allreduce(&localVolume, &volume, 1, MPI_DOUBLE, MPI_SUM, Petsc::Comm());
    return(volume);
}

void CBApplyPressure::WriteToFile(TFloat time)
{
    std::string defaultFilename = adapter_->GetSolver()->GetModel()->GetExporter()->GetExportDir() + "/ApplyPressure.dat";
    std::string filename = parameters_->Get<std::string>("Plugins.ApplyPressure.Filename", defaultFilename);
    
    if(filename != "" || GetAdapter()->GetSolver()->GetModel()->GetExporter()->GetExportDir() != "")
    {
        std::ofstream file;
        if(time == startTime_)
            file.open(filename.c_str());
        else
            file.open(filename.c_str(), std::ios::app);
        
        file << time << " " << currentPressure_ << " " << currentVolume_ << " " << volumeWork_ << "\n";
        file.close();
    }
}

void CBApplyPressure::ApplyToNodalForcesJacobian()
{
    for (auto s : targetElements_)
        s->CalcPressureJacobian(currentPressure_);
}
