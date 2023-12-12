/*
 * File: CBRuntimeEstimator.h
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



/// this class provides and displays estimation of a simulation's runtime
class CBRuntimeEstimator {
protected:
    CBTiming timing_;
    bool enableOutput_;
    PetscScalar wallStartTime_;
    PetscScalar wallCurrentTime_;
    PetscScalar wallLastOutputTime_;
    
public:
    CBRuntimeEstimator(CBTiming &timing) : timing_(timing) {
        enableOutput_ = true;
        wallStartTime_ = MPI_Wtime();
        wallCurrentTime_ = wallStartTime_;
        wallLastOutputTime_ = wallStartTime_;
    }
    
    /// allows to disable output by passing false. Output is enabled by default.
    void EnableOutput(bool enableOutput) {
        enableOutput_ = enableOutput;
    }
    
    /// convert time duration from (int)seconds to a text string
    std::string ConvertSecondsToHMMSS(TFloat i) {
        int hour    = int(i / 3600);
        int minutes = int((i - (hour * 3600)) / 60);
        int seconds = int(i - (hour * 3600) - (minutes * 60));
        
        return std::string(std::to_string(hour)+":"+ std::string((minutes < 10) ? "0" : "") + std::to_string(
                                                                                                             minutes)+":"+ std::string((seconds < 10) ? "0" : "") +std::to_string(seconds));
    }
    
    /// store the current wall clock time as start of the simulation
    void Tic() { wallStartTime_ = MPI_Wtime(); }
    
    /// set the current wall clock time
    void Toc() { wallCurrentTime_ = MPI_Wtime(); }
    
    /// show runtime estimation using simTime and the current wall clock time to estimate the progress of the simulation, Toc() needs to be called manually before
    void ShowEstimation(TFloat simTime) {
        Toc();
        if ((wallCurrentTime_ - wallLastOutputTime_ > 1) && (simTime != timing_.GetStartTime())) {
            PetscScalar percentDone                 = 100.0 / (timing_.GetStopTime() - timing_.GetStartTime()) *
            (simTime - timing_.GetStartTime());
            PetscScalar runningTime                 = (wallCurrentTime_ - wallStartTime_);
            PetscScalar estimatedTimeToRun          = (100 / percentDone) * runningTime;
            PetscScalar estimatedRemainingTimeToRun = estimatedTimeToRun - runningTime;
            if (enableOutput_) {
                DCCtrl::print << "\n\n" << percentDone << "% done !! Estimated time to run: " << ConvertSecondsToHMMSS(
                                                                                                                       estimatedTimeToRun) << ". Estimated remaining time: " << ConvertSecondsToHMMSS(estimatedRemainingTimeToRun) <<
                "      \n\n";
            }
            wallLastOutputTime_ = MPI_Wtime();
        }
    }
    
    /// show a final summary, to be called when the simulation finished successfully, needs Toc() to be called manually before
    void ShowFinalSummary() {
        Toc();
        if (enableOutput_) {
            DCCtrl::print << "\n\n\nSimulation finished successfully. [simulation time: " << ConvertSecondsToHMMSS(
                                                                                                                   wallCurrentTime_ - wallStartTime_)<< "] \n";
        } else {
            DCCtrl::print << "\n\n\nSimulation finished successfully.\n";
        }
    }
};
