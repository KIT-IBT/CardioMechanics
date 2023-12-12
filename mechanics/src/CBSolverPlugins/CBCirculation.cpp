/*
 * File: CBCirculation.cpp
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


#include "CBSolver.h"
#include "CBElementSurfaceT6.h"
#include "CBCirculation.h"
#include "CBCircWindkessel3.h"
#include "CBCircWindkessel3AsymptoticPressure.h"
#include "CBCircOneVentricle.h"
#include "CBCircTwoVentricles.h"
#include "CBCircLeftHeart.h"
#include "CBCircWholeHeart.h"
#include "CBCircWholeHeartValves.h"
#include "CBCircWholeHeartValvesDynamic.h"
#include "CBCircIsovolumetric.h"

CBCirculation::CBCirculation() : CBSolverPluginCavities() {}

CBCirculation::~CBCirculation() {
    for (auto &circModel : circModels_)
        delete circModel;
}

void CBCirculation::Init() {
    // cavities_ are handled by CBApplyPressureFromFunction, reduces duplicate code
    CBSolverPluginCavities::InitCavities("Plugins.Circulation");
    
    std::string defaultFilename = adapter_->GetSolver()->GetModel()->GetExporter()->GetExportDir() +
    "/CirculationIterations.dat";
    itsFilename_                   = parameters_->Get<std::string>("Plugins.Circulation.IterationsFile", defaultFilename);
    startTime_                     = parameters_->Get<TFloat>("Plugins.Circulation.StartTime",
                                                              std::numeric_limits<double>::lowest());
    if (startTime_ != std::numeric_limits<double>::lowest())
        status_ = CBStatus::DACCORD;
    
    tolerance_                     = parameters_->Get<TFloat>("Plugins.Circulation.Tolerance", 0.0001);
    maxIterations_                 = parameters_->Get<TInt>("Plugins.Circulation.MaxIterations", 100);
    nNonNewton_                    = parameters_->Get<TInt>("Plugins.Circulation.MaxNonNewtonIterations", 10);
    nSecant_                       = parameters_->Get<TInt>("Plugins.Circulation.SecantIterations", 3);
    perturbation_                  = parameters_->Get<bool>("Plugins.Circulation.Perturbation", true);
    preloadingTime1_               = parameters_->Get<TFloat>("Plugins.Circulation.PreloadingTime1", 0.0);
    preloadingTime2_               = parameters_->Get<TFloat>("Plugins.Circulation.PreloadingTime2", 0.0);
    dampingFactor_                 = parameters_->Get<double>("Plugins.Circulation.CouplingDamping.InitialFactor", 0.0);
    if (dampingFactor_ > 0.9)
        dampingFactor_ = 0.9;
    dampingActive_ = (dampingFactor_ < 1e-5) ? false : true;
    dampingDeclineFactor_          = parameters_->Get<double>("Plugins.Circulation.CouplingDamping.DeclineFactor", 0.0);
    
    std::string pressureUnit = parameters_->Get<std::string>("Plugins.Circulation.PressureUnit");
    if (pressureUnit == "mmHg")
        uP_ = 1.33322e2;
    else if (pressureUnit == "hPa")
        uP_ = 1e2;
    else if (pressureUnit == "Pa")
        uP_ = 1.0;
    else
        throw std::runtime_error("void CBCirculation::InitCircModels(): Pressure Unit " + pressureUnit + " unknown.");
    
    maxP_   = 1.33322e3  / uP_; // 10 mmHg
    deltaP_ = 1.33322    / uP_; // 0.01 mmHg
    dP_     = 1.33322e-6 / uP_; // 1e-8 mmHg
    
    std::string volumeUnit = parameters_->Get<std::string>("Plugins.Circulation.VolumeUnit");
    if (volumeUnit == "ml")
        uV_ = 1e6;
    else if (volumeUnit == "m3")
        uV_ = 1.0;
    else
        throw std::runtime_error("void CBCirculation::InitCircModels(): Volume Unit " + volumeUnit + " unknown.");
    
    itsFile_.open(itsFilename_.c_str());
    if (!itsFile_.good())
        throw std::runtime_error("CBCirculation::Init: Couldn't create " + itsFilename_ + ".");
    itsFile_.close();
    
    //    /////////////// Export file for circVolume, fineVolume and pressure in each newton iteration
    //    filename_ = adapter_->GetSolver()->GetModel()->GetExporter()->GetExportDir() + "/CirculationResidual.dat";;
    //    if(filename_ != "")
    //    {
    //        size_t pos = filename_.find_last_of("/\\");
    //        if(pos != std::string::npos)
    //        {
    //            std::string dir = filename_.substr(0, pos);
    //            frizzle::filesystem::CreateDirectory(dir);
    //        }
    //
    //        std::ofstream file(filename_.c_str());
    //        if(!file.good())
    //            throw std::runtime_error("CBCirculation::Init: Couldn't create " + filename_ + ".");
    //        file << "CurrentIter\tCircVolume1\tFineVolume1\tPressure1\tCircVolume2\tFineVolume2\tPressure2\tCircVolume3\tFineVolume3\tPressure3\tCircVolume4\tFineVolume4\tPressure4\n";
    //        file.close();
    //    }
    //    else
    //        throw std::runtime_error("CBCirculation::Init: Filename is empty.");
    //    ///////////////
    
    dt_ = Base::adapter_->GetSolver()->GetTiming().GetTimeStep();
    
    InitCircModels();
} // CBCirculation::Init

void CBCirculation::InitCircModels() {
    std::vector<std::string> circs = parameters_->GetChildNodes("Plugins.Circulation.Circs");
    
    for (auto &circ : circs) {
        if (parameters_->Get<bool>(circ + ".Active", false)) {
            CBCircModelBase *circModel = 0;
            
            std::string type = parameters_->Get<std::string>(circ + ".Type", "");
            if (type == "CircWindkessel3") {
                circModel = new CBCircWindkessel3(parameters_, circ);
                circModels_.push_back(circModel);
            } else if (type == "CircWindkessel3AsymptoticPressure") {
                circModel = new CBCircWindkessel3AsymptoticPressure(parameters_, circ);
                circModels_.push_back(circModel);
            } else if (type == "CircOneVentricle") {
                circModel = new CBCircOneVentricle(parameters_, circ);
                circModels_.push_back(circModel);
            } else if (type == "CircTwoVentricles") {
                circModel = new CBCircTwoVentricles(parameters_, circ);
                circModels_.push_back(circModel);
            } else if (type == "CircWholeHeart") {
                circModel = new CBCircWholeHeart(parameters_, circ);
                circModels_.push_back(circModel);
            } else if (type == "CircWholeHeartValves") {
                circModel = new CBCircWholeHeartValves(parameters_, circ);
                circModels_.push_back(circModel);
            } else if (type == "CircWholeHeartValvesDynamic") {
                circModel = new CBCircWholeHeartValvesDynamic(parameters_, circ);
                circModels_.push_back(circModel);
            } else if (type == "CircLeftHeart") {
                circModel = new CBCircLeftHeart(parameters_, circ);
                circModels_.push_back(circModel);
            } else if (type == "CircIsovolumetric") {
                circModel = new CBCircIsovolumetric(parameters_, circ);
                circModels_.push_back(circModel);
            } else {
                throw std::runtime_error("void CBCirculation::InitCircModels(): Circ Model type " + type + " is unknown.");
            }
        }
    }
    
    for (auto &circModel : circModels_) {
        TInt nC = circModel->GetNumberOfCavities();
        TInt *cavInd = new TInt[nC];
        circModel->GetCavitySurfaceIndices(cavInd);
        
        for (TInt i = 0; i < nC; i++) {
            if (std::find(cavitySurfaceIndices_.begin(), cavitySurfaceIndices_.end(),
                          cavInd[i]) != cavitySurfaceIndices_.end()) {
                throw std::runtime_error("void CBCirculation::InitCircModels(): Cavity " + std::to_string(
                                                                                                          cavInd[i]) + " is used more than once.");
            } else {
                cavitySurfaceIndices_.push_back(cavInd[i]);
            }
        }
        
        nC_ += nC;
        delete[] cavInd;
    }
    
    // check, whether a cavity exists for each surface index listed in the XML
    for (auto cavitySurfaceIndex : cavitySurfaceIndices_)
        if (cavities_.find(cavitySurfaceIndex) == cavities_.end()) {
            throw std::runtime_error("void CBCirculation::InitCircModels(): No cavity with surface index " +
                                     std::to_string(cavitySurfaceIndex) + " found.");
        }
    
    // allocate vectors and matrices (resize them to the total number of cavities)
    estPressure_.Resize(nC_);
    prevEstPressure_.Resize(nC_);
    pressure_.resize(N_);
    for (TInt i = 0; i < N_; i++)
        pressure_[i].Resize(nC_);
    fineVolume_.Resize(nC_);
    prevFineVolume_.Resize(nC_);
    circVolume_.Resize(nC_);
    prevCircVolume_.Resize(nC_);
    
    preloadingPressure1_.Resize(nC_);
    preloadingPressure2_.Resize(nC_);
    
    residual_.Resize(nC_);
    prevResidual_.Resize(nC_);
    
    fineCompli_.Resize(nC_);
    circCompli_.Resize(nC_);
    sysCompli_.Resize(nC_);
    sysElast_.Resize(nC_);
    sysElast_.SetToIdentityMatrix();
    
    // Get initial cavity volumes from FE model and check if cavity surfaces are closed
    std::string unclosedCavities = "";
    for (TInt i = 0; i < nC_; i++) {
        TInt cavInd = cavitySurfaceIndices_[i];
        fineVolume_(i) = uV_ * cavities_[cavInd]->CalcVolume();
        circVolume_(i) = fineVolume_(i);
        
        // Closed surface?
        if (!cavities_[cavInd]->ClosedSurfaceCheck())
            unclosedCavities += std::to_string(cavInd) + " ";
    }
    
    if (unclosedCavities != "") {
        throw std::runtime_error(
                                 "void CBCirculation::InitCircModels(): The following cavity surfaces are probably not closed: " +
                                 unclosedCavities);
    }
    
    // Set initial volumes of circModels and get preloading pressures and initial pressures from circModels
    TInt cavPos = 0;
    for (auto &circModel : circModels_) {
        TInt nC = circModel->GetNumberOfCavities();
        
        TFloat *preloadingPressure = new TFloat[2*nC];
        circModel->GetPreloadingPressures(preloadingPressure);
        
        TFloat *circVolume = new TFloat[nC];
        for (TInt i = 0; i < nC; i++) {
            circVolume[i] = fineVolume_(cavPos+i);
            preloadingPressure1_(cavPos+i) = preloadingPressure[i];
            preloadingPressure2_(cavPos+i) = preloadingPressure[nC+i];
        }
        
        circModel->SetInitialCavityVolumes(circVolume);
        
        cavPos += nC;
        delete[] preloadingPressure;
        delete[] circVolume;
    }
} // CBCirculation::InitCircModels

void CBCirculation::Prepare() {
    // do preparation
    status_ = CBStatus::PREPARING_SIMULATION;
}

double CBCirculation::GetPreparationProgress() {
    double progress = 1.0;
    
    // calculate preparation progress
    return progress;
}

bool CBCirculation::ExitCheck(TFloat time) {
    if (!hasStarted_)
        return false;
    
    bool allInSteadyState = true;
    for (auto &circModel : circModels_) {
        bool inSteadyState = circModel->SteadyStateCheck(time);
        if (allInSteadyState && !inSteadyState)
            allInSteadyState = false;
    }
    if (allInSteadyState)
        DCCtrl::print << "\nAll circulations in steady state." << std::endl;
    else
        DCCtrl::debug << "\nNot all circulations in steady state yet." << std::endl;
    
    return allInSteadyState;
}

void CBCirculation::WriteToFile(PetscScalar time) {
    if (!hasStarted_)
        return;
    
    if (time > preloadingTime1_+preloadingTime2_+dt_) {
        for (auto &circModel : circModels_)
            circModel->WriteResultsToFile(time);
        
        itsFile_.open(itsFilename_.c_str(), std::ios::app);
        if (!itsFile_.good())
            throw std::runtime_error("CBCirculation::WriteToFile: Couldn't create " + itsFilename_ + ".");
        itsFile_ << std::setprecision(3) << std::fixed << std::setw(14) << time << std::setw(14) << i_ << std::endl;
        itsFile_.close();
    }
}

void CBCirculation::GetCircVolume(VectorN<TFloat> *pressure, VectorN<TFloat> *volume) {
    TInt cavPos = 0;
    
    for (auto &circModel : circModels_) {
        TInt nC = circModel->GetNumberOfCavities();
        TFloat *circPressure = new TFloat[nC];
        TFloat *circVolume = new TFloat[nC];
        
        for (TInt i = 0; i < nC; i++)
            circPressure[i] = (*pressure)(cavPos+i);
        
        circModel->GetEstCavityVolumes(circPressure, dt_, circVolume);
        
        for (TInt i = 0; i < nC; i++)
            (*volume)(cavPos+i) = circVolume[i];
        
        cavPos += nC;
        delete[] circPressure;
        delete[] circVolume;
    }
}

void CBCirculation::ApplyToNodalForces() {
    if (!hasStarted_)
        return;
    
    for (TInt i = 0; i < nC_; i++) {
        TInt cavInd = cavitySurfaceIndices_[i];
        cavities_[cavInd]->ApplyPressure(uP_ * estPressure_(i));
    }
}

void CBCirculation::ApplyToNodalForcesJacobian() {
    //    for(TInt i = 0; i < nC_; i++)
    //    {
    //        TInt cavInd = cavitySurfaceIndices_[i];
    //        cavities_[cavInd]->ApplyPressureToJacobian(uP_ * estPressure_(i));
    //    }
}

void CBCirculation::Apply(PetscScalar time) {
    if (!hasStarted_) {
        if (time > startTime_) {
            hasStarted_ = true;
        } else {
            return;
        }
    }
    
    prevDt_ = dt_;
    dt_ = Base::adapter_->GetSolver()->GetTiming().GetTimeStep();
    
    if (preloadingFinished_) {
        if (lastPreloadingStep2_ && false) { // disabled to allow adjustment of initial velocity and acceleration via preloadingTime2
                                             // set nodal velocity and acceleration to zero after preloading phase 2
            Base::adapter_->GetSolver()->SetZeroVelocityAndAcceleration();
            DCCtrl::print << "Setting zero velocity and acceleration" << std::endl;
            lastPreloadingStep2_ = false;
        }
        
        if (dampingIteration_) {
            estPressure_ = dampingFactor_*pressure_[(n_-1)%N_] + (1.0-dampingFactor_)*estPressure_;
            str1_ << "Damping Iteration (damping factor: " << dampingFactor_ << ")";
            
            dampingFactor_ *= dampingDeclineFactor_;
            if (dampingFactor_ < 1e-5)
                dampingActive_ = false;
        } else if (status_ == CBStatus::SUCCESS) {
            if (!stepBack_) {
                str2_ << "n = " << n_ << " | t = " << time-dt_ << " | curr its = " << i_ << " | avg its = " <<
                std::setprecision(2) << std::fixed << (double)iSum_/n_ << std::endl
                << "----------------------------------------------------------" << std::endl;
                
                // accept estimated pressure from previous timestep
                for (auto &circModel : circModels_)
                    circModel->AcceptEstState();
                pressure_[n_%N_] = estPressure_;
                timeStep_[n_%N_] = prevDt_;
                n_++;
                i_ = 0;
                j_ = 0;
            }
            
            // make a new pressure estimate for the current, new timestep
            if (n_ > N_) // extrapolate pressure, if there are enough values in the buffer
                ExtrapolatePressure();
            
            else // simply increase the pressure a little, if there are not yet enough values in the buffer
                IncreasePressure();
        } else {
            if (i_ == 0) {
                if (n_ > N_) {
                    ExtrapolatePressure();
                } else {
                    IncreasePressure();
                }
            } else if (perturbation_ && (i_ <= nSecant_) ) {      // or if nC_ == 1
                EstimatePressureUsingSecant();
            } else if (perturbation_ && (i_ <= nSecant_+nC_) ) {
                PerturbPressure(i_-nSecant_-1);
            } else if ((j_ > nNonNewton_) && (j_ <= nNonNewton_+nC_) ) {
                PerturbPressure(j_-nNonNewton_-1);
            } else if ((!perturbation_ || (i_ > nSecant_+nC_+1)) && (j_ <= nNonNewton_) ) {
                EstimatePressureUsingQuasiNewton();
            }
            
            if (perturbation_ || (j_ > nNonNewton_+1)) {
                if ((i_ > nSecant_+1) && (i_ <= nSecant_+nC_+1))
                    AssembleComplianceMatrices(i_-nSecant_-2);
                
                else if (j_ > nNonNewton_+1)
                    AssembleComplianceMatrices(j_-nNonNewton_-2);
                
                if ((i_ == nSecant_+nC_+1) || (j_ == nNonNewton_+nC_+1))
                    EstimatePressureUsingNewton();
            }
        }
        
        // limit the change in pressure
        for (TInt i = 0; preloadingFinished_ && i < nC_; i++) {
            TFloat estPressureDiff = estPressure_(i)-prevEstPressure_(i);
            if (estPressureDiff > maxP_) {
                estPressure_(i) = prevEstPressure_(i) + maxP_;
                str1_ << " lim" << i << "+";
            } else if (estPressureDiff < -maxP_) {
                estPressure_(i) = prevEstPressure_(i) - maxP_;
                str1_ << " lim" << i << "-";
            }
        }
        
        stepBack_ = false;
        i_++;
        j_++;
        iSum_++;
    } else {
        if (lastPreloadingStep1_) {
            // set nodal velocity and acceleration to zero after preloading phase 1
            Base::adapter_->GetSolver()->SetZeroVelocityAndAcceleration();
            DCCtrl::print << "Setting zero velocity and acceleration" << std::endl;
            lastPreloadingStep1_ = false;
        }
        
        if (time < preloadingTime1_) {
            estPressure_ = -time/preloadingTime1_ * preloadingPressure1_;
            str1_ << "Preloading Phase 1";
        } else if (std::abs(time-preloadingTime1_) < 0.5*dt_) {
            estPressure_ = -1.0 * preloadingPressure1_;
            lastPreloadingStep1_ = true;
            str1_ << "Preloading Phase 1: Last step";
        } else if (time <= preloadingTime1_+preloadingTime2_) {
            estPressure_ = (time-preloadingTime1_)/preloadingTime2_ * preloadingPressure2_;
            str1_ << "Preloading Phase 2";
        } else {
            estPressure_ = preloadingPressure2_;
            pressure_[0] = estPressure_;
            lastPreloadingStep2_ = true;
            str1_ << "Preloading Phase 2: Last step";
        }
    }
    
    // Set Status to avoid accidental repeat when the solver diverges.
    // This mainly is a problem for CBCirculation.
    status_ = CBStatus::SUCCESS;
    
    str2_ << std::setw(5) << i_ << std::setw(10) << std::setprecision(5) << dt_;
    for (TInt i = 0; i < nC_; i++)
        str2_ << std::setw(12) << std::setprecision(5) << std::fixed << estPressure_(i);
    DCCtrl::print << str2_.str() << "\t" << str1_.str() << std::endl;
    str1_.str(std::string());
    str2_.str(std::string());
} // CBCirculation::Apply

void CBCirculation::AnalyzeResults() {
    if (!hasStarted_)
        return;
    
    if (!preloadingFinished_) {
        if (lastPreloadingStep2_) {
            // reset initial cavity volumes after preloading
            for (TInt i = 0; i < nC_; i++) {
                TInt cavInd = cavitySurfaceIndices_[i];
                fineVolume_(i) = uV_ * cavities_[cavInd]->CalcVolume();
                circVolume_(i) = fineVolume_(i);
            }
            TInt cavPos = 0;
            for (auto &circModel : circModels_) {
                TInt nC = circModel->GetNumberOfCavities();
                TFloat *circVolume = new TFloat[nC];
                for (TInt i = 0; i < nC; i++)
                    circVolume[i] = fineVolume_(cavPos+i);
                circModel->SetInitialCavityVolumes(circVolume);
                cavPos += nC;
                delete[] circVolume;
            }
            
            preloadingFinished_ = true;
            status_ = CBStatus::DACCORD;
        } else if (lastPreloadingStep1_) {
            // relax elements
            if (preloadingTime1_ != 0.0) {
                Base::adapter_->GetSolver()->RelaxElementsAndBases();
                DCCtrl::print << "Relaxing elements" << std::endl;
            }
            
            status_ = CBStatus::SUCCESS;
        } else {
            status_ = CBStatus::SUCCESS;
        }
        return;
    }
    
    // get volume from FE model
    prevFineVolume_ = fineVolume_;
    for (TInt i = 0; i < nC_; i++) {
        TInt cavInd = cavitySurfaceIndices_[i];
        fineVolume_(i) = uV_ * cavities_[cavInd]->CalcVolume();
    }
    
    // get volume from circulation
    prevCircVolume_ = circVolume_;
    GetCircVolume(&estPressure_, &circVolume_);
    
    status_ = CBStatus::SUCCESS;
    
    DCCtrl::debug << "\nCBCirculation::AnalyzeResults(): Residual V_FE - V_CIRC\n";
    
    if (dampingIteration_) {
        dampingIteration_ = false;
    } else {
        std::vector<TFloat> residual(nC_);
        
        //        std::ofstream file(filename_.c_str(), std::ios::app);
        //        file << i_ << " ";
        //
        for (TInt i = 0; i < nC_; i++) {
            residual[i] = std::abs(fineVolume_(i) - circVolume_(i));
            
            //            file << std::setprecision(12) << std::fixed << circVolume_(i) << " " << fineVolume_(i) << " " << estPressure_(i) << " ";
            //        }
            //        file << std::endl;
            //        file.close();
            //        for(TInt i = 0; i < nC_; i++)
            //        {
            if (residual[i] > tolerance_) {
                DCCtrl::debug << "\nResidual not low enough...";
                if (i_ > maxIterations_) {
                    status_ = CBStatus::FAILED;
                    DCCtrl::debug << "FAILED";
                } else {
                    status_ = CBStatus::REPEAT;
                    DCCtrl::debug << "REPEAT";
                    break;
                }
            }
            DCCtrl::debug << "\nCavity " << i << ": " << residual[i] << " ml";
        }
        
        if (dampingActive_ && (status_ == CBStatus::SUCCESS) ) {
            status_ = CBStatus::REPEAT;
            dampingIteration_ = true;
        }
    }
} // CBCirculation::AnalyzeResults

void CBCirculation::StepBack() {
    if (!hasStarted_)
        return;
    
    // DCCtrl::print << "CBCirculation::StepBack()" << std::endl;
    estPressure_ = pressure_[(n_-1)%N_];
    stepBack_ = true;
    i_ = 0;
    j_ = 0;
}

void CBCirculation::Extrapolate(std::vector<VectorN<TFloat>> *in, VectorN<TFloat> *out) {
    // extrapolate using 4th order Adams-Bashforth
    *out = (*in)[(n_-1)%N_] + dt_ / 24.0 *
    (55.0 * ( (*in)[(n_-1)%N_] - (*in)[(n_-2)%N_]) / timeStep_[(n_-1)%N_]
     -59.0 * ( (*in)[(n_-2)%N_] - (*in)[(n_-3)%N_]) / timeStep_[(n_-2)%N_]
     +37.0 * ( (*in)[(n_-3)%N_] - (*in)[(n_-4)%N_]) / timeStep_[(n_-3)%N_]
     - 9.0 * ( (*in)[(n_-4)%N_] - (*in)[(n_-5)%N_]) / timeStep_[(n_-4)%N_]);
}

void CBCirculation::ExtrapolatePressure() {
    prevEstPressure_ = estPressure_;
    Extrapolate(&pressure_, &estPressure_);
    if (dampingActive_)
        estPressure_ = (estPressure_ - dampingFactor_ * pressure_[(n_-1)%N_]) / (1.0 - dampingFactor_);
    str1_ << "Extrapolating";
}

void CBCirculation::IncreasePressure() {
    prevEstPressure_ = estPressure_;
    estPressure_ += deltaP_;
    str1_ << "Increasing";
}

void CBCirculation::PerturbPressure(TInt i) {
    prevEstPressure_(i) = estPressure_(i);
    estPressure_(i) += dP_;
    str1_ << "Perturbing cavity " << i;
}

void CBCirculation::AssembleComplianceMatrices(TInt i) {
    TFloat dP = estPressure_(i) - prevEstPressure_(i);
    
    fineCompli_.SetCol(i, (fineVolume_-prevFineVolume_)/dP);
    circCompli_.SetCol(i, (circVolume_-prevCircVolume_)/dP);
}

void CBCirculation::EstimatePressureUsingSecant() {
    // get slope of secant (dPdR)
    TFloat *dPdR = new TFloat[nC_];
    
    for (TInt i = 0; i < nC_; i++) {
        TFloat dP = estPressure_(i) - prevEstPressure_(i);
        TFloat estFineCompli = (fineVolume_(i) - prevFineVolume_(i)) / dP;
        TFloat estCircCompli = (circVolume_(i) - prevCircVolume_(i)) / dP;
        dPdR[i] = 1.0 / (estFineCompli - estCircCompli);
    }
    
    // estimate pressure using secant method
    prevEstPressure_ = estPressure_;
    for (TInt i = 0; i < nC_; i++) {
        TFloat R = fineVolume_(i) - circVolume_(i);
        estPressure_(i) -= dPdR[i] * R;
    }
    delete[] dPdR;
    str1_ << "Secant";
}

void CBCirculation::EstimatePressureUsingQuasiNewton() {
    // estimating new system compliance Matrix using "Good Broyden Update"
    
    prevResidual_ = residual_;
    residual_ = fineVolume_ - circVolume_;
    
    if ((i_ > 1) || (n_ < 10) ) { // use previous inverse Jacobian for first iteration of each new timestep
        VectorN<TFloat> deltaF = residual_-prevResidual_;
        VectorN<TFloat> deltaX = estPressure_-prevEstPressure_;
        VectorN<TFloat> vec1 = sysElast_ * deltaF;
        MatrixN<TFloat> correctorMat = DyadicProduct((deltaX-vec1)/(deltaX*vec1), deltaX) * sysElast_;
        sysElast_ += correctorMat;
    }
    
    prevEstPressure_ = estPressure_;
    estPressure_ -= sysElast_ * residual_;
    
    str1_ << "Quasi-Newton";
    
    //    str1_ << std::endl << "System compliance matrix:" << std::endl << sysElast_.GetInverse();
}

void CBCirculation::EstimatePressureUsingNewton() {
    // perform one regular Newton step when finished assembling
    
    residual_ = fineVolume_ - circVolume_;
    
    sysCompli_ = fineCompli_ - circCompli_;
    sysElast_ = sysCompli_.GetInverse();
    
    prevEstPressure_ = estPressure_;
    estPressure_ -= sysElast_ * residual_;
    
    j_ = 0;
    
    str1_ << "Newton";
    
    //    str1_ << std::endl << "System compliance matrix:" << std::endl << sysCompli_;
}
