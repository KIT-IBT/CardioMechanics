/*
 * File: CBCirculation.h
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


#ifndef CB_CIRCULATION
#define CB_CIRCULATION

#include <map>

#include "CBSolverPlugin.h"
#include "CBCirculationCavity.h"
#include "CBCircModel.h"
#include "MatrixN.h"

#include "CBSolverPluginCavities.h"

class CBCirculation : public CBSolverPluginCavities {
public:
    CBCirculation();
    ~CBCirculation();
    
    std::string GetName() override { return "Circulation"; }
    
    void Init() override;
    void Prepare() override;
    double GetPreparationProgress() override;
    
    CBStatus GetStatus() override { return status_; }
    
    bool ExitCheck(TFloat time) override;
    void WriteToFile(PetscScalar time) override;
    void Apply(PetscScalar time) override;
    void ApplyToNodalForces() override;
    void ApplyToNodalForcesJacobian() override;
    
    bool WantsToAnalyzeResults() override { return true; }
    
    void AnalyzeResults() override;
    void StepBack() override;
    
protected:
private:
    // void InitCavities(); // inherited from CBApplyPressureFromFunction
    // void SyncCavities(); // inerited from CBApplyPressureFromFunction
    void InitCircModels();
    void GetCircVolume(VectorN<TFloat> *pressure, VectorN<TFloat> *volume);
    void Extrapolate(std::vector<VectorN<TFloat>> *in, VectorN<TFloat> *out);
    void ExtrapolatePressure();
    void IncreasePressure();
    void PerturbPressure(TInt i);
    void AssembleComplianceMatrices(TInt i);
    void EstimatePressureUsingSecant();
    void EstimatePressureUsingQuasiNewton();
    void EstimatePressureUsingNewton();
    
    typedef CBSolverPlugin Base;
    
    // ----------
    
    std::ostringstream str1_, str2_;
    
    std::vector<CBCircModelBase *> circModels_;
    std::vector<TInt> cavitySurfaceIndices_;
    TInt nC_ = 0;               // total number of cavities, must be initialized to zero!
    
    CBStatus status_ = CBStatus::DACCORD;
    
    TFloat uV_, uP_;
    TFloat maxP_, deltaP_, dP_;
    
    bool stepBack_ = false;
    
    bool preloadingFinished_ = false;
    bool lastPreloadingStep1_ = false, lastPreloadingStep2_ = false;
    TFloat preloadingTime1_, preloadingTime2_;
    VectorN<TFloat> preloadingPressure1_, preloadingPressure2_;
    
    TFloat dt_, prevDt_;
    const static TInt N_ = 5;   // number of past values stored for extrapolation
    TFloat timeStep_[N_];
    std::vector<VectorN<TFloat>> pressure_;
    VectorN<TFloat> estPressure_, prevEstPressure_;
    VectorN<TFloat> fineVolume_, prevFineVolume_;
    VectorN<TFloat> circVolume_, prevCircVolume_;
    
    VectorN<TFloat> residual_, prevResidual_;
    
    MatrixN<TFloat> fineCompli_;
    MatrixN<TFloat> circCompli_;
    MatrixN<TFloat> sysCompli_;
    MatrixN<TFloat> sysElast_;
    
    TFloat tolerance_;
    TInt maxIterations_;
    TInt nNonNewton_;           // perform perturbations and a regular Newton step when the past nNonNewton_ iterations with the secant or quasi-Newton method did not achieve convergence
    TInt nSecant_;              // number of iterations to perform using the secant method before perturbation
    bool perturbation_;         // Shall the system compliance matrix initially be assembled by perturbing the individual cavity pressures at each new timestep?
    unsigned long int n_ = 1;   // time step counter
    unsigned int i_ = 0;        // iterations counter
    unsigned int j_ = 0;        // non-newton-iterations counter
    unsigned long int iSum_ = 0;
    
    std::string itsFilename_;
    std::ofstream itsFile_;
    std::string filename_;
    
    TFloat startTime_;
    bool hasStarted_ = false;
    
    bool dampingActive_;
    TFloat dampingFactor_;
    TFloat dampingDeclineFactor_;
    bool dampingIteration_ = false;
    
    // ----------
};

#endif // ifndef CB_CIRCULATION
