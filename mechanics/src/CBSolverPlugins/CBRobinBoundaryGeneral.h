/*
 * File: CBRobinBoundaryGeneral.h
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


#ifndef CB_ROBIN_BOUNDARY_GENERAL
#define CB_ROBIN_BOUNDARY_GENERAL

#include <map>
#include <memory>

#include "CBSolverPlugin.h"
#include "CBElementContactRobin.h"

using namespace math_pack;

class CBRobinBoundaryGeneral : public CBSolverPlugin {
public:
    CBRobinBoundaryGeneral();
    ~CBRobinBoundaryGeneral();
    
    std::string GetName() override { return "RobinBoundaryGeneral"; }
    
    /// fundamental CBSolverPlugin functions
    void Init() override;
    void Apply(PetscScalar time) override;
    void ApplyToNodalForces() override;
    void ApplyToNodalForcesJacobian() override;
    void StepBack() override;
    void Export(TFloat time) override;
    void WriteToFile(TFloat time) override;
    void Prepare() override;
    
    CBStatus GetStatus() override {return status_;}
    
protected:
private:
    TFloat dt_, prevDt_;
    typedef CBSolverPlugin Base;
    std::vector<CBElementContactRobin *> contactSurfaceElements_;
    std::vector<Vector3<TFloat>> initialPos_;
    std::vector<Vector3<TFloat>> displacement_;
    std::vector<Vector3<TFloat>> prevDisplacement_;
    std::vector<Vector3<TFloat>> velocity_;
    
    CBStatus status_ = CBStatus::DACCORD;
    
    bool isFirstStep_ = true;
    bool stepBack_ = false;
    bool hasStarted_ = false;
    
    /// functions
    void InitContactSurfaces();
    void CalcForceContributionOfElement(Vector3<TFloat> u, Vector3<TFloat> v, CBElementContactRobin *triangle,
                                        TFloat *nodalForces);
    
    /// xml parameters
    /// general options
    TFloat startTime_;
    bool export_;
    
    /// boundary model specific options
    TInt normalVectorSign_ = 1;
    PetscScalar alpha_;
    PetscScalar beta_;
    PetscInt surfaceIndex_;
    std::string filename_;
    std::ofstream file_;
    
    
    // ----- Values needed for export -----
    
    std::vector<Vector3<TFloat>> ContactForces_;
};
#endif // ifndef CB_ROBIN_BOUNDARY_GENERAL
