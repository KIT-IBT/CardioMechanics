/*
 * File: CBSolverActiveStressEstimator.h
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


#pragma once

#include <fstream>
#include <set>
#include <vector>

#include "CBSolver.h"
#include "CBSolverEquilibrium.h"
#include "CBContactHandling.h"

/// Inverse problem: estimate the scalar active tension in every element at each timestep so that the
/// deformed surface matches a prescribed target surface. A Gauss-Newton / Tikhonov iteration wraps the
/// static equilibrium solve; the target-vs-current surface distance is provided by the ContactHandling
/// plugin, and the estimated tension is injected into each element's TensionEstimator tension model.
class CBSolverActiveStressEstimator : public CBSolverEquilibrium {
public:
    CBSolverActiveStressEstimator() : CBSolverEquilibrium() {}
    std::string GetType() override { return "Active Stress Estimator (Static)"; }
    void Init(ParameterMap *parameters, CBModel *model) override;

protected:
    CBStatus SolverStep(PetscScalar time, bool forceJacobianAndDampingRecalculation = false) override;

private:
    void UpdateMasterNodesOfInterest();
    CBStatus EstimatorStep(PetscScalar time, int step);
    void GenerateElementLaplacian();
    void WriteToFile(PetscScalar time, int iter, TFloat norm);
    void WriteLcurveStuffToFile(PetscScalar time, int iter, TFloat resNorm, TFloat regNorm);
    void WriteHeaderToFile(std::string filename);

    CBContactHandling *contact_ = nullptr;

    IS   masterNodesIndices_    = nullptr;
    IS   nodesOfInterestIndices_ = nullptr;
    IS   elementsOfInterestIndices_ = nullptr;
    TInt numMasterNodesIndices_    = 0;
    TInt numNodesOfInterestIndices_ = 0;
    TInt numElementOfInterestIndices_ = 0;

    PetscInt *elementsOfInterestMapping_ = nullptr;
    PetscInt *masterNodesIndicesNodesOfInterestMapping_ = nullptr;
    PetscInt *masterNodesIndicesMapping_ = nullptr;
    PetscInt *nim_ = nullptr;

    std::vector<TInt> mat_;
    std::set<TInt>    mn_;

    // Nonzeros per row of dfdtau: one entry per solid element adjacent to the row's node.
    std::vector<PetscInt> dfdtauNnz_;

    Mat lTl_             = nullptr;
    Mat elementLaplacian_ = nullptr;
    Mat inv_             = nullptr;

    Vec dist_      = nullptr;
    Vec tmpNodes_  = nullptr;
    Vec ti_        = nullptr;  // T_i
    Vec ti1_       = nullptr;  // T_{i-1}
    Vec ti2_       = nullptr;  // T_{i-2}

    TFloat l1_ = 0, l2_ = 0, l3_ = 0, l4_ = 0, l5_ = 0, l6_ = 0, l7_ = 0;

    double lastTime_ = 0;
    TFloat lastEstimatorTimeStep_ = -1;

    std::ofstream file_;
    std::ofstream fileLcurve_;
    std::string   filename_;
    std::string   filenameLcurve_;
    bool          headerWritten_ = false;
    bool          headerLcurveWritten_ = false;
    TFloat        resNorm_ = -1;
    TFloat        regTermNorm_ = -1;
};
