/*
 * File: CBTensionEstimator.cpp
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


#include "CBTensionEstimator.h"

CBTensionEstimator::CBTensionEstimator(CBElementSolid *ele, ParameterMap *) {
    e_ = ele;
    assert(e_ != nullptr);
    Tmax_ = ele->GetMaterial()->GetProperties()->tensionMax_;
    activeTension_ = 0.0;
}

TFloat CBTensionEstimator::GetActiveTension() {
    return activeTension_ * Tmax_;
}

double CBTensionEstimator::CalcActiveTension(const math_pack::Matrix3<double> &, const double) {
    return activeTension_ * Tmax_;
}

math_pack::Matrix3<double> CBTensionEstimator::CalcActiveStress(const math_pack::Matrix3<double> &deformation, const double time) {
    // The inverse problem estimates the scalar fiber tension against target surfaces that were
    // generated with a pure fiber-direction active stress. Place the tension in the (0,0)/fiber
    // component only, without the deformation-dependent 1/sqrt(I4) scaling of the base class, so the
    // estimator recovers the same active stress that produced the targets.
    return CalcActiveTension(deformation, time) * math_pack::Matrix3<double> {1, 0, 0,  0, 0, 0,  0, 0, 0};
}

CBStatus CBTensionEstimator::SetActiveTensionAtQuadraturePoint(int, TFloat activeTension) {
    activeTension_ = activeTension;
    return CBStatus::SUCCESS;
}
