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

CBStatus CBTensionEstimator::SetActiveTensionAtQuadraturePoint(int, TFloat activeTension) {
    activeTension_ = activeTension;
    return CBStatus::SUCCESS;
}
