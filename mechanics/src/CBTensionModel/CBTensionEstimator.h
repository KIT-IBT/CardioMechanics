/*
 * File: CBTensionEstimator.h
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

#include "CBTensionModel.h"
#include "ParameterMap.h"

/// Mailbox tension model for the inverse problem. The active-stress estimator
/// solver writes the estimated fiber tension into each element via
/// SetActiveTensionAtQuadraturePoint; the element stress assembly and the
/// exporter read it back through CalcActiveTension / CalcActiveStress.
class CBTensionEstimator : public CBTensionModel {
public:
    CBTensionEstimator(CBElementSolid *ele, ParameterMap *parameters);

    double CalcActiveTension(const math_pack::Matrix3<double> &deformation, const double time) override;
    TFloat GetActiveTension() override;
    CBStatus SetActiveTensionAtQuadraturePoint(int indexQP, TFloat activeTension) override;

protected:
    CBElementSolid *e_;
    TFloat activeTension_ = 0.0;
};
