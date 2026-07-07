/*
 * File: CBPointsCtrl.h
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

#include "CBSolverPlugin.h"
#include "CBDataFromFile.h"

/// Prescribes the coordinates of a set of control nodes over time from a
/// coordinate-list file. Used to drive the target surface of the inverse
/// active-stress estimator: the listed nodes are pinned (Dirichlet) and moved
/// to their target positions each timestep.
class CBPointsCtrl : public CBSolverPlugin {
public:
    ~CBPointsCtrl() { if (pointsCoords_) delete pointsCoords_; }
    void Apply(TFloat time) override;
    std::string GetName() override { return std::string("PointsControl"); }
    void Init() override;

protected:
    CBDataFromFile   *pointsCoords_ = nullptr;
    std::vector<TInt> controllNodesLocalIndices_;
};
