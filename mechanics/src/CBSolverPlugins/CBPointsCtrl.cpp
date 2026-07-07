/*
 * File: CBPointsCtrl.cpp
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


#include <fstream>
#include <sstream>
#include <vector>

#include "CBPointsCtrl.h"
#include "CBSolver.h"

void CBPointsCtrl::Init() {
    std::string filename = parameters_->Get<std::string>("Plugins.PointsCtrl.Coordinates.File");

    pointsCoords_ = new CBDataFromFile();
    pointsCoords_->Init(filename);

    filename = parameters_->Get<std::string>("Plugins.PointsCtrl.Points.File", "");
    bool *bc = GetAdapter()->GetSolver()->GetNodesComponentsBoundaryConditionsGlobal();

    if (filename != "") {
        std::ifstream controllNodesFile(filename.c_str());
        if (!controllNodesFile.good())
            throw std::runtime_error("CBPointsCtrl::Init(): File " + filename + " does not exist");

        std::string str;
        while (getline(controllNodesFile, str)) {
            std::stringstream ss(str);
            int n;
            ss >> n;
            n--;
            bc[3*n]   = true;
            bc[3*n+1] = true;
            bc[3*n+2] = true;
            n = GetAdapter()->GetSolver()->GetModel()->GetForwardMapping(n);

            if ((n >= GetAdapter()->GetSolver()->GetLocalNodesFrom()) &&
                (n <= GetAdapter()->GetSolver()->GetLocalNodesTo()))
                controllNodesLocalIndices_.push_back(n - GetAdapter()->GetSolver()->GetLocalNodesFrom());
        }
    } else {
        std::vector<TInt> fromTo = parameters_->GetArray<TInt>("Plugins.PointsCtrl.Points.FromTo");
        if (fromTo.size() != 2)
            throw std::runtime_error("CBPointsCtrl::Init(): Plugins.PointsCtrl.Points.FromTo may only contain two values");

        for (int i = fromTo[0]-1; i <= fromTo[1]-1; i++) {
            bc[3*i]   = true;
            bc[3*i+1] = true;
            bc[3*i+2] = true;

            TInt n = GetAdapter()->GetSolver()->GetModel()->GetForwardMapping(i);

            if ((n >= GetAdapter()->GetSolver()->GetLocalNodesFrom()) &&
                (n <= GetAdapter()->GetSolver()->GetLocalNodesTo()))
                controllNodesLocalIndices_.push_back(n - GetAdapter()->GetSolver()->GetLocalNodesFrom());
        }
    }

    GetAdapter()->LinkNodesComponentsBoundaryConditionsGlobal(bc);
}

void CBPointsCtrl::Apply(TFloat time) {
    for (size_t i = 0; i < controllNodesLocalIndices_.size(); i++) {
        TFloat newCoords[3] = {pointsCoords_->Get(time, 3*i),
                               pointsCoords_->Get(time, 3*i+1),
                               pointsCoords_->Get(time, 3*i+2)};
        TInt newCoordsIndices[3] = {3*controllNodesLocalIndices_[i],
                                    3*controllNodesLocalIndices_[i]+1,
                                    3*controllNodesLocalIndices_[i]+2};

        // A target of exactly (0,0,0) marks a node without a valid target this step; skip it.
        if (!(newCoords[0] == 0.0 && newCoords[1] == 0.0 && newCoords[2] == 0.0))
            GetAdapter()->SetNodesCoords(3, newCoordsIndices, newCoords);
    }
}
