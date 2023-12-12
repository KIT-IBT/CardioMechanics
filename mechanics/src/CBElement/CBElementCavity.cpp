/*
 * File: CBElementCavity.cpp
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


#include "CBElementCavity.h"

TFloat CBElementCavity::CalcContributionToVolume()
{
    TFloat referenceCoords[3] = {0.0, 0.0, 0.0};
    return CalcContributionToVolume(referenceCoords);
}

Triangle<TFloat> CBElementCavity::GetTriangle()
{
    throw std::runtime_error("CBElementCavity::GetTriangle is not implemented for the requested type!");
    /// note: code is only implemented in CBElementSurfaceT3 [lb451]
}
