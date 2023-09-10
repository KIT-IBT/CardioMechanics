/*
 *  CBElementCavity.cpp
 *  CardioMechanics
 *
 *  Created by Thomas Fritz on 14.09.11.
 *  Copyright 2011 Institute for Biomedical Engineering, Karlsruhe Institute of Technology (KIT). All rights reserved.
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
