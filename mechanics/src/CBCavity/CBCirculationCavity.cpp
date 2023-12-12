/*
 * File: CBCirculationCavity.cpp
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


#include "DCCtrl.h"

#include "CBCirculationCavity.h"
#include "CBElementAdapter.h"
#include "Triangle.h"

#include <algorithm>

CBCirculationCavity::CBCirculationCavity() {}

CBCirculationCavity::~CBCirculationCavity() {
    for (auto &element : targetElementsVolume_)
        delete element;
    
    // this covers all elements, overset of targetElementsPressure_
}

void CBCirculationCavity::InsertNextElement(CBElementCavity *element) {
    // volume should be computed on all elements, pressure should be only applied on a certain subset
    
    targetElementsVolume_.push_back(element);
    
    // ignore certain surface indices / GetMaterialIndex()
    if (!initIgnoredSurfaceIndicesDone_) {
        throw std::runtime_error(
                                 "CBCirculationCavity::InsertNextElement(): ignoredSurfaceIndices_ was not initialized! You have to run 'SetIgnoredSurfaceIndices()' first.");
    }
    if (std::find(ignoredSurfaceIndices_.begin(), ignoredSurfaceIndices_.end(),
                  element->GetMaterialIndex()) != ignoredSurfaceIndices_.end() )
        return void();
    
    // ignore specific element indices
    if (!initIgnoredSurfaceElementsIndicesDone_) {
        throw std::runtime_error(
                                 "CBCirculationCavity::InsertNextElement(): ignoredSurfaceElementsIndices_ was not initialized! You have to run 'SetIgnoredSurfaceElementsIndices()' first.");
    }
    if (std::find(ignoredSurfaceElementsIndices_.begin(), ignoredSurfaceElementsIndices_.end(),
                  element->GetIndex()) != ignoredSurfaceElementsIndices_.end() )
        return void();
    
    targetElementsPressure_.push_back(element);
}

void CBCirculationCavity::InsertElementsFromSurfaceIndex(std::vector<CBElementCavity *> allElementsVec,
                                                         TInt surfaceIndex) {
    for (auto ele : allElementsVec)
        if (ele->GetSurfaceIndex() == surfaceIndex)
            InsertNextElement(ele);
}

void CBCirculationCavity::SetIndex(TInt index) {
    index_ = index;
}

TInt CBCirculationCavity::GetIndex() {
    return index_;
}

TFloat CBCirculationCavity::GetVolume() {
    return volume_;
}

TFloat CBCirculationCavity::CalcVolume() {
    TFloat referenceCoords[3] = {0.0, 0.0, 0.0};
    
    return CalcVolume(referenceCoords);
}

TFloat CBCirculationCavity::CalcVolume(TFloat *referenceCoords) {
    double localVolume = 0;
    double volume      = 0;
    
    for (auto &element : targetElementsVolume_)
        localVolume += element->CalcContributionToVolume(referenceCoords);
    
    MPI_Allreduce(&localVolume, &volume, 1, MPI_DOUBLE, MPI_SUM, Petsc::Comm());
    
    volume_ = volume;
    return volume_;
}

bool CBCirculationCavity::ClosedSurfaceCheck() {
    TFloat origin0[3] = {0., 0., 0.};
    TFloat origin1[3] = {1., 1., 1.};
    TFloat volume0 = 1e6 * CalcVolume(origin0); // volume in ml
    TFloat volume1 = 1e6 * CalcVolume(origin1);
    
    if (std::abs(volume1 - volume0) > 1e-10)
        return false;
    
    return true;
}

void CBCirculationCavity::ApplyPressure(TFloat pressure) {
    for (auto &element : targetElementsPressure_) {
        // t3,t6 without poked -> works by default
        // t3 with poked -> works with special code
        // t6 with poked -> does not work
        
        // t3,t6 without poked -> works by default
        if (pokedNodeIndices_.size() == 0) {
            // TODO: In a clean implementation of PokedPoints, the following is the only line that is actually needed in the whole function ApplyPressure().
            element->ApplyPressure(pressure);
        } else {
            if (element->GetNumberOfNodesIndices() == 3) {
                // t3 with poked -> works with special code
                
                //////////
                
                // What is the code below doing?
                // -> Probably copied from t3 ApplyPressure, but extended to respect (i.e. ignore) PokedNodeIndices_
                // Does not work correctly with T6 elements!
                
                
                if ((pokedNodeIndices_.size() != 0) && (element->GetNumberOfNodesIndices() == 3)) {
                    Triangle<TFloat> triangle = element->GetTriangle();
                    Vector3<TFloat> totalForce = -pressure *triangle.GetArea() * triangle.GetNormalVector();
                    
                    TInt nodesCoordsIndices[9];
                    TInt pressureMask[3] = {1, 1, 1};
                    for (unsigned int i = 0; i < 3; i++) {
                        nodesCoordsIndices[3 * i]     = 3 * element->GetNodeIndex(i);
                        nodesCoordsIndices[3 * i + 1] = 3 * element->GetNodeIndex(i) + 1;
                        nodesCoordsIndices[3 * i + 2] = 3 * element->GetNodeIndex(i) + 2;
                        
                        if (std::find(pokedNodeIndices_.begin(), pokedNodeIndices_.end(),
                                      element->GetNodeIndex(i)) != pokedNodeIndices_.end())
                            pressureMask[i] = 0;
                    }
                    TInt pressureMaskSum = pressureMask[0]+pressureMask[1]+pressureMask[2];
                    TFloat forceFactor = 1.0;
                    if (pressureMaskSum != 0)
                        forceFactor /= pressureMaskSum;
                    
                    bool bc[9];
                    element->GetAdapter()->GetNodesComponentsBoundaryConditions(9, nodesCoordsIndices, bc);
                    
                    TFloat forces[9];
                    for (unsigned int i = 0; i < 3; i++) {
                        for (unsigned int j = 0; j < 3; j++) {
                            if ((bc[3*i+j] != 0) || (pressureMask[i] == 0))
                                forces[3 * i + j] = 0;
                            else
                                forces[3 * i + j] = forceFactor * totalForce.Get(j);
                        }
                    }
                    element->GetAdapter()->AddNodalForcesComponents(9, nodesCoordsIndices, forces);
                }
                
                //////////
            } else {
                // t6 with poked -> does not work
                throw std::runtime_error(
                                         "CBCirculationCavity::ApplyPressure() ERROR PokedPoints is not supported for t6 elements (Plugins.PokedPoints.NodeIndices) ! You can use t3 elements with PokedPoints, t3/t6 without PokedPoints, or implement what is missing.");
            }
        }
    }
} // CBCirculationCavity::ApplyPressure

void CBCirculationCavity::ApplyPressureToJacobian(TFloat pressure) {
    for (auto &element : targetElementsPressure_)
        element->CalcPressureJacobian(pressure);
}

void CBCirculationCavity::SetIgnoredSurfaceIndices(std::vector<TInt> ignoredSurfaceIndices) {
    ignoredSurfaceIndices_ = ignoredSurfaceIndices;
    initIgnoredSurfaceIndicesDone_ = true;
}

void CBCirculationCavity::SetIgnoredSurfaceElementsIndices(std::vector<TInt> ignoredSurfaceElementsIndices) {
    ignoredSurfaceElementsIndices_ = ignoredSurfaceElementsIndices;
    initIgnoredSurfaceElementsIndicesDone_ = true;
}

void CBCirculationCavity::SetPokedNodeIndices(std::vector<TInt>
                                              pokedNodeIndices) {
    pokedNodeIndices_ = pokedNodeIndices;
    initPokedNodeIndicesDone_ = true;
}
