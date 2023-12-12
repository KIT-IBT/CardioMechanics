/*
 * File: CBCirculationCavity.h
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


#ifndef CB_CIRCULATION_CAVITY_H
#define CB_CIRCULATION_CAVITY_H

#include "CBElementCavity.h"

using namespace math_pack;

// TODO: What is the difference between CBCirculationCavity and CBCavity?

class CBCirculationCavity
{
public:
    CBCirculationCavity();
    ~CBCirculationCavity();
    /// inserts an element into the two internal target elements lists, after checking ignoredSurfaceIndices, ignoredElementIndices etc.
    void InsertNextElement(CBElementCavity* element);
    /// insert from allElementsVec, inserts only valid ones by respecting ignoredSurfaceIndices, ignoredEleIndices etc.
    void InsertElementsFromSurfaceIndex(std::vector<CBElementCavity*> allElementsVec, TInt surfaceIndex);
    /// set the specific index of this circulation cavity
    void SetIndex(TInt index);
    TInt GetIndex();
    /// return the volume stored in an internal variable usually CalcVolume has to be called first.
    TFloat GetVolume();
    /// Updates the internal volume variable by calculating the current volume.
    TFloat CalcVolume();
    TFloat CalcVolume(TFloat* referenceCoords);
    /// Check if the cavity surface is closed (true) or not closed (false).
    bool ClosedSurfaceCheck();
    void ApplyPressure(TFloat pressure);
    void ApplyPressureToJacobian(TFloat pressure);
    void SetAdapter(CBElementAdapter* adapter) { adapter_ = adapter; }
    CBElementAdapter* GetAdapter() { return adapter_; }
    
    void SetIgnoredSurfaceIndices(std::vector<TInt> ignoredSurfaceIndices);
    void SetIgnoredSurfaceElementsIndices(std::vector<TInt> ignoredSurfaceElementsIndices);
    void SetPokedNodeIndices(std::vector<TInt> pokedNodeIndices);
    
protected:
    
private:
    std::vector<CBElementCavity*>  targetElementsVolume_; //< all surface elements that are used to calculate the cavity volume, these should form a closed surface
    std::vector<CBElementCavity*>  targetElementsPressure_; //< all surface elements without 'ignoredSurfaceIndices' and without 'ignoredElementsIndices', usually a subset of targetElementsVolume_
    TInt                           index_ = 0;
    TFloat                         volume_ = 0;
    CBElementAdapter*              adapter_ = 0;
    
    bool initIgnoredSurfaceIndicesDone_ = false;
    bool initIgnoredSurfaceElementsIndicesDone_ = false;
    bool initPokedNodeIndicesDone_ = false;
    
    std::vector<TInt> ignoredSurfaceIndices_;
    std::vector<TInt> ignoredSurfaceElementsIndices_;
    std::vector<TInt> pokedNodeIndices_;
};

#endif
