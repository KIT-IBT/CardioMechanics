//
//  CBSolverPluginCavities.cpp
//  CardioMechanics
//
//  Created by Steffen Schuler on 01.03.16.
//
//

#include "CBSolverPluginCavities.h"
#include "CBSolver.h"
#include "CBElementCavity.h"
// #include "CBElementSurfaceT3.h"
// #include "CBElementSurfaceT6.h"


/// this allows to initialize cavities_ from a differently named derived SolverPlugin
/// the whole init stuff should go here. PluginRootKey allows to initialize cavities_ from a differently named derived SolverPlugin.
void CBSolverPluginCavities::InitCavities(std::string pluginRootKey)
{
    if (pluginRootKey == "")
        pluginRootKey = "Mesh";
    
    ignoredSurfaceIndices_         = parameters_->GetArray<TInt>(pluginRootKey + ".IgnoredSurfaceIndices", {});
    ignoredSurfaceElementsIndices_ = parameters_->GetArray<TInt>(pluginRootKey + ".IgnoredSurfaceElementsIndices", {});
    pokedNodeIndices_              = parameters_->GetArray<TInt>("Plugins.PokedPoints.NodeIndices", {});
    
    CreateCavities();
    SyncCavities();
}

/// creates cavities for each surface id found in the procesors local elements
void CBSolverPluginCavities::CreateCavities()
{
    DCCtrl::print << "create cavities..." << std::endl;
    // 1) create all (empty) cavities with correct cavity index,
    //    create a vector containing all elements of all cavities
    // 2) set necessary ignoredIndices, PokedNodes, ignoredMaterials etc.
    // 3) pass all elements to each, they decide on their own which elements to take
    
    std::vector<CBElementCavity*> allCavityElements;
    std::unordered_set<TInt> possibleCavityIds;
    for (auto eleIt : adapter_->GetElementVector())
    {
        // Iterate over all CBElements and check if it is of type CBElementCavity.
        // If so, get its surface index and add it to the cavity that is associated with the same surface index (this association is mapped in cavities_).
        CBElementCavity* element = dynamic_cast<CBElementCavity*>(eleIt);
        if (element != 0)
        {
            TInt index = element->GetSurfaceIndex();
            if (index <= 0)
                throw std::runtime_error("void CBSolverPluginCavities::InitCavities(): Surface indices must be > 0");
            
            allCavityElements.push_back(element);
            possibleCavityIds.insert(element->GetSurfaceIndex());
            
        }
    }
    
    // create empty cavities for all possible surface cavities
    for (auto index : possibleCavityIds)
    {
        DCCtrl::print << "possible index: " << index << std::endl;
        CBCirculationCavity* cavity = new CBCirculationCavity;
        cavity->SetAdapter(GetAdapter());
        cavity->SetIndex(index);
        cavities_.insert(std::pair<TInt, CBCirculationCavity*>(index, cavity));
    }
    
    // For each cavity, set the surface elements to ignore when applying pressure
    // and fill each cavity with the corresponding elements
    for (auto& cavity : cavities_)
    {
        cavity.second->SetIgnoredSurfaceIndices(ignoredSurfaceIndices_);
        cavity.second->SetIgnoredSurfaceElementsIndices(ignoredSurfaceElementsIndices_);
        cavity.second->SetPokedNodeIndices(pokedNodeIndices_);
        
        // insertFromAllElements automagically respect ignored surfaces, ignored indices etc.
        cavity.second->InsertElementsFromSurfaceIndex(allCavityElements, cavity.second->GetIndex());
    }
    
    DCCtrl::print << "create cavities... done" << std::endl;
}

void CBSolverPluginCavities::SyncCavities()
{
    // Syncing cavities is necessary since otherwise we will get a problem when we calculate whole volume of the cavity. In this case
    // cavities_ might be empty for one CPU and MPI_Allreduce produces a dead lock
    
    PetscInt* localCavities  = new PetscInt[DCCtrl::GetNumberOfProcesses()];
    PetscInt* globalCavities = new PetscInt[DCCtrl::GetNumberOfProcesses()];
    
    for (PetscInt i = 0; i < DCCtrl::GetNumberOfProcesses(); i++)
    {
        if (i == DCCtrl::GetProcessID())
            localCavities[i] = cavities_.size();
        else
            localCavities[i] = 0;
    }
    
    MPI_Allreduce(localCavities, globalCavities, DCCtrl::GetNumberOfProcesses(), MPI_INT, MPI_SUM, Petsc::Comm());
    
    PetscInt numIndices = 0;
    for (PetscInt i = 0; i < DCCtrl::GetNumberOfProcesses(); i++)
        numIndices += globalCavities[i];
    
    PetscInt* localCavitiesIndices = new PetscInt[numIndices];
    for (PetscInt i = 0; i < numIndices; i++)
        localCavitiesIndices[i] = 0;
    
    PetscInt* globalCavitiesIndices = new PetscInt[numIndices];
    
    PetscInt  pos = 0;
    for (PetscInt i = 0; i < DCCtrl::GetProcessID(); i++)
        pos += globalCavities[i];
    
    for (std::map<PetscInt, CBCirculationCavity*>::iterator it = cavities_.begin(); it != cavities_.end(); it++)
    {
        localCavitiesIndices[pos] = it->second->GetIndex();
        pos++;
    }
    
    MPI_Allreduce(localCavitiesIndices, globalCavitiesIndices, numIndices, MPI_INT, MPI_SUM, Petsc::Comm());
    
    for (PetscInt i = 0; i < numIndices; i++)
    {
        if (cavities_.find(globalCavitiesIndices[i]) == cavities_.end())
        {
            CBCirculationCavity* cavity = new CBCirculationCavity;
            cavity->SetAdapter(GetAdapter());
            cavity->SetIndex(globalCavitiesIndices[i]);
            
            cavities_.insert(std::pair<PetscInt, CBCirculationCavity*>(globalCavitiesIndices[i], cavity));
        }
    }
    
    for (auto it : cavities_)
        DCCtrl::print << "cavity: " << it.first << std::endl;
    
    delete[] localCavities;
    delete[] globalCavities;
    
    delete[] localCavitiesIndices;
    delete[] globalCavitiesIndices;
}

/// returns the map assigning surface indices to cavities
// std::map<TInt, CBCirculationCavity*>& CBSolverPluginCavities::GetCavities() 
// { 
//   return cavities_; 
// };

/// returns a single cavity corresponding to the surface index 'cavityIndex'
// CBCirculationCavity* CBSolverPluginCavities::GetCavity(int cavityIndex) 
// { 
//   return cavities_.at(cavityIndex); 
// };
