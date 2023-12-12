/*
 * File: CBLoadUnloadedState.cpp
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
#include "DCCtrlPETSc.h"

#include "CBSolver.h"
#include <vector>
#include "petscvec.h"

#include "CBLoadUnloadedState.h"

void CBLoadUnloadedState::Init() {
    DCCtrl::debug << "\t--- --- Load Unloaded State --- ---" << std::endl;
    
    settleDown_ = parameters_->Get<bool>("Plugins.LoadUnloadedState.SettleDown", false);
    
    InitPetscVectors();
    
    adapter_->GetSolver()->GetNodeCoordinates(originalCoords_);
    
    if (DCCtrl::IsProcessZero()) {
        std::string filename = parameters_->Get<std::string>("Plugins.LoadUnloadedState.unloadedNodes");
        LoadNodes(filename); // loads into coordsSeq_
    }
    
    // Change to unloaded configuration
    DCCtrl::debug << "\tChange to unloaded configuration...";
    if (DCCtrl::IsParallel()) {
        // Scatter from zero to all
        VecScatterBegin(coordsScatter_, coordsSeq_, unloadedCoords_, INSERT_VALUES, SCATTER_REVERSE);
        VecScatterEnd(coordsScatter_, coordsSeq_, unloadedCoords_, INSERT_VALUES, SCATTER_REVERSE);
        
        adapter_->GetSolver()->SetNodeCoordinates(unloadedCoords_);
        adapter_->GetSolver()->SetRefNodeCoordinates(unloadedCoords_);
    } else {
        adapter_->GetSolver()->SetNodeCoordinates(coordsSeq_);
        adapter_->GetSolver()->SetRefNodeCoordinates(coordsSeq_);
    }
    
    DCCtrl::debug << "DONE." << std::endl;
    DCCtrl::debug << "\tRelax elements and bases...";
    adapter_->GetSolver()->RelaxElementsAndBases();
    DCCtrl::debug << "DONE." << std::endl;
    DCCtrl::debug << "\tSet zero velocity and accelleration...";
    adapter_->GetSolver()->SetZeroVelocityAndAcceleration();
    VecSet(coordsSeq_, 0.0);
    DCCtrl::debug << "DONE." << std::endl;
    
    if (DCCtrl::IsProcessZero()) {
        std::string filename = parameters_->Get<std::string>("Plugins.LoadUnloadedState.inflatedNodes");
        LoadNodes(filename);
    }
    
    // Change to inflated configuration
    DCCtrl::debug << "\tChange to inflated configuration...";
    if (DCCtrl::IsParallel()) {
        // Scatter from zero to all
        VecScatterBegin(coordsScatter_, coordsSeq_, inflatedCoords_, INSERT_VALUES, SCATTER_REVERSE);
        VecScatterEnd(coordsScatter_, coordsSeq_, inflatedCoords_, INSERT_VALUES, SCATTER_REVERSE);
        
        adapter_->GetSolver()->SetNodeCoordinates(inflatedCoords_);
    } else {
        adapter_->GetSolver()->SetNodeCoordinates(coordsSeq_);
    }
    
    DCCtrl::debug << "DONE." << std::endl;
    DCCtrl::debug << "\tSet zero velocity and accelleration...";
    adapter_->GetSolver()->SetZeroVelocityAndAcceleration();
    DCCtrl::debug << "DONE." << std::endl;
    DCCtrl::debug << "\tSet zero cumulative displacement...";
    adapter_->GetSolver()->SetZeroDisplacement();
    DCCtrl::debug << "DONE." << std::endl;
    
    // Free memory
    DeInitPetscVectors();
    
    DCCtrl::debug << "\t--- --- Load Unloaded State --- ---" << std::endl;
} // CBLoadUnloadedState::Init

void CBLoadUnloadedState::Apply(TFloat time) {
    if (isFirstStep_ && settleDown_) {
        DCCtrl::debug << "\t--- --- Load Unloaded State --- ---" << std::endl;
        DCCtrl::debug << "\tSet zero velocity and accelleration...";
        adapter_->GetSolver()->SetZeroVelocityAndAcceleration();
        isFirstStep_ = false;
        DCCtrl::debug << "DONE." << std::endl;
        DCCtrl::debug << "\tSet zero cumulative displacement...";
        adapter_->GetSolver()->SetZeroDisplacement();
        DCCtrl::debug << "DONE." << std::endl;
        DCCtrl::debug << "\t--- --- Load Unloaded State --- ---" << std::endl;
    }
}

void CBLoadUnloadedState::InitPetscVectors() {
    VecCreate(PETSC_COMM_WORLD, &originalCoords_);
    VecSetSizes(originalCoords_, 3*adapter_->GetSolver()->GetNumberOfLocalNodes(), PETSC_DETERMINE);
    VecSetFromOptions(originalCoords_);
    VecZeroEntries(originalCoords_);
    
    VecDuplicate(originalCoords_, &unloadedCoords_);
    VecZeroEntries(unloadedCoords_);
    VecDuplicate(originalCoords_, &inflatedCoords_);
    VecZeroEntries(inflatedCoords_);
    
    VecScatterCreateToZero(originalCoords_, &coordsScatter_, &coordsSeq_);
}

void CBLoadUnloadedState::DeInitPetscVectors() {
    VecDestroy(&originalCoords_);
    VecScatterDestroy(&coordsScatter_);
    VecDestroy(&coordsSeq_);
    VecDestroy(&unloadedCoords_);
    VecDestroy(&inflatedCoords_);
}

// ToDo: Code duplication. Better use CBModelLoaderTetgen::LoadNodes
void CBLoadUnloadedState::LoadNodes(std::string filename) {
    // Format Mesh.node:
    // <amount of nodes> <Dimension (must be 3)> <Boundary Conditions> (0: no boundary condition, 1: Fixation)> <BoundaryMarker>
    // <Node Index> <x> <y> <z> <Boundary Condition>
    
    std::ifstream file(filename.c_str());
    
    if (!file.good()) {
        file.close();
        throw std::runtime_error("CBLoadUnloadedState::LoadNodes: File " + filename + " is not accessable!");
    } else {
        TInt               numNodes;
        TInt               currentNode;
        unsigned int       dim;
        std::string        line;
        std::istringstream linestream;
        
        getline(file, line);
        linestream.str(line);
        linestream >> numNodes;
        linestream >> dim;
        
        if (numNodes != adapter_->GetSolver()->GetNumberOfNodes()) {
            throw std::runtime_error(
                                     "CBLoadUnloadedState::LoadNodes: Number of nodes is nodes file != number of nodes in loaded model.");
        }
        if (dim != 3)
            throw std::runtime_error("CBLoadUnloadedState::LoadNodes: Only for 3D model !");
        
        std::string currentNodeAsString;
        TFloat x, y, z;
        int i = 0;
        while (getline(file, line)) {
            linestream.str(line);
            linestream >> currentNodeAsString;
            
            // check for comment lines and empty lines
            if ((std::string(&currentNodeAsString[0]) != std::string("#")) &&
                (std::string(&currentNodeAsString[0]) != std::string(""))) {
                if (i >= numNodes) {
                    throw std::runtime_error(
                                             "CBLoadUnloadedState::LoadNodes: Error, the expected number of nodes is smaller than the number of lines in " + filename +
                                             " !");
                }
                
                currentNode = std::atoi(currentNodeAsString.c_str());
                if (i+1 != currentNode) {
                    std::cout << line << std::endl;
                    throw std::runtime_error(
                                             "CBLoadUnloadedState::LoadNodes: Node indices are not consistent! Do they increase one by one, beginning from 1?");
                }
                
                // read x, y, z
                if (!(linestream >> x >> y >> z)) {
                    std::cout << line << std::endl;
                    throw std::runtime_error(
                                             "CBLoadUnloadedState::LoadNodes: Something went wrong with reading a line from the nodes file " + filename +
                                             "! Did you check if all lines have a sufficient number of coordinates?");
                }
                
                VecSetValue(coordsSeq_, 3*i,   1e-3*x, INSERT_VALUES);
                VecSetValue(coordsSeq_, 3*i+1, 1e-3*y, INSERT_VALUES);
                VecSetValue(coordsSeq_, 3*i+2, 1e-3*z, INSERT_VALUES);
                
                i++;
            }
        }
        if (i < numNodes) {
            throw std::runtime_error(
                                     "CBLoadUnloadedState::LoadNodes: Error, the expected number of nodes is larger than the number of lines in " + filename +
                                     " !");
        }
    }
    file.close();
} // CBLoadUnloadedState::LoadNodes
