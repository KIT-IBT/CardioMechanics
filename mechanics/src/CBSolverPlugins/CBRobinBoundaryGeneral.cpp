/* -------------------------------------------------------
 
 CBRobinBoundaryGeneral.cpp
 
 Ver. 1.1.0
 
 Created:       Tobias Gerach (06.09.2021)
 Last modified: Tobias Gerach (12.11.2021)
 
 Institute of Biomedical Engineering
 Karlsruhe Institute of Technology (KIT)
 
 http://www.ibt.kit.edu
 
 Copyright 2000-2009 - All rights reserved.
 
 ------------------------------------------------------ */

#include "CBSolver.h"
#include "CBRobinBoundaryGeneral.h"

CBRobinBoundaryGeneral::CBRobinBoundaryGeneral() : CBSolverPlugin() {}

CBRobinBoundaryGeneral::~CBRobinBoundaryGeneral() {
    for (auto it : contactSurfaceElements_) {
        delete it;
    }
}

void CBRobinBoundaryGeneral::Init() {
#warning Implementation only works for T3 surface elements at the moment
    
    /// read XML parameter input
    startTime_ =
    parameters_->Get<TFloat>("Plugins.RobinBoundaryGeneral.StartTime", std::numeric_limits<double>::lowest());
    export_ = parameters_->Get<bool>("Plugins.RobinBoundaryGeneral.Export", false);
    
    if (startTime_ != std::numeric_limits<double>::lowest())
        status_ = CBStatus::DACCORD;
    
    /// Depending on the mesh, we sometimes have to flip the surface normals
    bool flipSurfaceNormals = parameters_->Get<bool>("Plugins.RobinBoundaryGeneral.FlipSurfaceNormals", false);
    if (flipSurfaceNormals) {
        normalVectorSign_ = -1;
    }
    
    /// ID of the surface you want to use
    surfaceIndex_ = parameters_->Get<TInt>("Plugins.RobinBoundaryGeneral.SurfaceIndex");
    
    /// Parameters to adjust the force magnitude
    /// alpha_ : stiffness in Pa/m
    /// beta_ : dashpot viscosity in (Pa*s) / m
    alpha_             = parameters_->Get<TFloat>("Plugins.RobinBoundaryGeneral.Alpha", 1e5);
    beta_              = parameters_->Get<TFloat>("Plugins.RobinBoundaryGeneral.Beta", 5e3);
    
    /// solver timestep
    dt_ = Base::adapter_->GetSolver()->GetTiming().GetTimeStep();
    
    InitContactSurfaces();
}  // CBRobinBoundaryGeneral::Init

void CBRobinBoundaryGeneral::InitContactSurfaces() {
    DCCtrl::print << "\t\tcreate contact surfaces..." << std::endl;
    
    for (auto &eleIt : Base::GetAdapter()->GetElementVector()) {
        /// iterate over all CBElement and check if they are of type CBElementContact
        CBElementContactRobin *element = dynamic_cast<CBElementContactRobin *>(eleIt);
        if (element != 0) {
            /// material index
            TInt materialIndex = element->GetSurfaceIndex();
            if (materialIndex <= 0) {
                throw std::runtime_error("CBRobinBoundaryGeneral::InitContactSurfaces(): Surface indices must be > 0");
            }
            if (materialIndex == surfaceIndex_) {
                contactSurfaceElements_.push_back(element);
                ContactForces_.push_back(Vector3<TFloat>(0, 0, 0));
                initialPos_.push_back(element->GetCentroid());
                displacement_.push_back(Vector3<TFloat>(0, 0, 0));
                prevDisplacement_.push_back(Vector3<TFloat>(0, 0, 0));
                velocity_.push_back(Vector3<TFloat>(0, 0, 0));
            }
        }
    }
    DCCtrl::print << "\t\tContact surfaces initialized!" << std::endl;
}  // CBRobinBoundaryGeneral::InitContactSurfaces

void CBRobinBoundaryGeneral::StepBack() {
    for (int i = 0; i < contactSurfaceElements_.size(); i++) {
        displacement_.at(i) = prevDisplacement_.at(i);
    }
    
    stepBack_ = true;
}  // CBRobinBoundaryGeneral::StepBack()

void CBRobinBoundaryGeneral::Apply(TFloat time) {
    /// check if Plugin is actually started
    if (!hasStarted_) {
        if (time > startTime_) {
            hasStarted_ = true;
        } else {
            return;
        }
    }
    
    for (int i = 0; i < contactSurfaceElements_.size(); i++) {
        prevDisplacement_.at(i) = displacement_.at(i);
    }
    
    /// update variables
    prevDt_      = dt_;
    dt_          = Base::adapter_->GetSolver()->GetTiming().GetTimeStep();
    isFirstStep_ = false;
    stepBack_    = false;
}  // CBRobinBoundaryGeneral::Apply

void CBRobinBoundaryGeneral::ApplyToNodalForces() {
    if (!IsActive() || !hasStarted_) {
        return;
    }
    
    for (int i = 0; i < contactSurfaceElements_.size(); i++) {
        auto element = contactSurfaceElements_.at(i);
        TInt nodesCoordsIndices[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
        TFloat nodalForces[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
        bool   bc[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
        
        ContactForces_.at(i) = Vector3<TFloat>(0, 0, 0);
        displacement_.at(i)  = Vector3<TFloat>(0, 0, 0);
        velocity_.at(i)      = Vector3<TFloat>(0, 0, 0);
        
        /// add nodal forces to global vector
        for (int k = 0; k < 3; k++) {
            nodesCoordsIndices[3*k]   = 3.0 * element->GetNodeIndex(k);
            nodesCoordsIndices[3*k+1] = 3.0 * element->GetNodeIndex(k)+1;
            nodesCoordsIndices[3*k+2] = 3.0 * element->GetNodeIndex(k)+2;
        }
        adapter_->ApplyLocalToGlobalMapping(nodesCoordsIndices, 9);
        
        displacement_.at(i) = initialPos_.at(i) - element->GetCentroid();
        velocity_.at(i)     = (displacement_.at(i) - prevDisplacement_.at(i)) / dt_;
        
        CalcForceContributionOfElement(displacement_.at(i), velocity_.at(i), element, nodalForces);
        
        /// respect dirichlet boundary conditions
        Base::GetAdapter()->GetNodesComponentsBoundaryConditionsGlobal(9, nodesCoordsIndices, bc);
        for (int k = 0; k < 3; k++) {
            if (bc[3*k] != 0) {
                nodalForces[3*k] = 0;
            }
            if (bc[3*k+1] != 0) {
                nodalForces[3*k+1] = 0;
            }
            if (bc[3*k+2] != 0) {
                nodalForces[3*k+2] = 0;
            }
            ContactForces_.at(i)    += Vector3<TFloat>(nodalForces[3*k], nodalForces[3*k+1], nodalForces[3*k+2]);
        }
        Base::GetAdapter()->AddNodalForcesComponentsGlobal(9, nodesCoordsIndices, nodalForces);
    }
}  // CBRobinBoundaryGeneral::ApplyToNodalForces

void CBRobinBoundaryGeneral::ApplyToNodalForcesJacobian() {
    //  if (!IsActive() || !hasStarted_) {
    //    return;
    //  }
    //
    //  for (int i = 0; i < contactSurfaceElements_.size(); i++) {
    //    auto element = contactSurfaceElements_.at(i);
    //    TInt nodesCoordsIndices[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
    //    TFloat nodalForces[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
    //    bool   bc[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
    //    TFloat u[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
    //    TFloat v[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
    //
    //    ContactForces_.at(i) = Vector3<TFloat>(0, 0, 0);
    //    ContactDistances_.at(i) = Vector3<TFloat>(0, 0, 0);
    //
    //    /// add nodal forces to global vector
    //    for (int k = 0; k < 3; k++) {
    //      nodesCoordsIndices[3*k]   = 3.0 * element->GetNodeIndex(k);
    //      nodesCoordsIndices[3*k+1] = 3.0 * element->GetNodeIndex(k)+1;
    //      nodesCoordsIndices[3*k+2] = 3.0 * element->GetNodeIndex(k)+2;
    //    }
    //    adapter_->ApplyLocalToGlobalMapping(nodesCoordsIndices, 9);
    //
    //    VecDuplicate(GetSolver()->GetAbsDisplacementVector(), &displacement);
    //    VecCopy(GetSolver()->GetAbsDisplacementVector(), displacement);
    //    VecDuplicate(GetSolver()->GetVelocityVector(), &velocity);
    //    VecCopy(GetSolver()->GetVelocityVector(), velocity);
    //
    //    VecAssemblyBegin(displacement);
    //    VecAssemblyEnd(displacement);
    //    VecAssemblyBegin(velocity);
    //    VecAssemblyEnd(velocity);
    //
    //    VecGetValues(displacement, 9, nodesCoordsIndices, u);
    //    VecGetValues(velocity, 9, nodesCoordsIndices, v);
    //
    //    CalcForceContributionOfElement(u, v, element, nodalForces);
    //
    //    Base::GetAdapter()->GetNodesComponentsBoundaryConditionsGlobal(9, nodesCoordsIndices, bc);
    //
    //    TFloat epsilon = Base::GetAdapter()->GetFiniteDifferencesEpsilon();
    //    TFloat nodalForcesJacobian[9*9];
    //    TFloat nodalForces2[9];
    //    TFloat u2[9] = {0,0,0,0,0,0,0,0,0};
    //
    //    for (int j = 0; j < 9; j++) {
    //      for (int h = 0; h < 9; h++) {
    //        u2[h] = u[h];
    //      }
    //      u2[j] += epsilon;
    //
    //      CalcForceContributionOfElement(u2, v, element, nodalForces2);
    //
    //      for (int k = 0; k < 9; k++) {
    //        if ( (bc[k] == 0) && (bc[j] == 0) ) {
    //          nodalForcesJacobian[9*k+j] = (nodalForces2[k] - nodalForces[j]) / epsilon;
    //        } else {
    //          nodalForcesJacobian[9*k+j] = 0;
    //        }
    //      }
    //    }
    //
    //    Base::GetAdapter()->AddNodalForcesJacobianEntriesGlobal(9, nodesCoordsIndices, 9, nodesCoordsIndices, nodalForcesJacobian);
    //  }
} // CBRobinBoundaryGeneral::ApplyToNodalForcesJacobian

void CBRobinBoundaryGeneral::Export(TFloat time) {
    if (export_) {
        Vec contactForce;
        Vec contactPressure;
        
        Petsc::CreateVector(3*GetAdapter()->GetSolver()->GetNumberOfLocalElements(), PETSC_DETERMINE, &contactForce);
        Petsc::CreateVector(GetAdapter()->GetSolver()->GetNumberOfLocalElements(), PETSC_DETERMINE, &contactPressure);
        
        PetscInt from1, to1;
        PetscInt from2, to2;
        VecGetOwnershipRange(contactForce, &from1, &to1);
        VecGetOwnershipRange(contactPressure, &from2, &to2);
        
        VecZeroEntries(contactPressure);
        VecZeroEntries(contactForce);
        
        for (int i = 0; i < contactSurfaceElements_.size(); i++) {
            auto e = contactSurfaceElements_.at(i);
            Vector3<TFloat> cf = ContactForces_.at(i);
            Vector3<TFloat> u = displacement_.at(i);
            Vector3<TFloat> v = velocity_.at(i);
            
            PetscInt indices[3] =
            {from1 + 3*e->GetLocalIndex(), from1 +  3*e->GetLocalIndex()+1, from1 + 3*e->GetLocalIndex()+2};
            PetscScalar f[3] = {-cf(0), -cf(1), -cf(2)}; // negative, to make the vectors point in the right direction
            VecSetValues(contactForce, 3, indices, f, INSERT_VALUES);
            
            PetscScalar p = cf.Norm() / e->GetArea();
            VecSetValue(contactPressure, from2 + e->GetLocalIndex(), p, INSERT_VALUES);
        }
        
        GetAdapter()->GetSolver()->ExportElementsVectorData("ContactForce", contactForce);
        GetAdapter()->GetSolver()->ExportElementsScalarData("ContactPressure", contactPressure);
        
        // free memory
        VecDestroy(&contactPressure);
        VecDestroy(&contactForce);
    }
} // CBRobinBoundaryGeneral::Export

void CBRobinBoundaryGeneral::WriteToFile(TFloat time) {}

void CBRobinBoundaryGeneral::Prepare() {}

void CBRobinBoundaryGeneral::CalcForceContributionOfElement(Vector3<TFloat> u, Vector3<TFloat> v,
                                                            CBElementContactRobin *triangle, TFloat *nodalForces) {
    /// we use a one point quadrature rule for the integration on the linear triangle element e
    /// therefore, the force f at node I is given with
    /// f_i = - A_e * sum_i^n[ W * N_i(l1, l2, l3) * p * normalVec ]
    /// n = 1
    /// W = 1
    /// l1 = l2 = l3 = 1/3
    /// p = alpha * u  + beta * v
    TFloat area = triangle->GetTriangle().GetArea();
    TFloat W    = 1;
    TFloat scaling = triangle->GetSurfaceTractionScaling();
    
    /// Shape functions T3 element
    std::function<double(double, double, double)> Ni[3];
    
    Ni[0] = [](double l1, double l2, double l3) {return l1; };
    Ni[1] = [](double l1, double l2, double l3) {return l2; };
    Ni[2] = [](double l1, double l2, double l3) {return l3; };
    
    /// iterate over nodes
    for (unsigned int i = 0; i < 3; i++) {
        Vector3<TFloat> forceMagnitude = scaling * alpha_ * u + scaling * beta_ * v;
        Vector3<TFloat> f = -area * forceMagnitude * W * Ni[i](1.0/3.0, 1.0/3.0, 1.0/3.0);
        
        nodalForces[3*i]   = f.X();
        nodalForces[3*i+1] = f.Y();
        nodalForces[3*i+2] = f.Z();
    }
}  // CBRobinBoundaryGeneral::CalcForceContributionOfElement
