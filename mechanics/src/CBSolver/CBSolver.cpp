/*
 *  CBSolver.cpp
 *  CardioMechanics
 *
 *  Created by Thomas Fritz on 13.09.11.
 *  Copyright 2011 Institute for Biomedical Engineering, Karlsruhe Institute of Technology (KIT). All rights reserved.
 *
 */
#include "CBSolver.h"

#include "CBFileManager.h"
#include "CBElementSurfaceT6.h"
#include "CBElementSurfaceT3.h"
#include "CBRuntimeEstimator.h"
#include "CBModelLoader.h"
#include "CBModelLoaderTetgen.h"
#include <iostream>
#include <fstream>

#undef CHKERRQ
#define CHKERRQ(n) do {\
if (PetscUnlikely(n)) {                  \
PetscError(PETSC_COMM_SELF, __LINE__, PETSC_FUNCTION_NAME, __FILE__, (PetscErrorCode)n, PETSC_ERROR_REPEAT, " "); \
throw std::runtime_error("PETSc error"); \
}\
} while (0)

void CBSolver::InitParameters() {
    if (parameters_ == 0) {
        throw std::runtime_error(
                                 "CBSolver::InitSolverOptions(): void CBSolver::Init(ParameterMap* parameters, CBModel* model) has to be run before");
    }
    precision_   = parameters_->Get<TFloat>("Solver.Precision", 1e-8);
    epsilon_     = parameters_->Get<TFloat>("Solver.Epsilon", 1e-12);
    adapter_->SetFiniteDifferencesEpsilon(epsilon_);
    ignoreMaxIt_ = parameters_->Get<bool>("Solver.IgnoreMaxIt", false);
    maxSnesIts_  = parameters_->Get<TInt>("Solver.MaxIterations", 50);
    maxFunEval_  = parameters_->Get<TInt>("Solver.MaxFunEval", 1000);
    numNonZeros_ = parameters_->Get<TInt>("Solver.NonZeros", 2000);
    timing_.Init(*parameters_);
    exportSNESMatrix_ = parameters_->Get<bool>("Solver.ExportSNESMatrix", false);
    
    // snesStepsFilename_ = parameters_->Get<std::string>("Solver.SNES.Export","");  // unused
    
    if (timing_.GetMaxSimTimeStep() > model_->GetExporter()->GetExportTimeStep()) {
        DCCtrl::print
        <<
        "CBSolver::InitSolverOptions(): Max simulation time step is bigger than the output time step, Does this make sense ? I don't think so ... "
        << "Will set export time step to match max simulation time step " << timing_.GetMaxSimTimeStep() << "."
        << std::endl;
        model_->GetExporter()->SetExportTimeStep(timing_.GetMaxSimTimeStep());
    }
    stepUpIterations_ = parameters_->Get<int>("Solver.StepUpIterations", 1000); // don't increase time step when snes took too many iterations (no longer used)
    enableTiming_ = parameters_->Get<bool>("Solver.Timing", true);  // show output of the runtime estimation
} // CBSolver::InitParameters

void CBSolver::InitFormulation() {
    auto formulationType = parameters_->Get<std::string>("Solver.Formulation", "TotalLagrangian");
    
    if (formulationType != "TotalLagrangian")
        throw std::runtime_error("CBSolver::InitFormulation(): Type of formulation: " + formulationType + " is unkown.");
    
    formulation_ = new CBFormulationTotalLagrangian(this);
    isInitFormulationDone_ = true;
}

void CBSolver::ExtractSolidElements() {
    if (!isInitElementsDone_)
        throw std::runtime_error("CBSolver::ExtractSolidElements(): InitElements() has to be run first");
    
    for (auto &it : elements_) {
        auto *element = dynamic_cast<CBElementSolid *>(it);
        if (element) {
            solidElements_.push_back(element);
        }
    }
    isExtractSolidElementsDone_ = true;
}

void CBSolver::ExtractSurfaceElements() {
    if (!isInitElementsDone_)
        throw std::runtime_error("CBSolver::ExtractSurfaceElements(): InitElements() has to be run first");
    
    for (auto &it : elements_) {
        auto *element = dynamic_cast<CBElementSurface *>(it);
        if (element) {
            surfaceElements_.push_back(element);
        }
    }
}

void CBSolver::InitMaterials() {
    if (!isExtractSolidElementsDone_)
        throw std::runtime_error("void CBSolver::InitMaterials(): ExtractSolidElements() has to be run first");
    std::set<TInt> materialIndexes;
    
    for (auto &it : solidElements_) {
        auto matIndex = it->GetMaterialIndex();
        if (materialIndexes.insert(matIndex).second) {
            auto *mat = materialFactory_.New(parameters_, matIndex);
            it->SetMaterial(mat);
            materials_.insert(std::pair<TInt, CBMaterial *>(matIndex, mat));
        } else {
            CBMaterial *mat = materials_.find(matIndex)->second;
            it->SetMaterial(mat);
        }
    }
}

void CBSolver::InitActiveStress() {
    // unscaled active Stress
    std::cout << activeStressData_ << std::endl;
    
    if (activeStressData_) {
        PetscInt numActiveStressTensorComponents = 0;
        activeStressTensorIndices_    = new int[1];
        activeStressTensorIndices_[0] = 0;
        
        if (DCCtrl::GetNumberOfProcesses() == 1)
            VecCreateSeq(PETSC_COMM_SELF, model_->GetElements().size(), &activeStress_);
        else
            VecCreateMPI(Petsc::Comm(), numLocalElements_, PETSC_DECIDE, &activeStress_);
        
        VecGetOwnershipRange(activeStress_, &activeStressLowerIndex_, &activeStressUpperIndex_);
        
        activeStressTensorComponents_        = new PetscScalar[numLocalElements_];
        activeStressTensorComponentsIndices_ = new PetscInt[numLocalElements_];
        numActiveStressTensorComponents      = 1;
        adapter_->LinkActiveStress(activeStress_, numActiveStressTensorComponents, activeStressTensorIndices_);
    } else {
        activeStress_ = 0;
        adapter_->LinkActiveStress(0, 0, 0);
    }
    
    // exportActiveStress (scaled)
    
    if (activeStressData_) {
        PetscInt numActiveStressTensorComponents = 0;
        activeStressTensorIndices_    = new int[1];
        activeStressTensorIndices_[0] = 0;
        
        if (DCCtrl::GetNumberOfProcesses() == 1)
            VecCreateSeq(PETSC_COMM_SELF, model_->GetElements().size(), &exportActiveStress_);
        else
            VecCreateMPI(Petsc::Comm(), numLocalElements_, PETSC_DECIDE, &exportActiveStress_);
        
        VecGetOwnershipRange(exportActiveStress_, &activeStressLowerIndex_, &activeStressUpperIndex_);
        
        activeStressTensorComponents_        = new PetscScalar[numLocalElements_];
        activeStressTensorComponentsIndices_ = new PetscInt[numLocalElements_];
        numActiveStressTensorComponents      = 1;
        adapter_->LinkActiveStress(exportActiveStress_, numActiveStressTensorComponents, activeStressTensorIndices_);
    } else {
        exportActiveStress_ = 0;
        adapter_->LinkActiveStress(0, 0, 0);
    }
} // CBSolver::InitActiveStress

void CBSolver::InitActiveStressData() {
    if (!isInitNodesDone_ || !isInitElementsDone_)
        throw std::runtime_error("CBSolver::InitActiveStress(): InitNodes() and InitElements() has to be finished first");
    
    activeStressData_ = new CBDataPerMaterial(parameters_);
    
    InitActiveStress();
    
    isInitActiveStressDataDone_ = true;
}

void CBSolver::InitTensionModels() {
    // init array for ExportActiveStressData_ (which is activeStress_,
    // but could actually be initialized using InitActiveStressVectors() or so)
    if (!activeStress_) {
        // PetscInt numActiveStressTensorComponents = 0;
        activeStressTensorIndices_    = new int[1];
        activeStressTensorIndices_[0] = 0;
        
        if (DCCtrl::GetNumberOfProcesses() == 1)
            VecCreateSeq(PETSC_COMM_SELF, model_->GetElements().size(), &activeStress_);
        else
            VecCreateMPI(Petsc::Comm(), numLocalElements_, PETSC_DECIDE, &activeStress_);
        
        VecGetOwnershipRange(activeStress_, &activeStressLowerIndex_, &activeStressUpperIndex_);
        
        activeStressTensorComponents_        = new PetscScalar[numLocalElements_];
        activeStressTensorComponentsIndices_ = new PetscInt[numLocalElements_];
    } else {
        activeStress_ = 0;
    }
    
    if (!isExtractSolidElementsDone_)
        throw std::runtime_error("void CBSolver::InitActiveStress(): ExtractSolidElements() has to be run first");
    
    // Initialize FileManager
    fileManager_.Init(parameters_, &timing_);
    tensionFactory_.Init(parameters_, &timing_, &fileManager_);
    
    // The lat offset which is added to the activation time of every element. (Given in the xml-file)
    float latOffset = fileManager_.GetLatOffset();
    
    // generate necessary CBData* objects using a factory and assign them to the corresponding element
    for (auto solidElement : solidElements_) {
        solidElement->SetTensionModel(tensionFactory_.New(solidElement) );
        
        // set individual activation times for each element, if a LAT file was given.
        if (fileManager_.GetProcessLAT()) {
            solidElement->GetTensionModel()->AddToActivationTime(fileManager_.GetActivationTimesVector().at(solidElement->
                                                                                                            GetIndex()).second);
        } else {
            solidElement->GetTensionModel()->AddToActivationTime(0.0 + latOffset);
        }
    }
} // CBSolver::InitTensionModels

void CBSolver::Init(ParameterMap *parameters, CBModel *model) {
    parameters_ = parameters;
    model_  = model;
    
    adapter_ = new CBElementAdapter();
    adapter_->SetSolver(this);
    
    InitParameters();
    model_->SetMaxNnz(numNonZeros_);
    LoadMesh();
    
    MPI_Barrier(Petsc::Comm());
    
    InitFormulation();
    
    // DCCtrl::print << "\xd\t\tInitializing: Active stress data ...                                                                                ";
    // InitActiveStressData();
    // InitActiveStress();
    DCCtrl::print <<
    "\xd\t\tInitializing: Tension models ...                                                                                    ";
    InitTensionModels();
    DCCtrl::print <<
    "\xd\t\tInitializing: Vectors ...                                                                                           ";
    InitVectors();
    DCCtrl::print <<
    "\xd\t\tInitializing: Matrices ...                                                                                          ";// ( Hint: if this takes too long preallocate more memory for the matrices [solver.nonZeros] ...";
    InitMatrices();
    DCCtrl::print <<
    "\xd\t\tInitializing: Data exporter ...                                                                                     ";
    InitExporter();
    DCCtrl::print <<
    "\xd\t\tInitializing: PETSc ...                                                                                             ";
    InitPETScSolver();
    DCCtrl::print <<
    "\xd\t\tInitializing: Plugins ...                                                                                           "
    << std::endl;
    InitPlugins();
    DCCtrl::print <<
    "\xd\t\tInitializing solver finished !                                                                                      ";
    DCCtrl::print << "\n\n\t\t\tActive Plugins:\n";
    if (plugins_.size() == 0)
        DCCtrl::print << "\t\t\t - " << "None" << "\n";
    for (auto it : plugins_)
        DCCtrl::print << "\t\t\t - " << it->GetName() << "\n";
    
    if (LoadedModel_ != 0) {
        LoadVelocityAndAcceleration();
        LoadActiveStress();
    }
    
    MPI_Barrier(Petsc::Comm());
} // CBSolver::Init

void CBSolver::LoadVelocityAndAcceleration() {
    if (LoadedModel_->IsVelocity())
        SetVelocity(LoadedModel_->GetVelocity());
    if (LoadedModel_->IsAcceleration())
        SetAcceleration(LoadedModel_->GetAcceleration());
}

void CBSolver::LoadActiveStress() {
    if (activeStressData_ && LoadedModel_->IsActiveStress()) {
        std::vector<TFloat> AS = LoadedModel_->GetActiveStress();
        
        for (auto &it : solidElements_) {
            int iL = it->GetLocalIndex();
            int iG = it->GetIndex();
            TFloat as = AS[iG];
            VecSetValue(activeStress_, iL, as, INSERT_VALUES);
        }
        
        VecAssemblyBegin(activeStress_);
        VecAssemblyEnd(activeStress_);
    }
}

CBStatus CBSolver::CalcNodalForces() {
    CBStatus rc = CBStatus::FAILED;
    PetscErrorCode prc1 = VecZeroEntries(nodalForces_);
    
    adapter_->LinkNodalForces(nodalForces_);
    UpdateGhostNodesAndLinkToAdapter();
    rc = formulation_->CalcNodalForces();
    PetscErrorCode prc2 = VecAssemblyBegin(nodalForces_);
    PetscErrorCode prc3 = VecAssemblyEnd(nodalForces_);
    if ((prc1 == 0) && (prc2 == 0) && (prc3 == 0) && (rc == CBStatus::SUCCESS) )
        rc = CBStatus::SUCCESS;
    return rc;
}

CBStatus CBSolver::CalcNodalForcesJacobian() {
    CBStatus rc = CBStatus::FAILED;
    PetscErrorCode prc0 = MatZeroEntries(nodalForcesJacobian_);
    
    CreateNodesJacobianAndLinkToAdapter();
    UpdateGhostNodesAndLinkToAdapter();
    
    formulation_->CalcNodalForcesJacobian();
    
    for (auto &it : plugins_)
        it->ApplyToNodalForcesJacobian();
    PetscErrorCode prc1 = MatAssemblyBegin(nodalForcesJacobian_, MAT_FINAL_ASSEMBLY);
    PetscErrorCode prc2 = MatAssemblyEnd(nodalForcesJacobian_, MAT_FINAL_ASSEMBLY);
    PetscErrorCode prc3 = MatAXPY(nodalForcesJacobian_, 1, boundaryConditionsNodalForcesJacobianDiagonalComponents_,
                                  DIFFERENT_NONZERO_PATTERN);
    if ((prc0 == 0) && (prc1 == 0) && (prc2 == 0) && (prc3 == 0))
        rc = CBStatus::SUCCESS;
    return rc;
}

const std::vector<CBElement *> & CBSolver::GetElementVector() {
    return elements_;
}

const std::vector<CBElementSolid *> & CBSolver::GetSolidElementVector() {
    return solidElements_;
}

const std::vector<CBSolverPlugin *> & CBSolver::GetPlugInVector() {
    return plugins_;
}

void CBSolver::GetDomainDecomposition(ParameterMap *parameters, CBModel *model, std::vector<std::pair<int,
                                      int>> &mapping) {
    DCCtrl::print << "\tDomain Decomposition ... \n\n";
    DCCtrl::print << "\t\tNumber of processes: " << DCCtrl::GetNumberOfProcesses() << "." << "\n\n";
    DCCtrl::print << "\t\tInitializing: ";
    
    Init(parameters, model);
    
    DCCtrl::print <<
    "\xd\t\tInitializing: done !!!                                                                                             \n";
    
    CreateNodesJacobianAndLinkToAdapter();
    DCCtrl::print << "\t\tGenerating connectivity matrix ... ";
    
    Vec displacement;
    VecDuplicate(nodes_, &displacement);
    
    CalcNodalForcesJacobian(displacement, nodalForcesJacobian_);
    
    Mat adj;
    IS is, isg;
    IS isGlobal, isgGlobal;
    
    MatConvert(nodalForcesJacobian_, MATMPIADJ, MAT_INITIAL_MATRIX, &adj);
    
    MatPartitioning part;
    
    MatPartitioningCreate(Petsc::Comm(), &part);
    MatPartitioningSetAdjacency(part, adj);
    MatPartitioningSetFromOptions(part);
    MatPartitioningApply(part, &is);
    
    ISPartitioningToNumbering(is, &isg);
    
    ISAllGather(is, &isGlobal);
    ISAllGather(isg, &isgGlobal);
    
    PetscInt size;
    ISGetSize(is, &size);
    
    const PetscInt *j;
    const PetscInt *k;
    
    ISGetIndices(isgGlobal, &j);
    ISGetIndices(isGlobal, &k);
    
    std::vector<PetscInt> m;
    m.resize(size, 0);
    
    for (int i = 0; i < size; i++)
        m.at(j[i]) = i;
    
    mapping.resize(size/3, std::pair<int, int>(0, 0));
    
    int n = 0;
    
    for (int i = 0; i < size; i++) {
        if (m.at(i)%3 == 0) {
            std::pair<int, int> p(n, k[m.at(i)]);
            mapping.at(m.at(i)/3) = p;
            n++;
        }
    }
    
    MatPartitioningDestroy(&part);
    MatDestroy(&adj);
    ISDestroy(&is);
    ISDestroy(&isg);
    DCCtrl::print << " done.";
} // CBSolver::GetDomainDecomposition

void CBSolver::SetNodesComponentsBoundaryConditionsGlobal(bool *globalNodesComponentsBoundaryConditions) {
    globalNodesComponentsBoundaryConditions_ = globalNodesComponentsBoundaryConditions;
    for (int i = 0; i < numLocalNodes_; i++) {
        bool bc[3];
        PetscInt indices[3] = {3*i, 3*i+1, 3*i+2};
        adapter_->GetNodesComponentsBoundaryConditions(3, indices, bc);
        PetscScalar val[9] = {
            (bc[0] == true) ? 1.0 : 0.0,
            0,
            0,
            0,
            (bc[1] == true) ? 1.0 : 0.0,
            0,
            0,
            0,
            (bc[2] == true) ? 1.0 : 0.0};
        MatSetValuesLocal(boundaryConditionsNodalForcesJacobianDiagonalComponents_, 3, indices, 3, indices, val,
                          INSERT_VALUES);
    }
    
    MatAssemblyBegin(boundaryConditionsNodalForcesJacobianDiagonalComponents_, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(boundaryConditionsNodalForcesJacobianDiagonalComponents_, MAT_FINAL_ASSEMBLY);
}

void CBSolver::InitNodalForces() {
    if (!isInitNodesDone_)
        throw std::runtime_error("CBSolver::InitNodes() has to be run before CBSolver::InitNodalForces()");
    
    if (DCCtrl::IsParallel())
        VecCreateMPI(Petsc::Comm(), 3 * numLocalNodes_, PETSC_DECIDE, &nodalForces_);
    else
        VecCreateSeq(PETSC_COMM_SELF, 3 * numNodes_, &nodalForces_);
    
    VecSetLocalToGlobalMapping(nodalForces_, nodesIndicesMapping_);
    adapter_->LinkNodalForces(nodalForces_);
}

void CBSolver::InitNodalForcesJacobian() {
    if (!isInitNodesDone_)
        throw std::runtime_error("CBSolver::InitNodes() has to be run before CBSolver::InitNodalForcesJacobian()");
    
    CreateNodesJacobianAndLinkToAdapter();
    UpdateGhostNodesAndLinkToAdapter();
    
    adapter_->LinkNodalForcesJacobian(nodalForcesJacobian_);
    formulation_->CalcNodalForcesJacobian();
    
    isInitNodalForcesJacobianDone_ = true;
}

void CBSolver::UpdateGhostNodesAndLinkToAdapter() {
    if (DCCtrl::IsParallel()) {
        if (nodesSeq_) {
            VecDestroy(&nodesSeq_);
            nodesSeq_ = 0;
        }
        
        if (refNodesSeq_) {
            VecDestroy(&refNodesSeq_);
            refNodesSeq_ = 0;
        }
        
        VecGhostUpdateBegin(nodes_, INSERT_VALUES, SCATTER_FORWARD);
        VecGhostUpdateEnd(nodes_, INSERT_VALUES, SCATTER_FORWARD);
        VecGhostGetLocalForm(nodes_, &nodesSeq_);
        adapter_->LinkNodes(nodesSeq_);
        
        VecGhostUpdateBegin(refNodes_, INSERT_VALUES, SCATTER_FORWARD);
        VecGhostUpdateEnd(refNodes_, INSERT_VALUES, SCATTER_FORWARD);
        VecGhostGetLocalForm(refNodes_, &refNodesSeq_);
        adapter_->LinkRefNodes(refNodesSeq_);
    } else {
        adapter_->LinkNodes(nodes_);
        adapter_->LinkRefNodes(refNodes_);
    }
} // CBSolver::UpdateGhostNodesAndLinkToAdapter

void CBSolver::InitPlugins() {
    UpdateGhostNodesAndLinkToAdapter();
    adapter_->LinkNodalForcesJacobian(nodalForcesJacobian_);
    pluginFactory_ = new CBSolverPluginFactory();
    pluginFactory_->LoadAllPlugins(plugins_, parameters_);
    
    for (auto &it : plugins_) {
        it->SetAdapter(adapter_);
        it->Init();
        it->ApplyToNodalForcesJacobian();
    }
}

Matrix3fVector CBSolver::GetDeformationTensors() {
    Matrix3fVector tensors;
    
    for (auto &it : solidElements_) {
        Matrix3f f;
        it->GetDeformationTensor(f);
        tensors.push_back(f);
    }
    
    return tensors;
}

void CBSolver::UpdateActiveStress(PetscScalar time) {
    //    PetscInt cnt = 0;
    //  if(activeStressData_) {
    //    for(auto& it : solidElements_) {
    
    // todo: this should be able to depend dynamically on everything instead of being set just once, similar to the constitutive model
    
    //        activeStressData_->SetTime(time);
    // it->GetMaterial()->GetConstitutiveModel()->SetTemplateForce(activeStressData_->Get(time, it->GetIndex(), it->GetMaterialIndex()));
    // active stress currently gets scaled by the constitutive model, but should rather be in ActiveStressModel
    
    // todo: instead of setting this here to zero, it should be remove from CalcNodalForcesHelperFunction in CBElementT4, CBElementT10, (partially done) and so on...
    // activeStressTensorComponents_[cnt]        =  activeStressData_->Get(time, it->GetIndex(), it->GetMaterialIndex() );
    // activeStressTensorComponentsIndices_[cnt] = it->GetLocalIndex() + activeStressLowerIndex_;
    // cnt++;
    //        }
    
    //        VecSetValues(activeStress_, cnt, activeStressTensorComponentsIndices_, activeStressTensorComponents_, INSERT_VALUES);
    //        VecAssemblyBegin(activeStress_);
    //        VecAssemblyEnd(activeStress_);
    //    }
}

void CBSolver::Export(TFloat timeStep) {
    UpdateGhostNodesAndLinkToAdapter();
    
    model_->SetCurrentTime(timeStep);
    
    ExportNodes();
    ExportBasisAtQuadraturePoint();
    
    // global data
    
    UpdateExportGlobalData("ModelVolume", GetVolumeOfSolidElements());
    UpdateExportGlobalData("DeformationEnergy", GetDeformationEnergy());
    if (UpdateExportGlobalData("SnesIterations", snesIterations_))
        snesIterations_ = 0;
    if (UpdateExportGlobalData("KspIterations", kspIterations_))
        kspIterations_ = 0;
    
    // (global?) plugin data
    for (auto &it : plugins_)
        it->Export(timeStep);
    
    // per-element data
    
    std::string str = "ActiveStress";
    if (model_->GetExporter()->GetExportOption(str, true) && activeStress_) {
        UpdateExportActiveStress(model_->GetCurrentTime(), activeStress_);
        ExportElementsScalarData(str, activeStress_);
    }
    
    // str = "ActiveStress";
    // if(model_->GetExporter()->GetExportOption(str, true) && exportActiveStress_) {
    //   UpdateExportActiveStress(model_->GetCurrentTime(), exportActiveStress_);
    //   ExportElementsScalarData(str, exportActiveStress_);
    // }
    
    ExportFiber();
    ExportJacobian();
    ExportSheet();
    ExportSheetNormal();
    ExportBasesAtAllQuadraturePoints();
    ExportDeformation();
    ExportLambda();
    ExportCauchy();
    ExportPK2Stress();
    ExportGreenLagrangeStrain();
    ExportLocalActivationTime();
} // CBSolver::Export

/// Takes the active tension (a scalar value, that is basically AS(0,0) ) for
/// the latest time step and puts it into a vector that can later be exported by
/// ExportElementsScalarData() .
void CBSolver::UpdateExportActiveStress(PetscScalar time, Vec &stressVec) {
    // deprecated: does basically the same as updateActiveStress(time), but uses the function
    // scaled with Tmax and depending on lambda which is located falsely in
    // CBConstitutiveModel, but should rather be in a CBStressModel/CBTensionModel
    PetscInt cnt = 0;
    
    // if (activeStressData_) {
    if (activeStressTensorComponents_) {
        for (auto &it : solidElements_) {
            Matrix3<TFloat> deformationTensor;
            it->GetDeformationTensor(deformationTensor);
            
            TFloat as = it->GetTensionModel()->CalcActiveTension(deformationTensor, time);
            activeStressTensorComponents_[cnt]        = as;
            activeStressTensorComponentsIndices_[cnt] = it->GetLocalIndex() + activeStressLowerIndex_;
            cnt++;
        }
        VecSetValues(stressVec, cnt, activeStressTensorComponentsIndices_, activeStressTensorComponents_, INSERT_VALUES);
        VecAssemblyBegin(stressVec);
        VecAssemblyEnd(stressVec);
    }
} // CBSolver::UpdateExportActiveStress

bool CBSolver::UpdateExportGlobalData(std::string identifier, TFloat value) {
    bool result = model_->GetExporter()->GetExportOption(identifier, true);
    
    if (result)
        model_->SetGlobalData(identifier, value);
    return result;
}

/// This exports sigma = J^-1 F S F^T
void CBSolver::ExportCauchy() {
    if (model_->GetExporter()->GetExportOption("Cauchy", false)) {
        Vec cauchyDiag;
        Vec cauchyOff;
        
        PetscInt from, to;
        
        if (DCCtrl::IsParallel())
            VecCreateMPI(Petsc::Comm(), 3 * numLocalElements_, PETSC_DETERMINE, &cauchyDiag);
        else
            VecCreateSeq(PETSC_COMM_SELF, 3*numElements_, &cauchyDiag);
        
        VecGetOwnershipRange(cauchyDiag, &from, &to);
        
        VecDuplicate(cauchyDiag, &cauchyOff);
        
        VecZeroEntries(cauchyDiag);
        VecZeroEntries(cauchyOff);
        
        for (auto &it : solidElements_) {
            Matrix3f cauchy;
            it->GetCauchyStress(cauchy);
            int i = it->GetLocalIndex();
            PetscScalar diag[3] = {cauchy(0, 0), cauchy(1, 1), cauchy(2, 2)};
            PetscScalar off[3] = {cauchy(0, 1), cauchy(0, 2), cauchy(1, 2)};
            PetscInt index[3] = {from+3*i, from +3*i+1, from +3*i+2};
            VecSetValues(cauchyDiag, 3, index, diag, INSERT_VALUES);
            VecSetValues(cauchyOff, 3, index, off, INSERT_VALUES);
        }
        
        VecAssemblyBegin(cauchyDiag);
        VecAssemblyBegin(cauchyOff);
        
        VecAssemblyEnd(cauchyDiag);
        VecAssemblyEnd(cauchyOff);
        
        ExportElementsVectorData("Sigma_ii", cauchyDiag);
        ExportElementsVectorData("Sigma_ij", cauchyOff);
        
        VecDestroy(&cauchyDiag);
        VecDestroy(&cauchyOff);
    }
} // CBSolver::ExportCauchy

/// export J = Det(F) to keep track of volume changes / incompressability constraint
void CBSolver::ExportJacobian() {
    if (model_->GetExporter()->GetExportOption("Jacobian", false)) {
        Vec jacobian;
        PetscInt from, to;
        if (DCCtrl::IsParallel())
            VecCreateMPI(Petsc::Comm(), numLocalElements_, PETSC_DETERMINE, &jacobian);
        else
            VecCreateSeq(PETSC_COMM_SELF, numElements_, &jacobian);
        
        VecGetOwnershipRange(jacobian, &from, &to);
        VecZeroEntries(jacobian);
        
        for (auto &it : solidElements_) {
            Matrix3f F;
            it->GetDeformationTensor(F);
            PetscScalar J = F.Det();
            
            int i = it->GetLocalIndex();
            PetscInt index = from+i;
            VecSetValue(jacobian, index, J, INSERT_VALUES);
        }
        VecAssemblyBegin(jacobian);
        VecAssemblyEnd(jacobian);
        ExportElementsScalarData("J", jacobian);
        VecDestroy(&jacobian);
    }
} // CBSolver::ExportJacobian

/// write fiber vector, respects deformation, which may alter its orientation and length
void CBSolver::ExportFiber() {
    if (model_->GetExporter()->GetExportOption("Fiber", true)) {
        Vec fiberRotated;
        PetscInt from, to;
        
        if (DCCtrl::IsParallel())
            VecCreateMPI(Petsc::Comm(), 3 * numLocalElements_, PETSC_DETERMINE, &fiberRotated);
        else
            VecCreateSeq(PETSC_COMM_SELF, 3*numElements_, &fiberRotated);
        
        VecGetOwnershipRange(fiberRotated, &from, &to);
        
        VecZeroEntries(fiberRotated);
        
        for (auto &it : solidElements_) {
            Matrix3f F;
            it->GetDeformationTensor(F);
            Vector3<TFloat> fr = F * it->GetBasis()->GetCol(0);  // fiber rotated
            
            int i = it->GetLocalIndex();
            PetscScalar fibrot[3] = {fr(0), fr(1), fr(2)};
            PetscInt index[3] = {from+3*i, from +3*i+1, from +3*i+2};
            VecSetValues(fiberRotated, 3, index, fibrot, INSERT_VALUES);
        }
        
        VecAssemblyBegin(fiberRotated);
        VecAssemblyEnd(fiberRotated);
        ExportElementsVectorData("f", fiberRotated);
        VecDestroy(&fiberRotated);
    }
} // CBSolver::ExportFiber

/// write sheet vector
void CBSolver::ExportSheet() {
    if (model_->GetExporter()->GetExportOption("Sheet", true)) {
        Vec sheetRotated;
        PetscInt from, to;
        
        if (DCCtrl::IsParallel())
            VecCreateMPI(Petsc::Comm(), 3 * numLocalElements_, PETSC_DETERMINE, &sheetRotated);
        else
            VecCreateSeq(PETSC_COMM_SELF, 3*numElements_, &sheetRotated);
        
        VecGetOwnershipRange(sheetRotated, &from, &to);
        
        VecZeroEntries(sheetRotated);
        
        for (auto &it : solidElements_) {
            Matrix3f F;
            it->GetDeformationTensor(F);
            Vector3<TFloat> sr = F * it->GetBasis()->GetCol(1);  // sheet rotated
            
            int i = it->GetLocalIndex();
            PetscScalar sheetrot[3] = {sr(0), sr(1), sr(2)};
            PetscInt index[3] = {from+3*i, from +3*i+1, from +3*i+2};
            VecSetValues(sheetRotated, 3, index, sheetrot, INSERT_VALUES);
        }
        
        VecAssemblyBegin(sheetRotated);
        VecAssemblyEnd(sheetRotated);
        ExportElementsVectorData("s", sheetRotated);
        VecDestroy(&sheetRotated);
    }
} // CBSolver::ExportSheet

/// write sheetnormal vector
void CBSolver::ExportSheetNormal() {
    if (model_->GetExporter()->GetExportOption("Normal", true)) {
        Vec sheetnormalRotated;
        PetscInt from, to;
        
        if (DCCtrl::IsParallel())
            VecCreateMPI(Petsc::Comm(), 3 * numLocalElements_, PETSC_DETERMINE, &sheetnormalRotated);
        else
            VecCreateSeq(PETSC_COMM_SELF, 3*numElements_, &sheetnormalRotated);
        
        VecGetOwnershipRange(sheetnormalRotated, &from, &to);
        
        VecZeroEntries(sheetnormalRotated);
        
        for (auto &it : solidElements_) {
            Matrix3f F;
            it->GetDeformationTensor(F);
            Vector3<TFloat> snr = F * it->GetBasis()->GetCol(2);  // sheetnormal rotated
            
            int i = it->GetLocalIndex();
            PetscScalar sheetnormalrot[3] = {snr(0), snr(1), snr(2)};
            PetscInt index[3] = {from+3*i, from +3*i+1, from +3*i+2};
            VecSetValues(sheetnormalRotated, 3, index, sheetnormalrot, INSERT_VALUES);
        }
        
        VecAssemblyBegin(sheetnormalRotated);
        VecAssemblyEnd(sheetnormalRotated);
        ExportElementsVectorData("n", sheetnormalRotated);
        VecDestroy(&sheetnormalRotated);
    }
} // CBSolver::ExportSheetNormal

void CBSolver::ExportBasesAtAllQuadraturePoints() {
    if (model_->GetExporter()->GetExportOption("TetgenBases", false)) {
        TInt nQ = model_->GetSolidElements().at(0)->GetNumberOfQuadraturePoints();
        for (auto &e : model_->GetSolidElements()) {
            if (e->GetNumberOfQuadraturePoints() != nQ)
                throw std::runtime_error("ERROR: Tetgen export is currently only implemented for all T4 or all T10 elements.");
        }
        
        for (int q = 0; q < nQ; q++) {
            Vec fVec, sVec, snVec;
            PetscInt from1, to1;
            PetscInt from2, to2;
            PetscInt from3, to3;
            
            if (DCCtrl::IsParallel())
                VecCreateMPI(Petsc::Comm(), 3 * numLocalElements_, PETSC_DETERMINE, &fVec);
            else
                VecCreateSeq(PETSC_COMM_SELF, 3*numElements_, &fVec);
            VecDuplicate(fVec, &sVec);
            VecDuplicate(fVec, &snVec);
            
            VecGetOwnershipRange(fVec, &from1, &to1);
            VecGetOwnershipRange(sVec, &from2, &to2);
            VecGetOwnershipRange(snVec, &from3, &to3);
            
            VecZeroEntries(fVec);
            VecZeroEntries(sVec);
            VecZeroEntries(snVec);
            
            for (auto &it : solidElements_) {
                Matrix3f F;
                it->GetDeformationTensor(F);
                Vector3<TFloat> f  = F * it->GetBasisAtQuadraturePoint(q)->GetCol(0);
                Vector3<TFloat> s  = F * it->GetBasisAtQuadraturePoint(q)->GetCol(1);
                Vector3<TFloat> sn = F * it->GetBasisAtQuadraturePoint(q)->GetCol(2);
                
                // rescale f to length = 1
                f.Normalize();
                
                // Gram Schmidt orthogonalize s wrt. f
                s.Normalize();
                s = s - (f*s) * f;
                
                // sn needs to be orthogonal to the others
                sn.Normalize();
                sn = CrossProduct(f, s);
                
                int idx = it->GetLocalIndex();
                PetscInt index1[3] = {from1+3*idx, from1+3*idx+1, from1+3*idx+2};
                VecSetValues(fVec, 3, index1, f.GetArray(), INSERT_VALUES);
                PetscInt index2[3] = {from2+3*idx, from2+3*idx+1, from2+3*idx+2};
                VecSetValues(sVec, 3, index2, s.GetArray(), INSERT_VALUES);
                PetscInt index3[3] = {from3+3*idx, from3+3*idx+1, from3+3*idx+2};
                VecSetValues(snVec, 3, index3, sn.GetArray(), INSERT_VALUES);
            }
            
            VecAssemblyBegin(fVec);  VecAssemblyEnd(fVec);
            VecAssemblyBegin(sVec);  VecAssemblyEnd(sVec);
            VecAssemblyBegin(snVec); VecAssemblyEnd(snVec);
            std::string qString = std::to_string(q+1);
            ExportElementsVectorData("TetgenBasesFiberAtQuadraturePoint" + qString, fVec);
            ExportElementsVectorData("TetgenBasesSheetAtQuadraturePoint" + qString, sVec);
            ExportElementsVectorData("TetgenBasesSheetnormalAtQuadraturePoint" + qString, snVec);
            VecDestroy(&fVec);
            VecDestroy(&sVec);
            VecDestroy(&snVec);
        }
    }
} // CBSolver::ExportBasesAtAllQuadraturePoints

void CBSolver::ExportPK2Stress() {
    if (model_->GetExporter()->GetExportOption("PK2Stress", true)) {
        Vec diagonalElements;
        Vec offdiagonalElementes;
        
        PetscInt from, to;
        
        if (DCCtrl::IsParallel()) {
            VecCreateMPI(Petsc::Comm(), 3 * numLocalElements_, PETSC_DETERMINE, &diagonalElements);
        } else {
            VecCreateSeq(PETSC_COMM_SELF, 3*numElements_, &diagonalElements);
        }
        
        VecGetOwnershipRange(diagonalElements, &from, &to);
        
        VecDuplicate(diagonalElements, &offdiagonalElementes);
        
        VecZeroEntries(diagonalElements);
        VecZeroEntries(offdiagonalElementes);
        VecSet(diagonalElements, 0);
        VecSet(offdiagonalElementes, 0);
        
        for (auto &it : solidElements_) {
            Matrix3f S;
            S = it->GetPK2Stress();
            
            int i = it->GetLocalIndex();
            PetscScalar diag[3] = {S(0, 0), S(1, 1), S(2, 2)};
            PetscScalar off[3] = {S(0, 1), S(0, 2), S(1, 2)};
            
            PetscInt index[3] = {from+3*i, from +3*i+1, from +3*i+2};
            
            VecSetValues(diagonalElements, 3, index, diag, INSERT_VALUES);
            VecSetValues(offdiagonalElementes,  3, index, off,  INSERT_VALUES);
        }
        
        VecAssemblyBegin(diagonalElements);
        VecAssemblyBegin(offdiagonalElementes);
        VecAssemblyEnd(diagonalElements);
        VecAssemblyEnd(offdiagonalElementes);
        
        ExportElementsVectorData("S_ii", diagonalElements);
        ExportElementsVectorData("S_ij", offdiagonalElementes);
        
        VecDestroy(&diagonalElements);
        VecDestroy(&offdiagonalElementes);
    }
} // CBSolver::ExportPK2Stress

void CBSolver::ExportGreenLagrangeStrain() {
    if (model_->GetExporter()->GetExportOption("GreenLagrangeStrain", true)) {
        Vec diagonalElements;
        Vec offdiagonalElementes;
        
        PetscInt from, to;
        
        if (DCCtrl::IsParallel()) {
            VecCreateMPI(Petsc::Comm(), 3 * numLocalElements_, PETSC_DETERMINE, &diagonalElements);
        } else {
            VecCreateSeq(PETSC_COMM_SELF, 3*numElements_, &diagonalElements);
        }
        
        VecGetOwnershipRange(diagonalElements, &from, &to);
        
        VecDuplicate(diagonalElements, &offdiagonalElementes);
        
        VecZeroEntries(diagonalElements);
        VecZeroEntries(offdiagonalElementes);
        VecSet(diagonalElements, 1);
        VecSet(offdiagonalElementes, 1);
        
        for (auto &it : solidElements_) {
            Matrix3f F, C, E, I;
            it->GetDeformationTensor(F);
            I = Matrix3<TFloat>::Identity();
            C = F.GetTranspose() * F;
            E = 0.5 * (C - I);
            
            int i = it->GetLocalIndex();
            PetscScalar diag[3] = {E(0, 0), E(1, 1), E(2, 2)};
            PetscScalar off[3] = {E(0, 1), E(0, 2), E(1, 2)};
            
            PetscInt index[3] = {from+3*i, from +3*i+1, from +3*i+2};
            
            VecSetValues(diagonalElements, 3, index, diag, INSERT_VALUES);
            VecSetValues(offdiagonalElementes,  3, index, off,  INSERT_VALUES);
        }
        
        VecAssemblyBegin(diagonalElements);
        VecAssemblyBegin(offdiagonalElementes);
        VecAssemblyEnd(diagonalElements);
        VecAssemblyEnd(offdiagonalElementes);
        
        ExportElementsVectorData("E_ii", diagonalElements);
        ExportElementsVectorData("E_ij", offdiagonalElementes);
        
        VecDestroy(&diagonalElements);
        VecDestroy(&offdiagonalElementes);
    }
} // CBSolver::ExportGreenLagrangeStrain

void CBSolver::ExportDeformation() {
    if (model_->GetExporter()->GetExportOption("Deformation", false)) {
        Vec deformDiag;
        Vec deformOff;
        
        PetscInt from, to;
        
        if (DCCtrl::IsParallel()) {
            VecCreateMPI(Petsc::Comm(), 3 * numLocalElements_, PETSC_DETERMINE, &deformDiag);
        } else {
            VecCreateSeq(PETSC_COMM_SELF, 3*numElements_, &deformDiag);
        }
        
        VecGetOwnershipRange(deformDiag, &from, &to);
        
        VecDuplicate(deformDiag, &deformOff);
        
        VecZeroEntries(deformDiag);
        VecZeroEntries(deformOff);
        VecSet(deformDiag, 1);
        VecSet(deformOff, 1);
        
        for (auto &it : solidElements_) {
            Matrix3f F;
            it->GetDeformationTensor(F);
            
            F = F.GetTranspose() * F;
            Matrix3f &C = F;
            
            int i = it->GetLocalIndex();
            
            PetscScalar diag[3] = {C(0, 0), C(1, 1), C(2, 2)};
            PetscScalar off[3] = {C(0, 1), C(0, 2), C(1, 2)};
            
            PetscInt index[3] = {from+3*i, from +3*i+1, from +3*i+2};
            
            VecSetValues(deformDiag, 3, index, diag, INSERT_VALUES);
            VecSetValues(deformOff,  3, index, off,  INSERT_VALUES);
        }
        
        VecAssemblyBegin(deformDiag);
        VecAssemblyBegin(deformOff);
        VecAssemblyEnd(deformDiag);
        VecAssemblyEnd(deformOff);
        
        ExportElementsVectorData("C_ii", deformDiag);
        ExportElementsVectorData("C_ij", deformOff);
        
        VecDestroy(&deformDiag);
        VecDestroy(&deformOff);
    }
} // CBSolver::ExportDeformation

void CBSolver::ExportLambda() {
    if (model_->GetExporter()->GetExportOption("Lambda", true)) {
        Vec lambda;
        
        PetscInt from, to;
        
        if (DCCtrl::IsParallel()) {
            VecCreateMPI(Petsc::Comm(), 3 * numLocalElements_, PETSC_DETERMINE, &lambda);
        } else {
            VecCreateSeq(PETSC_COMM_SELF, 3*numElements_, &lambda);
        }
        
        VecGetOwnershipRange(lambda, &from, &to);
        
        VecSet(lambda, 1);
        
        for (auto &it : solidElements_) {
            Vector3<TFloat> l;
            
            Matrix3f F;
            it->GetDeformationTensor(F);
            
            for (int i = 0; i < 3; i++)
                l(i) = F.GetCol(i).Norm();
            
            int i = it->GetLocalIndex();
            
            PetscScalar lam[3] = {l(0), l(1), l(2)};
            
            PetscInt index[3] = {from+3*i, from +3*i+1, from +3*i+2};
            
            VecSetValues(lambda,     3, index, lam,  INSERT_VALUES);
        }
        
        VecAssemblyBegin(lambda);
        VecAssemblyEnd(lambda);
        
        ExportElementsVectorData("Lambda", lambda);
        
        VecDestroy(&lambda);
    }
} // CBSolver::ExportLambda

/// Method to write the local activation times in a file (mainly for debugging purposes).

/* @author Armin Mueller
 * @return void
 */
void CBSolver::writeLATToFile(const std::vector<int> &indices, const std::vector<double> &values) {
    std::ofstream file;
    
    file.open("lat-output.txt", std::ios::out | std::ios::trunc);
    
    // info text which leads the document
    file << "This file contains the IDs and the internally used values of the local activation "
    << "time (LAT) for each element. The LATs contain also the LAT-offset given in the input xml file.\n" <<
    std::endl;
    
    int size = indices.size();
    file << "Number of Elements: " << size << "\n";
    for (int i = 0; i < size; i++) {
        // adding 1 shifts the ids so that they match the ids in the input lat - file (starting with id 1 instead of id 0)
        file << (indices.at(i) + 1) << " " << values.at(i) << "\n";
    }
    
    file.close();
}

/// Method to export the local activation times to view them later

/* @author Armin Mueller
 * @return void
 */
void CBSolver::ExportLocalActivationTime() {
    // Export LAT values only in the first timestep. Since they won't change during simulation this should be sufficient.
    if (model_->GetExporter()->GetExportOption("LocalActivationTime", false) && (model_->GetCurrentTime() == 0.0)) {
        // save the LAT indices and the LAT values (aka the start times) separately. But latIndices[i] holds the index of latValues[i]
        std::vector<int> latIndices;
        std::vector<double> latValues;
        
        // Create the Petsc vector to store the LAT values in (it is needed for the export-method)
        Vec localActivationTime;
        if (DCCtrl::IsParallel()) {
            VecCreateMPI(Petsc::Comm(), numLocalElements_, PETSC_DETERMINE, &localActivationTime);
        } else {
            VecCreateSeq(PETSC_COMM_SELF, numElements_, &localActivationTime);
        }
        
        // get the element ranges of each thread from the Petsc LAT vector
        PetscInt latLowerIndex, latUpperIndex;
        VecGetOwnershipRange(localActivationTime, &latLowerIndex, &latUpperIndex);
        VecZeroEntries(localActivationTime);
        
        // fill both std-vectors with their data
        for (auto &it : solidElements_) {
            // adding latLowerIndex to the local index results in the correct IDs since the elements are distributed over all threads
            // (normal sorting doesn't work here)
            int elementIndex = it->GetLocalIndex() + latLowerIndex;
            double lat = it->GetTensionModel()->GetActivationTime();
            
            latIndices.push_back(elementIndex);
            latValues.push_back(lat);
        }
        
        // sync local std-vector to become global std-vectors that contain all information from all and everything
        if (DCCtrl::IsParallel()) {
            // collect all lat values
            DCCtrl::GatherToZero(latIndices);
            DCCtrl::GatherToZero(latValues);
        }
        
        // copy values to the Petsc Vector, since the export-method can only work on Petsc vectors
        VecSetValues(localActivationTime, latIndices.size(), latIndices.data(), latValues.data(), INSERT_VALUES);
        
        
        // assemble the vector, export it and then destroy it to retrieve the occupied space
        // (assembly isn't compulsorily needed here, since it is also done by the ExportElementsScalarData-method)
        VecAssemblyBegin(localActivationTime);
        VecAssemblyEnd(localActivationTime);
        ExportElementsScalarData("LocalActivationTime", localActivationTime);
        VecDestroy(&localActivationTime);
        
        // if debugging is on and the LAT reader is really used (proxy-condition: if there is a LAT-file given), write the internal used lat values to a txt-file
        if (DCCtrl::GetDebug() && !parameters_->Get<std::string>("Mesh.LAT.FilePath", "").empty()) {
            DCCtrl::debug << "Writing LAT-debug file (\"lat-output.txt\") ...";
            writeLATToFile(latIndices, latValues);
            DCCtrl::debug << " DONE!\n" << std::endl;
        }
    }
} // CBSolver::ExportLocalActivationTime

void CBSolver::DeInit() {
    if (formulation_ != 0)
        delete formulation_;
    if (activeStressData_ != 0)
        delete activeStressData_;
    
    if (DCCtrl::IsParallel())
        DeInitExporter();
    if (adapter_ != 0) {
        delete adapter_;
        adapter_ = 0;
    }
    
    if (globalNodesComponentsBoundaryConditions_ != 0) {
        delete[] globalNodesComponentsBoundaryConditions_;
        globalNodesComponentsBoundaryConditions_ = 0;
    }
    
    if (activeStressTensorComponents_) {
        delete[] activeStressTensorComponents_;
        activeStressTensorComponents_ = 0;
    }
    if (activeStressTensorComponentsIndices_) {
        delete[] activeStressTensorComponentsIndices_;
        activeStressTensorComponentsIndices_ = 0;
    }
    if (activeStressTensorIndices_) {
        delete[] activeStressTensorIndices_;
        activeStressTensorComponentsIndices_ = 0;
    }
    
    MatDestroy(&boundaryConditionsNodalForcesJacobianDiagonalComponents_);
    VecDestroy(&nodes_);
    
    if (nodesSeq_)
        VecDestroy(&nodesSeq_);
    
    if (refNodesSeq_)
        VecDestroy(&refNodesSeq_);
    
    VecDestroy(&activeStress_);
    VecDestroy(&nodalForces_);
    if (corruptElements_ != 0)
        VecDestroy(&corruptElements_);
} // CBSolver::DeInit

void CBSolver::LoadMesh() {
    numElements_ = model_->GetElements().size();
    numNodes_    = model_->GetNodes().size();
    std::cout.width(10);
    
    DCCtrl::print <<
    "\xd\t\tInitializing: Mesh: Determining nodes ranges...                                                ";
    DetermineNodesRanges();
    
    DCCtrl::print <<
    "\xd\t\tInitializing: Mesh: Elements ...                                                               ";
    InitElements();
    
    DCCtrl::print <<
    "\xd\t\tInitializing: Mesh: Extracting solid elements ...                                              ";
    ExtractSolidElements();
    
    DCCtrl::print <<
    "\xd\t\tInitializing: Mesh: Extracting surface elements ...                                            ";
    ExtractSurfaceElements();
    
    DCCtrl::print <<
    "\xd\t\tInitializing: Mesh: Initializing materials ...                                                 ";
    InitMaterials();
    
    DCCtrl::print <<
    "\xd\t\tInitializing: Mesh: Nodes indices mapping ...                                                  ";
    InitNodesIndicesMapping();
    InitNodesIndicesMappingNonGhosted();
    
    DCCtrl::print <<
    "\xd\t\tInitializing: Mesh: Nodes ...                                                                  ";
    InitNodes();
    model_->ResetRefNodes();
    
    DCCtrl::print <<
    "\xd\t\tInitializing: Mesh: Boundary conditions ...                                                    ";
    InitNodesComponentsBoundaryConditions();
    
    DCCtrl::print <<
    "\xd\t\tInitializing: Mesh: Preparing elements ...                                                     ";
    PrepareElements();
    
    DCCtrl::print <<
    "\xd\t\tInitializing: Mesh: Calculating NNZ ...                                                        ";
    model_->CalcNonZeroes();
    
    // HERE LOAD MESH STATE AFTER THE PREPARE ELEMENTS, WHERE THE SHAPE FUNCTIONS ARE CALCULATED
    if (parameters_->IsAvailable("MeshLoadedState.Format")) {
        DCCtrl::print <<
        "\xd\t\tInitializing: Mesh: Init Loaded State ...                                                      ";
        InitLoadedState();
    }
} // CBSolver::LoadMesh

void CBSolver::InitLoadedState() {
    CBModelLoader *modelLoader;
    std::string Tag = "MeshLoadedState.Format";
    std::string format = parameters_->Get<std::string>(Tag);
    
    if (format == "Tetgen")
        modelLoader = new CBModelLoaderTetgen(parameters_); // deleted in deconstructor
    else
        throw std::runtime_error("Mesh Format: " + format + " is unkown !");
    LoadedModel_ = modelLoader->Load("MeshLoadedState");
    
    if (LoadedModel_->GetNodes().size() != model_->GetNodes().size())
        throw std::runtime_error("Number of nodes in the loaded model is not the same as in the unloaded model!");
    if (LoadedModel_->GetElements().size() != model_->GetElements().size())
        throw std::runtime_error("Number of elements in the loaded model is not the same as in the unloaded model!");
    
    InitNodes(LoadedModel_, false);
    if (parameters_->Get<bool>(Tag+".VTK.ResetStresses", false)) {
        RelaxElementsAndBases();
        DCCtrl::print<<"RESET STRESSES\n";
    }
}

void CBSolver::DetermineNodesRanges() {
    nodesRanges_ = model_->GetNodesRanges();
    
    if (nodesRanges_.size() == 0) {
        /* Determine the distribution of the nodes to the processes
         *
         * Example: 19 nodes, 4 processes
         *
         * nodesRanges[0] = 0
         * nodesRanges[1] = 4
         * nodesRanges[2] = 9
         * nodesRanges[3] = 14
         * nodesRanges[4] = 19
         *
         * process 1:  0  1  2  3
         * process 2:  4  5  6  7  8
         * process 3:  9 10 11 12 13
         * process 4: 14 15 16 17 18
         *
         */
        nodesRanges_.resize(DCCtrl::GetNumberOfProcesses() + 1);
        nodesRanges_[0] = 0;
        PetscInt n = numNodes_;
        for (unsigned int i = 1; i < DCCtrl::GetNumberOfProcesses(); i++) {
            PetscInt m = n / ((DCCtrl::GetNumberOfProcesses() - i) + 1);
            nodesRanges_[i] = nodesRanges_[i - 1] + m;
            n -= m;
        }
        
        nodesRanges_[DCCtrl::GetNumberOfProcesses()] = numNodes_;
    }
    
    localNodesFrom_        = nodesRanges_[DCCtrl::GetProcessID()];                // Example: process 2: localNodesFrom = 4
    localNodesTo_          = nodesRanges_[DCCtrl::GetProcessID() + 1] - 1;        //                                       localNodesTo = 9-1 = 8
    numLocalNodes_         = (localNodesTo_ - localNodesFrom_) + 1; //           numLocalNodes = 8 - 4 + 1 = 5
}

/*!
 * Distributes the elements from the model object and distributes
 * them to the processes.
 * Furthermore the function return a std::set<PetscInt> which contains the needed ghost nodes.
 * This set is used later on to generate the PETSc ghost vector
 */
void CBSolver::InitElements() {
    // Determine for each element the process [bestProcess] on which most of its nodes are located.
    if (nodesRanges_.size() == 0)
        throw std::runtime_error("CBSolver::DetermineNodesRanges() has to be run before CBSolver::InitElements() ");
    
    std::vector<unsigned int> processHistogramm;
    
    std::set<PetscInt>       ghostNodesIndicesSet; // A set is used for the ghost nodes to prevent duplicates
    
    for (PetscInt i = 0; i < numElements_; i++) {
        processHistogramm.assign(DCCtrl::GetNumberOfProcesses(), 0);
        
        CBElement *element = model_->GetElements().at(i);
        
        for (unsigned int j = 0; j < element->GetNumberOfNodesIndices(); j++) {
            PetscInt node = element->GetNodeIndex(j);
            for (unsigned int k = 0; k < DCCtrl::GetNumberOfProcesses(); k++) {
                if ((node >= nodesRanges_[k]) && (node < nodesRanges_[k + 1]))
                    processHistogramm[k]++;
            }
        }
        unsigned int maxValue    = 0;
        unsigned int bestProcess = 0;
        for (unsigned int j = 0; j < DCCtrl::GetNumberOfProcesses(); j++) {
            if (processHistogramm[j] > maxValue) {
                maxValue    = processHistogramm[j];
                bestProcess = j;
            }
        }
        
        // Determine the ghost nodes for the element which are located on the other processes and store the element
        if (bestProcess == DCCtrl::GetProcessID()) {
            for (unsigned int j = 0; j < element->GetNumberOfNodesIndices(); j++) {
                PetscInt node = element->GetNodeIndex(j);
                if ((node < localNodesFrom_) || (node > localNodesTo_) ) {
                    ghostNodesIndicesSet.insert(node);
                }
            }
            
            //            CBElement* solverElement = elementFactory_.New(element->GetType());
            //            solverElement->SetIndex(element->GetIndex());
            //            PetscInt matIndex = element->GetMaterialIndex();
            //            solverElement->SetMaterialIndex(matIndex);
            //            solverElement->SetBasis(element->GetBasis()->GetConvertedMatrix<PetscScalar>());
            //
            //            for(unsigned int k = 0;k < element->GetNumberOfNodesIndices();k++)
            //            {
            //                PetscInt index = element->GetNodeIndex(k);
            //                solverElement->SetNodeIndex(k, index);
            //            }
            
            CBElement *solverElement = element->Clone();
            solverElement->SetLocalIndex(elements_.size());
            solverElement->SetAdapter(adapter_);
            
            // Set Indices to global indices !!! //
            
            
            elements_.push_back(solverElement);
        }
    }
    
    // Convert elements global indices to local indices !!!!
    
    for (std::vector<CBElement *>::iterator element = elements_.begin(); element != elements_.end(); element++)
        for (unsigned int i = 0; i < (*element)->GetNumberOfNodesIndices(); i++) {
            PetscInt index = (*element)->GetNodeIndex(i);
            
            /*
             * #warning Optimize !!! -> only use if index is outside local nodes range !!!!
             *
             *  std::set<PetscInt>::iterator it = ghostNodesIndicesSet.find(index);
             *
             * if( it == ghostNodesIndicesSet.end())
             * (*element)->SetNodeIndex(i, ( (*element)->GetNodeIndex(i) - localNodesFrom));
             * else
             * {
             * (*element)->SetNodeIndex(i, (numLocalNodes + distance(ghostNodesIndicesSet.begin(),it)));
             * }
             */
            if ((index <= localNodesTo_) && (index >= localNodesFrom_) ) {
                (*element)->SetNodeIndex(i, (index - localNodesFrom_));
            } else {
                std::set<PetscInt>::iterator it = ghostNodesIndicesSet.find(index);
                if (it == ghostNodesIndicesSet.end()) {
                    throw std::runtime_error("Big problem with the ghostnodes -> Fix it, now !!!!");
                } else {
                    (*element)->SetNodeIndex(i, (numLocalNodes_ + distance(ghostNodesIndicesSet.begin(), it)));
                }
            }
        }
    
    numLocalElements_ = elements_.size();
    numGhostNodes_ = ghostNodesIndicesSet.size();
    ghostNodes_ = new PetscInt[numGhostNodes_];
    
    PetscInt i = 0;
    for (std::set<PetscInt>::iterator it = ghostNodesIndicesSet.begin(); it != ghostNodesIndicesSet.end(); it++) {
        ghostNodes_[i] = *it;
        i++;
    }
    isInitElementsDone_ = true;
} // CBSolver::InitElements

void CBSolver::InitNodes(CBModel *fromModel, bool isRefNode) {
    if (fromModel == nullptr)
        fromModel = model_;
    
    if (!isInitElementsDone_)
        throw std::runtime_error("CBSolver::InitElements() has to be run before CBSolver::InitNodes()");
    
    // Use VecCreateGhost to create a new PETSc vector and load local nodes to it;
    
    if (DCCtrl::IsParallel()) {
        PetscInt *ghostNodesCoordsIndices = new PetscInt[3 * numGhostNodes_];
        for (PetscInt i = 0; i < numGhostNodes_; i++) {
            PetscInt n = ghostNodes_[i];
            
            ghostNodesCoordsIndices[3 * i]     = 3 * n;
            ghostNodesCoordsIndices[3 * i + 1] = 3 * n + 1;
            ghostNodesCoordsIndices[3 * i + 2] = 3 * n + 2;
        }
        VecCreateGhost(Petsc::Comm(), 3 * numLocalNodes_, PETSC_DECIDE, 3 * numGhostNodes_, ghostNodesCoordsIndices,
                       &nodes_);
        VecCreateGhost(
                       Petsc::Comm(), 3 * numLocalNodes_, PETSC_DECIDE, 3 * numGhostNodes_, ghostNodesCoordsIndices, &refNodes_);
        
        // delete[] ghostNodesCoordsIndices;
    } else {
        VecCreateSeq(Petsc::Comm(), 3 * numNodes_, &nodes_);
        VecCreateSeq(Petsc::Comm(), 3 * numNodes_, &refNodes_);
    }
    
    
    for (PetscInt i = localNodesFrom_; i <= localNodesTo_; i++) {
        PetscInt pos[3];
        
        pos[0] = 3 * i;
        pos[1] = 3 * i + 1;
        pos[2] = 3 * i + 2;
        
        PetscScalar          val[3];
        
        Vector3<PetscScalar> vec = fromModel->GetNodes().at(i);
        val[0] = (PetscScalar)vec.X();
        val[1] = (PetscScalar)vec.Y();
        val[2] = (PetscScalar)vec.Z();
        
        VecSetValues(nodes_, 3, pos, val, INSERT_VALUES);
        VecSetValues(refNodes_, 3, pos, val, INSERT_VALUES);
    }
    
    VecAssemblyBegin(nodes_);
    VecAssemblyEnd(nodes_);
    VecAssemblyBegin(refNodes_);
    VecAssemblyEnd(refNodes_);
    
    if (DCCtrl::IsParallel()) {
        VecGhostUpdateBegin(nodes_, INSERT_VALUES, SCATTER_FORWARD);
        VecGhostUpdateEnd(nodes_, INSERT_VALUES, SCATTER_FORWARD);
        
        if (isRefNode) {
            VecGhostUpdateBegin(refNodes_, INSERT_VALUES, SCATTER_FORWARD);
            VecGhostUpdateEnd(refNodes_, INSERT_VALUES, SCATTER_FORWARD);
        }
    }
    
    adapter_->LinkRefNodes(refNodes_);
    
    isInitNodesDone_ = true;
} // CBSolver::InitNodes

void CBSolver::InitNodesComponentsBoundaryConditions() {
    if (!isInitNodesDone_) {
        throw std::runtime_error(
                                 "CBSolver::InitNodes() has to be run before CBSolver::InitNodesComponentsBoundaryConditions()");
    }
    
    globalNodesComponentsBoundaryConditions_ = new bool[3 * numNodes_];
    for (PetscInt i = 0; i < numNodes_; i++) {
        PetscInt bc = model_->GetNodeBoundaryConditions().at(i);
        
        if (bc & 1)
            globalNodesComponentsBoundaryConditions_[3 * i] = 1;
        else
            globalNodesComponentsBoundaryConditions_[3 * i] = 0;
        
        
        if (bc & 2)
            globalNodesComponentsBoundaryConditions_[3 * i + 1] = 1;
        else
            globalNodesComponentsBoundaryConditions_[3 * i + 1] = 0;
        
        if (bc & 4)
            globalNodesComponentsBoundaryConditions_[3 * i + 2] = 1;
        else
            globalNodesComponentsBoundaryConditions_[3 * i + 2] = 0;
    }
    
    if (DCCtrl::IsParallel()) {
        MatCreateAIJ(
                     Petsc::Comm(), 3 * numLocalNodes_, 3 * numLocalNodes_, PETSC_DETERMINE, PETSC_DETERMINE, 3, PETSC_NULL, 1, PETSC_NULL,
                     &boundaryConditionsNodalForcesJacobianDiagonalComponents_);
        MatSetLocalToGlobalMapping(boundaryConditionsNodalForcesJacobianDiagonalComponents_, nodesIndicesMapping_,
                                   nodesIndicesMapping_);
    } else {
        MatCreateSeqAIJ(
                        Petsc::Comm(), 3 * numNodes_, 3 * numNodes_, 3, PETSC_NULL,
                        &boundaryConditionsNodalForcesJacobianDiagonalComponents_);
        MatSetLocalToGlobalMapping(boundaryConditionsNodalForcesJacobianDiagonalComponents_, nodesIndicesMapping_,
                                   nodesIndicesMapping_);
    }
    
    adapter_->LinkNodesComponentsBoundaryConditionsGlobal(globalNodesComponentsBoundaryConditions_);
} // CBSolver::InitNodesComponentsBoundaryConditions

void CBSolver::InitNodesIndicesMapping() {
    if (!isInitElementsDone_)
        throw std::runtime_error("CBSolver::InitElements() has to be run before CBSolver::InitNodesIndicesMapping()");
    
    PetscInt  numTotalNodes = numLocalNodes_ + numGhostNodes_;
    
    PetscInt *globalNodesCoordsIndices = new PetscInt[3 * numTotalNodes];
    
    for (PetscInt i = 0; i < numTotalNodes; i++) {
        if (i < numLocalNodes_) {
            PetscInt n = (localNodesFrom_ + i);
            globalNodesCoordsIndices[3 * i]     = 3 * n;
            globalNodesCoordsIndices[3 * i + 1] = 3 * n + 1;
            globalNodesCoordsIndices[3 * i + 2] = 3 * n + 2;
        } else {
            PetscInt n = (ghostNodes_[i - numLocalNodes_]);
            globalNodesCoordsIndices[3 * i]     = 3 * n;
            globalNodesCoordsIndices[3 * i + 1] = 3 * n + 1;
            globalNodesCoordsIndices[3 * i + 2] = 3 * n + 2;
        }
    }
    
    ISLocalToGlobalMappingCreate(
                                 Petsc::Comm(), 1, 3 * numTotalNodes, globalNodesCoordsIndices, PETSC_COPY_VALUES, &nodesIndicesMapping_);
    
    adapter_->LinkLocalToGlobalMapping(nodesIndicesMapping_);
    delete[] globalNodesCoordsIndices;
} // CBSolver::InitNodesIndicesMapping

void CBSolver::InitNodesIndicesMappingNonGhosted() {
    if (!isInitElementsDone_) {
        throw std::runtime_error(
                                 "CBSolver::InitElements() has to be run before CBSolver::InitNodesIndicesMappingNonGhosted()");
    }
    
    PetscInt *globalNodesCoordsIndices = new PetscInt[3 * numLocalNodes_];
    
    for (PetscInt i = 0; i < numLocalNodes_; i++) {
        PetscInt n = (localNodesFrom_ + i);
        globalNodesCoordsIndices[3 * i]     = 3 * n;
        globalNodesCoordsIndices[3 * i + 1] = 3 * n + 1;
        globalNodesCoordsIndices[3 * i + 2] = 3 * n + 2;
    }
    
    ISLocalToGlobalMappingCreate(
                                 Petsc::Comm(), 1, 3 * numLocalNodes_, globalNodesCoordsIndices, PETSC_COPY_VALUES, &nodesIndicesMappingNonGhosted_);
    
    adapter_->LinkLocalToGlobalMappingNonGhosted(nodesIndicesMappingNonGhosted_);
    delete[] globalNodesCoordsIndices;
}

PetscScalar CBSolver::CalcFiniteDifferencesEpsilon(Vec a) {
    // adapted form snesj.c [petsc-3.3-p1]
    
    PetscReal      dx_min = 1.e-16,
    dx_par = 1.e-1;
    PetscScalar unorm;
    
    VecNorm(a, NORM_2, &unorm);
    PetscScalar dx = 1.0 + unorm;
    if (PetscAbsScalar(dx) < dx_min)
        dx = (PetscRealPart(dx) < 0. ? -1. : 1.) * dx_par;
    
    dx *= epsilon_;
    
    return dx;
}

void CBSolver::CreateNodesJacobianAndLinkToAdapter() {
    if (nodalForcesJacobian_ != 0)
        MatDestroy(&nodalForcesJacobian_);
    
    if (DCCtrl::IsParallel()) {
        MatCreateAIJ(
                     Petsc::Comm(), 3 * numLocalNodes_, 3 * numLocalNodes_, PETSC_DETERMINE, PETSC_DETERMINE, 0,
                     model_->GetNodeNeighborsForNnz().data() + localNodesFrom_*3, 0,
                     model_->GetNodeNeighborsForNnz().data() + localNodesFrom_*3, &nodalForcesJacobian_);
        MatSetLocalToGlobalMapping(nodalForcesJacobian_, nodesIndicesMapping_, nodesIndicesMapping_);
    } else {
        PetscErrorCode ierr = MatCreateSeqAIJ(Petsc::Comm(), 3 * numNodes_, 3 * numNodes_, 0,
                                              model_->GetNodeNeighborsForNnz().data(), &nodalForcesJacobian_);
        ierr = MatSetLocalToGlobalMapping(nodalForcesJacobian_, nodesIndicesMapping_, nodesIndicesMapping_);
    }
    adapter_->LinkNodalForcesJacobian(nodalForcesJacobian_);
}

CBStatus CBSolver::PrepareSimulation() {
    status_ = CBStatus::PREPARING_SIMULATION;
    auto p = plugins_.begin();
    bool initializePlugins = true;
    int initStep = 0;
    
    std::string initializingPlugin;
    DCCtrl::print << "\n\n\n-------------------------- Initializing Simulation... --------------------------\n\n";
    TFloat tt = 0;
    timing_.SetCurrentTime(0);
    timing_.ResetForInit();
    
    while (initializePlugins) {
        DCCtrl::verbose << "Initialization Step " << initStep << " [time step: " << timing_.GetTimeStep() <<"]\n";
        CBStatus rc = CBStatus::NOTHING_DONE;
        
        CreateNodesJacobianAndLinkToAdapter(); // Is this needed by the plugins? [lb451]
        
        while (p != plugins_.end()) {
            if ((*p)->GetStatus() == CBStatus::WAITING)
                (*p)->Prepare();
            
            if ((*p)->GetStatus() == CBStatus::PREPARING_SIMULATION) {
                initializingPlugin = (*p)->GetName();
                break;
            } else if ((*p)->GetStatus() == CBStatus::DACCORD) {
                DCCtrl::print << "\nPlugin " << (*p)->GetName() << " has been initialized successfully \n";
                initStep = 0;
                p++;
            } else {
                throw std::runtime_error("CBSolver::Run(): Initialization of the plugin [" + (*p)->GetName() + "] has failed");
            }
        }
        if (p == plugins_.end()) {
            DCCtrl::print << "\nInitialization finished " << "\n";
            initializePlugins = false;
        }
        
        // === init ===============================
        UpdateActiveStress(0);
        rc = SolverStep(0, true);
        timing_.SetSolverStatus(rc);
        UpdateGhostNodesAndLinkToAdapter();
        
        // ========================================
        
        if (rc != CBStatus::SUCCESS) {
            for (auto &it : plugins_)
                it->StepBack();
            
            tt -= timing_.GetTimeStep();
        }
        
        tt += timing_.GetTimeStep();
        if (timing_.ComputeNextTimeStep() != CBStatus::SUCCESS) {
            status_ = CBStatus::FAILED;
            shallExit_ = true;
            break;
        }
        
        if (initializePlugins) {
            if ((*p)->GetPreparationProgress() * 100.0 < 99.9) {
                DCCtrl::print << "\xdPlugin: " << initializingPlugin << " - Initialization progress " <<
                (*p)->GetPreparationProgress() * 100.0<< "%\n";
            } else {
                DCCtrl::print <<"\xdPlugin: " << initializingPlugin << " - Initialization progress almost finished \n";
            }
            DCCtrl::verbose << "\n";
            DCCtrl::debug << "\n";
            initStep++;
        }
    }
    
    // TODO: write as own plugin that can be turned on or off
    SetZeroVelocityAndAcceleration();
    status_ = CBStatus::DACCORD;
    return status_;
} // CBSolver::PrepareSimulation

void CBSolver::Run() {
    if (!(status_ == CBStatus::WAITING)) {
        throw std::runtime_error(
                                 "CBSolver::Run(): void CBSolver::Init(ParameterMap* parameters, CBModel* model,CBModelExporter* modelExporter) has to be invoked first !!!");
    }
    if (!model_->HasExporter())
        throw std::runtime_error("CBSolver::Run(): No Exporter initialized !");
    
    status_ = CBStatus::RUNNING;
    
    CBRuntimeEstimator runtimeEstimator(timing_);
    runtimeEstimator.EnableOutput(enableTiming_); // disable runtime output, can be useful for debugging logfiles using diff
    PetscScalar lastExportTime      = timing_.GetStartTime();
    PetscScalar time                = timing_.GetStartTime();
    bool        isFirstStep         = true;
    
    PrepareSimulation();
    
    DCCtrl::print << "\n\n\n---------------------------- Running Simulation... -----------------------------\n\n";
    timing_.ResetForSim();
    runtimeEstimator.Tic();
    
    while (!shallExit_ && timing_.IsValid(time)) { // 1e-14 prevents more runs after 100% done
        CBStatus rc = CBStatus::NOTHING_DONE;
        
        bool repeat = true;
        while (repeat) { // this loop iterates until one timestep could be computed successfully
            CreateNodesJacobianAndLinkToAdapter();
            
            // the variable 'time' contains always the last successful (!) timestep
            // the variable 't' contains the next guessed time where the problem will be tried to be solved
            
            PetscScalar t = time + timing_.GetTimeStep();
            timing_.SetCurrentTime(t);
            
            if (isFirstStep) {
                time = timing_.GetStartTime() - timing_.GetTimeStep();     // do the first pass with t=0, to export the initial state
                t = timing_.GetStartTime();
            }
            
            DCCtrl::verbose << "\n--- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---\n";
            DCCtrl::verbose << "\nStep: " << time << " + " << timing_.GetTimeStep() << ", [time step size: " <<
            timing_.GetTimeStep() << " s]\n\n";
            
            // === sim ================================
            UpdateActiveStress(t);
            rc = SolverStep(t);
            timing_.SetSolverStatus(rc);
            UpdateGhostNodesAndLinkToAdapter();
            
            // ========================================
            
            if (rc == CBStatus::SUCCESS) {
                repeat = false;
                status_ = CBStatus::SUCCESS;
                time = t;
            } else {
                repeat = true;
            }
            
            // reset plugins as needed
            if (repeat) {
                if (rc == CBStatus::FAILED) {
                    for (auto &p : plugins_)
                        p->StepBack();
                } else if (rc == CBStatus::REPEAT) {
                    for (auto &p : plugins_) {
                        if ((p->GetStatus() != CBStatus::REPEAT) && (p->GetName() != "acCELLerate"))
                            p->StepBack();
                    }
                }
            }
            
            // this might change the result of GetTimeStep()
            if (timing_.ComputeNextTimeStep() != CBStatus::SUCCESS) {
                status_ = CBStatus::FAILED;
                shallExit_ = true;
                break;
            }
        } // end while(repeat)
        
        
        if (isFirstStep) {
            // initial solution t=0 should not count, but found timestep should stay
            isFirstStep = false;
            timing_.ResetCounter();
        }
        
        // EXPORT STUFF //
        
        if ((time == timing_.GetStartTime()) ||
            (time - lastExportTime >= (1-1e-12)*model_->GetExporter()->GetExportTimeStep())) { // 1e-15 ~double machine epsilon
            Export(time);
            if (DCCtrl::IsProcessZero()) {
                // TODO: model should take care of how many processes can export
                DCCtrl::print << "Exporting time step " << time << " s \n";
                model_->GetExporter()->WriteToFile();
            }
            if (exportSNESMatrix_ && (time != timing_.GetStartTime())) {
                this->ExportSNESMatrix(time);
            }
            
            lastExportTime = time;
        }
        
        // Write to file from Plugins
        for (auto &p : plugins_)
            // note: process zero check should better happen inside the plugins to allow parallel writing
            if (DCCtrl::IsProcessZero()) {
                p->WriteToFile(time);
            }
        
        // check if a plugin decided that no further timesteps need to be computed
        for (auto &p : plugins_) {
            if (p->ExitCheck(time)) {
                DCCtrl::print << "Plugin " << p->GetName() << " commanded to exit." << std::endl;
                shallExit_ = true;
            }
        }
        
        runtimeEstimator.ShowEstimation(time);
    }  // end of SUCCESSFUL time steps loop
    
    if (status_ == CBStatus::SUCCESS) {
        runtimeEstimator.ShowFinalSummary();
    } else {
        Export(time);
        if (DCCtrl::IsProcessZero()) {
            DCCtrl::print << "Exporting FAILED time step \n";
            model_->GetExporter()->SetFailedSimCounter();
            model_->GetExporter()->WriteToFile();
        }
        DCCtrl::print << "\n\n\n ######## SIMULATION FAILED #########\n";
    }
} // CBSolver::Run

void CBSolver::ExportNodesVectorData(std::string key, Vec data) {
    VecAssemblyBegin(data);
    VecAssemblyEnd(data);
    
    model_->nodesVectorData.PrepareDataMap(key, math_pack::Vector3<double>(0, 0, 0));
    auto it = model_->nodesVectorData.GetData(key).begin();
    std::lock_guard<std::mutex>(model_->GetLock());
    if (DCCtrl::IsParallel()) {
        VecCopy(data, nodesVectorData_);
        VecScatterBegin(nodesVectorDataCtx_, nodesVectorData_, nodesVectorDataSeq_, INSERT_VALUES, SCATTER_FORWARD);
        VecScatterEnd(nodesVectorDataCtx_, nodesVectorData_, nodesVectorDataSeq_, INSERT_VALUES, SCATTER_FORWARD);
        
        if (DCCtrl::IsProcessZero()) {
            for (PetscInt i = 0; i < numNodes_; i++) {
                PetscInt    pos[3];
                PetscScalar values[3];
                pos[0] = 3 * i;
                pos[1] = 3 * i + 1;
                pos[2] = 3 * i + 2;
                VecGetValues(nodesVectorDataSeq_, 3, pos, values);
                Vector3<PetscScalar> n((PetscScalar)values[0], (PetscScalar)values[1], (PetscScalar)values[2]);
                *it = n;
                it++;
            }
        }
    } else {
        for (PetscInt i = 0; i < numNodes_; i++) {
            PetscInt    pos[3];
            PetscScalar values[3];
            pos[0] = 3 * i;
            pos[1] = 3 * i + 1;
            pos[2] = 3 * i + 2;
            VecGetValues(data, 3, pos, values);
            Vector3<PetscScalar> n((PetscScalar)values[0], (PetscScalar)values[1], (PetscScalar)values[2]);
            *it = n;
            it++;
        }
    }
} // CBSolver::ExportNodesVectorData

TFloat CBSolver::GetVolumeOfSolidElements() {
    TFloat modelVolume_ = 0;
    TFloat localModelVolume = 0;
    
    for (auto &it : adapter_->GetSolver()->GetSolidElementVector())
        localModelVolume += it->GetVolume();
    MPI_Allreduce(&localModelVolume, &modelVolume_, 1, MPI_DOUBLE, MPI_SUM, Petsc::Comm());
    return modelVolume_;
}

TFloat CBSolver::GetDeformationEnergy() {
    TFloat deformationEnergy_ = 0;
    TFloat localDeformationEnergy = 0;
    
    for (auto &it : adapter_->GetSolver()->GetSolidElementVector())
        localDeformationEnergy += it->GetDeformationEnergy();
    MPI_Allreduce(&localDeformationEnergy, &deformationEnergy_, 1, MPI_DOUBLE, MPI_SUM, Petsc::Comm());
    return deformationEnergy_;
}

void CBSolver::ExportNodesScalarData(std::string key, Vec data) {
    VecAssemblyBegin(data);
    VecAssemblyEnd(data);
    
    model_->nodesScalarData.PrepareDataMap(key, 0);
    auto it = model_->nodesScalarData.GetData(key).begin();
    std::lock_guard<std::mutex>(model_->GetLock());
    
    auto dataRef = data;
    if (DCCtrl::IsParallel()) {
        VecCopy(data, nodesScalarData_);
        VecScatterBegin(nodesScalarDataCtx_, nodesScalarData_, nodesScalarDataSeq_, INSERT_VALUES, SCATTER_FORWARD);
        VecScatterEnd(nodesScalarDataCtx_, nodesScalarData_, nodesScalarDataSeq_, INSERT_VALUES, SCATTER_FORWARD);
        
        if (!DCCtrl::IsProcessZero())
            return;
        
        dataRef = nodesScalarDataSeq_;
    }
    
    for (PetscInt i = 0; i < numNodes_; i++) {
        PetscScalar value;
        PetscInt    pos = i;
        VecGetValues(dataRef, 1, &pos, &value);
        *it = value;
        it++;
    }
} // CBSolver::ExportNodesScalarData

void CBSolver::ExportElementsVectorData(std::string key, Vec data) {
    VecAssemblyBegin(data);
    VecAssemblyEnd(data);
    
    auto &elementVectorData = model_->elementsVectorData;
    if (!elementVectorData.IsData(key))
        elementVectorData.InitDataMap(key, Vector3<TFloat>(0, 0, 0));
    
    auto it = elementVectorData.GetData(key).begin();
    
    std::lock_guard<std::mutex>(model_->GetLock());
    if (DCCtrl::IsParallel()) {
        VecCopy(data, elementsVectorData_);
        VecScatterBegin(elementsVectorDataCtx_, elementsVectorData_, elementsVectorDataSeq_, INSERT_VALUES,
                        SCATTER_FORWARD);
        VecScatterEnd(elementsVectorDataCtx_, elementsVectorData_, elementsVectorDataSeq_, INSERT_VALUES, SCATTER_FORWARD);
        
        if (DCCtrl::IsProcessZero()) {
            for (PetscInt i = 0; i < numElements_; i++) {
                PetscInt    pos[3];
                PetscScalar values[3];
                pos[0] = 3 * i;
                pos[1] = 3 * i + 1;
                pos[2] = 3 * i + 2;
                VecGetValues(elementsVectorDataSeq_, 3, pos, values);
                Vector3<PetscScalar>                                      n((PetscScalar)values[0], (PetscScalar)values[1],
                                                                            (PetscScalar)values[2]);
                auto it2 = it + exportElementsMapping_[i];
                *it2 = n;
            }
        }
    } else {
        for (PetscInt i = 0; i < numElements_; i++) {
            PetscInt    pos[3];
            PetscScalar values[3];
            pos[0] = 3 * i;
            pos[1] = 3 * i + 1;
            pos[2] = 3 * i + 2;
            VecGetValues(data, 3, pos, values);
            Vector3<PetscScalar> n((PetscScalar)values[0], (PetscScalar)values[1], (PetscScalar)values[2]);
            *it = n;
            it++;
        }
    }
} // CBSolver::ExportElementsVectorData

void CBSolver::ExportElementsScalarData(std::string key, Vec data) {
    PetscErrorCode ierr;
    
    ierr = VecAssemblyBegin(data); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(data); CHKERRQ(ierr);
    
    auto &elementsScalarData = model_->elementsScalarData;
    elementsScalarData.PrepareDataMap(key, 0);
    
    auto it = model_->elementsScalarData.GetData(key).begin();
    
    std::lock_guard<std::mutex>(model_->GetLock());
    if (DCCtrl::IsParallel()) {
        ierr = VecCopy(data, elementsScalarData_); CHKERRQ(ierr);
        ierr = VecScatterBegin(elementsScalarDataCtx_, elementsScalarData_, elementsScalarDataSeq_, INSERT_VALUES,
                               SCATTER_FORWARD); CHKERRQ(ierr);
        ierr = VecScatterEnd(elementsScalarDataCtx_, elementsScalarData_, elementsScalarDataSeq_, INSERT_VALUES,
                             SCATTER_FORWARD); CHKERRQ(ierr);
        
        if (DCCtrl::IsProcessZero()) {
            for (PetscInt i = 0; i < numElements_; i++) {
                PetscScalar                                               value;
                PetscInt                                                  pos = i;
                ierr = VecGetValues(elementsScalarDataSeq_, 1, &pos, &value); CHKERRQ(ierr);
                auto it2 = it + exportElementsMapping_[i];
                *it2 = value;
            }
        }
    } else {
        for (PetscInt i = 0; i < numElements_; i++) {
            PetscScalar value;
            PetscInt    pos = i;
            ierr = VecGetValues(data, 1, &pos, &value); CHKERRQ(ierr);
            *it = value;
            it++;
        }
    }
} // CBSolver::ExportElementsScalarData

void CBSolver::InitExporter() {
    if (!isInitNodesDone_ || !isInitElementsDone_)
        throw std::runtime_error("CBSolver::InitExporter(): InitNodes() and InitElements() has to be finished first");
    
    if (DCCtrl::IsParallel()) {
        VecCreateMPI(PETSC_COMM_WORLD, numLocalNodes_, PETSC_DECIDE, &nodesScalarData_);
        VecCreateSeq(PETSC_COMM_SELF, numNodes_, &nodesScalarDataSeq_);
        VecScatterCreateToZero(nodesScalarData_, &nodesScalarDataCtx_, &nodesScalarDataSeq_);
        
        VecCreateMPI(PETSC_COMM_WORLD, 3 * numLocalNodes_, PETSC_DECIDE, &nodesVectorData_);
        VecCreateSeq(PETSC_COMM_SELF, numNodes_, &nodesVectorDataSeq_);
        VecScatterCreateToZero(nodesVectorData_, &nodesVectorDataCtx_, &nodesVectorDataSeq_);
        
        VecCreateMPI(PETSC_COMM_WORLD, numLocalElements_, PETSC_DECIDE, &elementsScalarData_);
        VecCreateSeq(PETSC_COMM_SELF, numElements_, &elementsScalarDataSeq_);
        VecScatterCreateToZero(elementsScalarData_, &elementsScalarDataCtx_, &elementsScalarDataSeq_);
        
        VecCreateMPI(PETSC_COMM_WORLD, 3 * numLocalElements_, PETSC_DECIDE, &elementsVectorData_);
        VecCreateSeq(PETSC_COMM_SELF, numElements_, &elementsVectorDataSeq_);
        VecScatterCreateToZero(elementsVectorData_, &elementsVectorDataCtx_, &elementsVectorDataSeq_);
        
        PetscInt  numLocalIndices = elements_.size();
        PetscInt *localIndices    = new PetscInt[numLocalIndices];
        
        //  Elements Mapping for export. Elements have there own indices due to domain decomposition
        int       i = 0;
        
        for (auto &it : elements_) {
            localIndices[i++] = it->GetIndex();
            
            // i++;
        }
        
        PetscInt  numElementsMapping;
        const PetscInt *elementsMapping;
        
        IS localIs;
        IS globalIs;
        
        ISCreateGeneral(PETSC_COMM_WORLD, numLocalIndices, localIndices, PETSC_COPY_VALUES, &localIs);
        ISAllGather(localIs, &globalIs);
        
        ISGetSize(globalIs, &numElementsMapping);
        ISGetIndices(globalIs, &elementsMapping);
        
        
        if (DCCtrl::IsProcessZero()) {
            exportElementsMapping_ = new PetscInt[numElementsMapping];
            memcpy(exportElementsMapping_, elementsMapping, numElementsMapping * sizeof(PetscInt));
        }
        
        delete[] localIndices;
        ISDestroy(&localIs);
        ISDestroy(&globalIs);
    }
} // CBSolver::InitExporter

void CBSolver::ExportNodes() {
    std::lock_guard<std::mutex>(model_->GetLock());
    if (DCCtrl::IsParallel()) {
        VecCopy(nodes_, nodesVectorData_);
        VecScatterBegin(nodesVectorDataCtx_, nodesVectorData_, nodesVectorDataSeq_, INSERT_VALUES, SCATTER_FORWARD);
        VecScatterEnd(nodesVectorDataCtx_, nodesVectorData_, nodesVectorDataSeq_, INSERT_VALUES, SCATTER_FORWARD);
        
        if (DCCtrl::IsProcessZero()) {
            for (PetscInt i = 0; i < numNodes_; i++) {
                PetscInt    pos[3];
                PetscScalar values[3];
                pos[0] = 3 * i;
                pos[1] = 3 * i + 1;
                pos[2] = 3 * i + 2;
                VecGetValues(nodesVectorDataSeq_, 3, pos, values);
                Vector3<PetscScalar> n((PetscScalar)values[0], (PetscScalar)values[1], (PetscScalar)values[2]);
                model_->UpdateNode(i, n);
            }
        }
    } else {
        for (PetscInt i = 0; i < numNodes_; i++) {
            PetscInt    pos[3];
            PetscScalar values[3];
            pos[0] = 3 * i;
            pos[1] = 3 * i + 1;
            pos[2] = 3 * i + 2;
            VecGetValues(nodes_, 3, pos, values);
            Vector3<PetscScalar> n((PetscScalar)values[0], (PetscScalar)values[1], (PetscScalar)values[2]);
            model_->UpdateNode(i, n);
        }
    }
} // CBSolver::ExportNodes

void CBSolver::ExportBasisAtQuadraturePoint() {
    std::lock_guard<std::mutex>(model_->GetLock());
    for (int i = 0; i < numLocalElements_; i++)
        for (int j = 0; j < elements_[i]->GetNumberOfQuadraturePoints(); j++)
            model_->GetElements().at(elements_[i]->GetIndex())->SetBasisAtQuadraturePoint(j,
                                                                                          *(elements_[i]->
                                                                                            GetBasisAtQuadraturePoint(j)));
}

void CBSolver::DeInitExporter() {
    if (DCCtrl::IsParallel()) {
        VecScatterDestroy(&nodesScalarDataCtx_);
        VecScatterDestroy(&nodesVectorDataCtx_);
        VecScatterDestroy(&elementsScalarDataCtx_);
        VecScatterDestroy(&elementsVectorDataCtx_);
        
        VecDestroy(&nodesScalarData_);
        VecDestroy(&nodesScalarDataSeq_);
        VecDestroy(&nodesVectorData_);
        VecDestroy(&nodesVectorDataSeq_);
        VecDestroy(&elementsScalarData_);
        VecDestroy(&elementsScalarDataSeq_);
        VecDestroy(&elementsVectorData_);
        VecDestroy(&elementsVectorDataSeq_);
        delete[] exportElementsMapping_;
    }
}

void CBSolver::PrepareElements() {
    PrepareElements(solidElements_);
}

/// reset stress s.t. deformation tensor of each element equals I (?)
void CBSolver::PrepareElements(std::vector<CBElementSolid *> elements) {
    if (!isInitNodesDone_ || !isInitElementsDone_)
        throw std::runtime_error("CBSolver::PrepareElements(): InitNodes() and InitElements() has be finished first");
    
    UpdateGhostNodesAndLinkToAdapter();
    adapter_->LinkLocalToGlobalMapping(nodesIndicesMapping_);
    
    for (auto &it : elements) {
        it->CheckNodeSorting();
        it->UpdateShapeFunctions();
    }
}

void CBSolver::RelaxElementsAndBases() {
    RelaxElementsAndBases(solidElements_);
}

/// reset elements to reflect deformation=I and fiber basis = current deformed fibers, in contrast to prepare elements, rotation matrix is resetted to the current state as well
void CBSolver::RelaxElementsAndBases(std::vector<CBElementSolid *> elements) {
    // set current rotated basis as new basis
    for (auto &e : elements) {
        for (int i = 0; i < e->GetNumberOfQuadraturePoints(); i++) {
            Matrix3<TFloat> F;
            e->GetDeformationTensor(F);
            Matrix3<TFloat> *b = e->GetBasisAtQuadraturePoint(i);
            
            Vector3<TFloat> f  = b->GetCol(0);
            Vector3<TFloat> s  = b->GetCol(1);
            Vector3<TFloat> sn = b->GetCol(2);
            
            f  = F * f;
            s  = F * s;
            sn = F * sn;
            
            // rescale f to length = 1
            f.Normalize();
            
            // Gram Schmid orthogonalize s wrt. f
            s.Normalize();
            s = s - (f*s) * f;
            
            // sn needs to be orthogonal to the others
            sn.Normalize();
            sn = CrossProduct(f, s);
            
            Matrix3<TFloat> b2(f, s, sn);
            e->SetBasisAtQuadraturePoint(i, b2);
        }
    }
    
    // reset internal stress
    PrepareElements(elements);
} // CBSolver::RelaxElementsAndBases

void CBSolver::RelaxElementsAndBasesByMaterialIndex(TInt mat) {
    std::vector<CBElementSolid *> elements;
    
    for (auto &it : solidElements_)
        if (it->GetMaterialIndex() == mat)
            elements.push_back(it);
    
    RelaxElementsAndBases(elements);
}

void CBSolver::RelaxElementsAndBasesByIndex(TInt index) {
    std::vector<CBElementSolid *> elements;
    
    for (auto &it : solidElements_)
        if (it->GetIndex() == index)
            elements.push_back(it);
    
    RelaxElementsAndBases(elements);
}

Vector3fVector CBSolver::GetCurrentNodeCoordinates() {
    Vector3fVector retNodes;
    
    if (DCCtrl::IsParallel()) {
        VecCopy(nodes_, nodesVectorData_);
        VecScatterBegin(nodesVectorDataCtx_, nodesVectorData_, nodesVectorDataSeq_, INSERT_VALUES, SCATTER_FORWARD);
        VecScatterEnd(nodesVectorDataCtx_, nodesVectorData_, nodesVectorDataSeq_, INSERT_VALUES, SCATTER_FORWARD);
        
        if (DCCtrl::IsProcessZero()) {
            for (PetscInt i = 0; i < numNodes_; i++) {
                PetscInt    pos[3];
                PetscScalar values[3];
                pos[0] = 3 * i;
                pos[1] = 3 * i + 1;
                pos[2] = 3 * i + 2;
                VecGetValues(nodesVectorDataSeq_, 3, pos, values);
                Vector3<PetscScalar> n((PetscScalar)values[0], (PetscScalar)values[1], (PetscScalar)values[2]);
                retNodes.push_back(n);
            }
        }
    } else {
        for (PetscInt i = 0; i < numNodes_; i++) {
            PetscInt    pos[3];
            PetscScalar values[3];
            pos[0] = 3 * i;
            pos[1] = 3 * i + 1;
            pos[2] = 3 * i + 2;
            VecGetValues(nodes_, 3, pos, values);
            Vector3<PetscScalar> n((PetscScalar)values[0], (PetscScalar)values[1], (PetscScalar)values[2]);
            retNodes.push_back(n);
        }
    }
    return retNodes;
} // CBSolver::GetCurrentNodeCoordinates

void CBSolver::SetCurrentNodeCoordinates(Vector3fVector newNodes) {
    for (PetscInt i = 0; i < numNodes_; i++) {
        PetscInt    pos[3];
        PetscScalar values[3];
        pos[0] = 3 * i;
        pos[1] = 3 * i + 1;
        pos[2] = 3 * i + 2;
        values[0] = newNodes[i](0);
        values[1] = newNodes[i](1);
        values[2] = newNodes[i](2);
        VecSetValues(nodes_, 3, pos, values, INSERT_VALUES);
    }
}
