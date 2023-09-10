/*
 *  CBSolver.h
 *  CardioMechanics
 *
 *  Created by Thomas Fritz on 23.02.10.
 *  Copyright 2010 Institute for Biomedical Engineering, Karlsruhe Institute of Technology (KIT). All rights reserved.
 *
 */

#ifndef CB_SOLVER
#define CB_SOLVER

#include "CBFileManager.h"
#include "CBTensionFactory.h"
#include "CBStatus.h"
#include "CBModel.h"
#include "CBElement.h"
#include "CBElementSolid.h"
#include "CBElementFactory.h"
#include "CBElementAdapter.h"
#include "CBMaterialFactory.h"
#include "CBFormulationTotalLagrangian.h"
#include "CBSolverPlugin.h"
#include "CBSolverPluginFactory.h"
#include "CBDataPerMaterial.h"
#include "CBDataCtrl.h"
#include "CBStatus.h"
#include "CBSolver.h"
#include "CBSolverPluginFactory.h"
#include "CBTiming.h"
#include "CBData.h"

#include "ParameterMap.h"
#include "DCCtrl.h"

#include <vector>
#include <functional>
#include <petscsnes.h>

using namespace math_pack;

class CBElement;
class CBElementSolid;
class CBModel;
class CBElementFactory;
class CBSolverPlugin;
class CBSolverPluginFactory;

using Vector3f       = Vector3<TFloat>;
using Vector3fVector = std::vector<Vector3f>;
using Matrix3f       = Matrix3<TFloat>;
using Matrix3fVector = std::vector<Matrix3f>;

class CBSolver {
public:
    CBSolver()
    {}
    
    virtual ~CBSolver() {
        DeInit();
    }
    
    virtual void        Init(ParameterMap *parameters, CBModel *model);
    virtual void        DeInit();
    virtual void        Run();
    virtual CBStatus    PrepareSimulation();
    virtual std::string GetType() = 0;
    
    virtual void SetZeroVelocityAndAcceleration() {}
    
    virtual void SetZeroDisplacement() {}
    
    virtual void SetVelocity(Vector3fVector vel) {}
    
    virtual void SetAcceleration(Vector3fVector acc) {}
    
    void LoadVelocityAndAcceleration();
    void LoadActiveStress();
    
    CBStatus GetStatus() {return status_; }
    
    TInt GetNumberOfNodes() {return numNodes_; }
    
    TInt GetNumberOfElements() {return numElements_; }
    
    TInt GetNumberOfSolidElements() {return solidElements_.size(); }
    
    void GetNodeCoordinates(Vec &nodes) {VecDuplicate(nodes_, &nodes); VecCopy(nodes_, nodes); }
    
    void SetNodeCoordinates(Vec &nodes) {VecCopy(nodes,  nodes_); UpdateGhostNodesAndLinkToAdapter(); }
    
    void SetRefNodeCoordinates(Vec &nodes) {VecCopy(nodes, refNodes_); UpdateGhostNodesAndLinkToAdapter(); }
    
    CBTiming & GetTiming() {return timing_; }
    
    CBModel *GetModel() {return model_; }
    
    bool IsLoadedModel() {return LoadedModel_ != 0; }
    
    CBModel       *GetLoadedModel() {return LoadedModel_; }
    
    TFloat         GetVolumeOfSolidElements();
    Matrix3fVector GetDeformationTensors();
    
    CBData *GetActiveStressDataSource() {return activeStressData_; }
    
    void SetActiveStressDataSource(CBData *activeStressData) {activeStressData_ = activeStressData; InitActiveStress(); }
    
    PetscInt GetLocalNodesFrom() {return localNodesFrom_; }
    
    PetscInt GetLocalNodesTo() {return localNodesTo_; }
    
    TInt GetNumberOfLocalNodes() {return numLocalNodes_; }
    
    TInt GetNumberOfLocalElements() {return numLocalElements_; }
    
    TInt GetNumberOfGhostNodes() {return numGhostNodes_; }
    
    bool *GetNodesComponentsBoundaryConditionsGlobal() {return globalNodesComponentsBoundaryConditions_; }
    
    void SetNodesComponentsBoundaryConditionsGlobal(bool *globalNodesComponentsBoundaryConditions);
    void GetDomainDecomposition(ParameterMap *parameters, CBModel *model, std::vector<std::pair<int,
                                int>> &mapping);
    virtual CBStatus CalcNodalForcesJacobian(Vec displacement, Mat jacobian) = 0;
    virtual CBStatus CalcNodalForces();
    virtual CBStatus CalcNodalForcesJacobian();
    
    virtual TFloat GetKineticEnergy() {return 0; }
    
    virtual TFloat GetDeformationEnergy();
    
    virtual TFloat GetDampingEnergyDissipation() {return 0; }
    
    CBElementAdapter *GetAdapter() {return adapter_; }
    
    void UpdateActiveStress(PetscScalar time);
    void UpdateExportActiveStress(PetscScalar time, Vec &stressVec);
    void ExportNodes();
    void ExportElementsVectorData(std::string key, Vec);
    void ExportElementsScalarData(std::string key, Vec);
    void ExportNodesVectorData(std::string key, Vec);
    void ExportNodesScalarData(std::string key, Vec);
    void PrepareElements(std::vector<CBElementSolid *>);
    void PrepareElements();
    void RelaxElementsAndBases();
    void RelaxElementsAndBases(std::vector<CBElementSolid *>);
    void RelaxElementsAndBasesByMaterialIndex(TInt mat);
    void RelaxElementsAndBasesByIndex(TInt index);
    
    CBMaterial                           *GetMaterial(int i) {return materials_[i]; }
    
    const std::vector<CBElement *>      & GetElementVector();
    const std::vector<CBElementSolid *> & GetSolidElementVector();
    const std::vector<CBSolverPlugin *> & GetPlugInVector();
    
    virtual void ExportSNESMatrix(TFloat time) {}
    
    Vector3fVector GetCurrentNodeCoordinates();
    void           SetCurrentNodeCoordinates(Vector3fVector newNodes);
    
protected:
    virtual void InitFormulation();
    virtual void InitParameters();
    virtual void InitElements();
    virtual void InitNodes(CBModel *fromModel = nullptr, bool isRefNode = true);
    virtual void InitNodesComponentsBoundaryConditions();
    virtual void InitNodalForces();
    virtual void InitNodalForcesJacobian();
    virtual void InitActiveStress();
    virtual void InitTensionModels();
    virtual void InitExporter();
    virtual void InitPlugins();
    virtual void InitLoadedState();
    void         UpdateGhostNodesAndLinkToAdapter();
    void         CreateNodesJacobianAndLinkToAdapter();
    
    void ActivateAllElements() {
        for (auto &it : solidElements_)
            it->Activate();
    }
    
    void         InitNodesIndicesMapping();
    void         InitNodesIndicesMappingNonGhosted();
    void         LoadMesh();
    void         DetermineNodesRanges();
    virtual void InitVectors()     = 0;
    virtual void InitMatrices()    = 0;
    virtual void InitPETScSolver() = 0;
    void         DeInitExporter();
    virtual void Export(TFloat timeStep);
    void         ExportFiber();
    void         ExportSheet();
    void         ExportSheetNormal();
    void         ExportBasesAtAllQuadraturePoints();
    void         ExportDeformation();
    void         ExportJacobian();
    void         ExportLambda();
    void         ExportCauchy();
    void         ExportPK2Stress();
    void         ExportGreenLagrangeStrain();
    bool         UpdateExportGlobalData(std::string identifier, TFloat value);
    
    /// Method to export the local activation times to view them later (i.e. in paraview)
    void ExportLocalActivationTime();
    
    /// Helper-method. Used to debug the calculated lat values, writes the values in a txt-file
    static void         writeLATToFile(const std::vector<int> &indices, const std::vector<double> &values);
    virtual CBStatus    SolverStep(PetscScalar time, bool forceJacobianAndDampingRecalculation = false) = 0;
    virtual PetscScalar CalcFiniteDifferencesEpsilon(Vec a);
    void                ExtractSolidElements();
    void                ExtractSurfaceElements();
    void                InitActiveStressData();
    void                InitMaterials();
    void                ExportBasisAtQuadraturePoint();
    
    // Exporter Options
    
    CBModel *model_                 = 0;     // No Ownership
    CBModel *LoadedModel_           = 0;
    ParameterMap *parameters_       = 0;     // No Ownership
    ParameterMap *pluginparameters_ = 0;     // No Ownership
    
    CBTiming timing_;
    
    CBFileManager fileManager_;       /// manages multiple tension file objects
    CBTensionFactory tensionFactory_;  /// provides/creates tensionModel objects for each material
    
    CBFormulation *formulation_ = 0;
    std::map<TInt, CBMaterial *> materials_;
    
    std::vector<CBElement *> elements_; // Factory takes control over elements -> delete is not needed !!
    std::vector<CBElementSolid *> solidElements_;
    std::vector<CBElementSurface *> surfaceElements_;
    
    CBMaterialFactory materialFactory_;
    CBSolverPluginFactory *pluginFactory_ = 0;
    
    TInt numElements_ = 0;
    TInt numNodes_    = 0;
    
    TInt jacobianComponentCriteria_ = 0;
    PetscInt lagJacobian_           = -2; // unused
    PetscInt lagPreconditioner_     = -2; // unused
    
    CBStatus status_                     = CBStatus::NOT_INITIALIZED;
    SNESConvergedReason convergedReason_ = SNES_CONVERGED_FNORM_ABS;
    int snesIterations_                  = 0;
    int kspIterations_                   = 0;
    int stepUpIterations_                = -1;
    int refSnesIts_                      = 10;
    int maxSnesIts_                      = 25;
    int maxFunEval_                      = 1000;
    
    bool isInitElementsDone_            = false;
    bool isInitNodesDone_               = false;
    bool isInitNodalForcesJacobianDone_ = false;
    bool isInitFormulationDone_         = false;
    bool isExtractSolidElementsDone_    = false;
    bool isInitActiveStressDataDone_    = false;
    bool ignoreMaxIt_                   = false;
    bool exportSNESMatrix_              = false;
    
    /// Indicates whether or not the debugging of the internal lat-values was already done
    bool executedAlready = false;
    
    // ------- PETSc vectors and matrics -----------
    CBElementAdapter *adapter_ = 0;
    
    PetscInt numNonZeros_ = 0;
    std::vector<PetscInt> nodesRanges_;
    
    //! Global indices, of course ;-)
    PetscInt localNodesFrom_ = 0;
    PetscInt localNodesTo_   = 0;
    
    //! Global indices, of course ;-)
    
    PetscInt *ghostNodes_    = 0;
    PetscInt  numGhostNodes_ = 0;
    
    PetscInt numLocalElements_                            = 0;
    PetscInt numLocalNodes_                               = 0;
    PetscInt *nonZeroElementsDiagonal_                    = 0;
    PetscInt *nonZeroElementsOffDiagonal_                 = 0;
    ISLocalToGlobalMapping nodesIndicesMapping_           = 0;
    ISLocalToGlobalMapping nodesIndicesMappingNonGhosted_ = 0;
    PetscInt *elementsIndicesMapping_                     = 0;
    
    Vec nodes_               = 0;
    Vec refNodes_            = 0;
    Vec nodesSeq_            = 0;
    Vec refNodesSeq_         = 0;
    Vec nodalForces_         = 0;
    Mat nodalForcesJacobian_ = 0;
    
    PetscInt activeStressLowerIndex_                               = 0;
    PetscInt activeStressUpperIndex_                               = 0;
    PetscScalar *activeStressTensorComponents_                     = 0;
    PetscInt *activeStressTensorComponentsIndices_                 = 0;
    PetscInt *activeStressTensorIndices_                           = 0;
    Vec activeStress_                                              = 0;  // unscaled active stress
    Vec exportActiveStress_                                        = 0;  // scaled active stress for export
    Vec corruptElements_                                           = 0;
    bool *globalNodesComponentsBoundaryConditions_                 = 0;
    Mat   boundaryConditionsNodalForcesJacobianDiagonalComponents_ = 0;
    
    bool enableTiming_ = true;
    
    PetscScalar lastInterruptionTime_ = 0;
    
    // ---------- Parallel Data Exporter -----------
    //    VecScatter and Sequential vectors are used to export data to the model holder ( only on process zero !!!)
    
    Vec elementsScalarData_           = 0;
    Vec elementsScalarDataSeq_        = 0;
    VecScatter elementsScalarDataCtx_ = 0;
    
    Vec elementsVectorData_           = 0;
    Vec elementsVectorDataSeq_        = 0;
    VecScatter elementsVectorDataCtx_ = 0;
    
    Vec nodesVectorData_           = 0;
    Vec nodesVectorDataSeq_        = 0;
    VecScatter nodesVectorDataCtx_ = 0;
    
    Vec nodesScalarData_           = 0;
    Vec nodesScalarDataSeq_        = 0;
    VecScatter nodesScalarDataCtx_ = 0;
    
    PetscInt *exportElementsMapping_ = 0;
    
    // std::string snesStepsFilename_;  // unused
    
    // ---------------------------------------------
    
    CBData *activeStressData_ = 0;
    TFloat  precision_        = 0;
    TFloat  epsilon_          = 0;
    std::vector<CBSolverPlugin *> plugins_;
    
    int ExportCounter_ = 0;
    
    Vec absDisplacement_ = 0;
    Vec velocity_     = 0;
    
private:
    bool shallExit_ = false;
};  // class CBSolver
#endif  // ifndef CB_SOLVER
