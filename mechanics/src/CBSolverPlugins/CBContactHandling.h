/* -------------------------------------------------------
 
 CBContactHandling.cpp
 
 Ver. 1.0.0
 
 Created:       Thomas Fritz   (19.09.2011)
 Last modified: Tobias Gerach  (01.02.2021)
 
 Institute of Biomedical Engineering
 Karlsruhe Institute of Technology (KIT)
 
 http://www.ibt.kit.edu
 
 Copyright 2000-2009 - All rights reserved.
 
 ------------------------------------------------------ */

#ifndef CB_CONTACT_HANDLING_H
#define CB_CONTACT_HANDLING_H

#include <map>
#include <set>
#include "Matrix3.h"

#include "CBSolverPlugin.h"

#include "CBElementContactMaster.h"
#include "CBElementContactSlave.h"

using namespace math_pack;

class CBContactHandling : public CBSolverPlugin {
public:
    CBContactHandling();
    virtual ~CBContactHandling() {
        for (auto it : masterElements_)
            delete it;
        
        for (auto it : slaveElements_)
            delete it;
        
        delete jacobianBuffer_;
    }
    
    void Init() override;
    void Apply(PetscScalar time) override;
    void ApplyToNodalForces() override;
    void ApplyToNodalForcesJacobian() override;
    void StepBack() override;
    void Export(TFloat time) override;
    void WriteToFile(TFloat time) override;
    
    void SetAlpha(TFloat a) {alpha_ = a; }
    
    TFloat GetAlpha() {return alpha_; }
    
    void           GetMasterNodesDistancesToSlaveElements(Vec *d);
    std::set<TInt> GetMasterNodesLocalIndices();
    std::set<TInt> GetMasterWithSlaveNodesLocalIndices();
    void           Prepare() override;
    double         GetPreparationProgress() override;
    
    std::string GetName() override {return "ContactHandling"; }
    
protected:
private:
    friend class CBParameterEstimator;
    void DetermineSlaveNodes();
    void DetermineInitialSlaveElementsAtGaussPoints();
    bool CheckIfSlave(TFloat *slaveNodes, int slaveInd, Vector3<TFloat> *p, Vector3<TFloat> *nv, TFloat &dist);
    int  SearchForSlave(TFloat *slaveNodes, int oldSlave, Vector3<TFloat> *p, Vector3<TFloat> *nv, TFloat &dist);
    void DetermineSlaveElementsAtGaussPoints();
    void DetermineSlaveElementsAtVertices();
    void LoadMasterElements();
    void LoadSlaveElements();
    void UpdateDistancesMasterSlave();
    void CalcContributionToContactForceAtGaussPoint(const Triangle<TFloat> &masterTriangle,
                                                    const Triangle<TFloat> &slaveTriangle, TInt gaussPointIndex,
                                                    TFloat *nodalForces, TFloat *distances, TFloat scaling);
    Vector3<TFloat> CalculateDistanceAtGaussPoint(const Triangle<TFloat> &masterTriangle,
                                                  const Triangle<TFloat> &slaveTriangle, TInt gaussPointIndex);
    Vector3<TFloat> CalculateDistanceAtVertex(const Triangle<TFloat> &masterTriangle,
                                              const Triangle<TFloat> &slaveTriangle, TInt index);
    std::vector<CBElementContactMaster *> masterElements_;
    std::vector<CBElementContactSlave *>  slaveElements_;
    
    // -----
    std::map<int, std::vector<int> *> slaveNeighbors_;
    int maxDepth_ = 10;
    
    // -----
    
    Vec slaveElementsNodes_;
    Vec masterNodesDistToSlaveElements_;
    std::set<TInt> masterNodesLocalIndices_;
    std::vector<TInt> slaveElementsSurfaceIndices_;
    Vec slaveElementsNodesSeq_;
    
    std::vector<TInt> slaveElementsNodesIndicesGlobal_;
    
    PetscInt numGlobalSlaveElements_;
    
    VecScatter scatter_;
    PetscScalar *jacobianBuffer_;
    
    std::string filename_;
    std::string initType_;
    std::string surfaceNormalDirection_;
    TFloat dt_, prevDt_;
    PetscScalar maxDistanceToSlave_;
    PetscScalar transitionDistance_;
    PetscScalar maxAngle_;
    PetscScalar alpha_;
    PetscScalar beta_;
    PetscScalar maxAlpha_;
    PetscScalar cntMax_;
    PetscInt    cnt_;
    PetscScalar lastTime_  = 0;
    PetscScalar time_      = 0;
    TInt normalVectorSign_ = 1;
    TInt oneWayForce_;
    
    TFloat startTime_;
    bool   hasStarted_ = false;
    
    bool useContactHandlingToFitPeri_ = false;
    
    bool isFirstStep_ = true;
    bool stepBack_    = false;
    
    // ----- Values needed for export -----
    
    std::vector<Vector3<TFloat>> masterContactForces_;
    std::vector<Vector3<TFloat>> masterCorrespondingSlaveNormal_;
    std::vector<Vector3<TFloat>> masterContactDistances_;
    std::vector<TInt> masterCorrespondingSlaveFound_;
    PetscScalar averageDist_;
    PetscScalar averageContactPressure_;
    PetscScalar globalAverageDist_;  // only process zero !!!
    PetscScalar globalAverageContactPressure_;  // only process zero !!!
    bool export_;
    
    typedef CBSolverPlugin Base;
};
#endif  // ifndef CB_CONTACT_HANDLING_H
