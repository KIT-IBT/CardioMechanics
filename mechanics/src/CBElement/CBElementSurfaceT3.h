/*
 *  CBElementSurfaceT3.h
 *  CardioMechanics
 *
 *  Created by Thomas Fritz on 31.03.11.
 *  Copyright 2011 Institute for Biomedical Engineering, Karlsruhe Institute of Technology (KIT). All rights reserved.
 *
 */

#ifndef CB_ELEMENT_SURFACE_T3_H
#define CB_ELEMENT_SURFACE_T3_H

#include "CBElementCavity.h"
#include "Vector3.h"
#include "Triangle.h"
#include <array>

using namespace math_pack;


class CBElementSurfaceT3 : public CBElementCavity {
public:
    CBElementSurfaceT3() : CBElementCavity() {}
    
    CBElementSurfaceT3(CBElementSurfaceT3 &other);
    virtual ~CBElementSurfaceT3() {}
    
    static CBElement *New();
    CBElement *Clone() override;
    void SetNodeIndex(unsigned int i, TInt j) override;
    TInt GetNodeIndex(unsigned int i) override;
    
    TInt GetNumberOfQuadraturePoints() override {return 1;}
    
    std::string GetType() override {return std::string("T3"); }
    
    virtual void SetWeight(unsigned int i, TFloat w) override;
    virtual TFloat GetWeight(unsigned int i) override;
    
    unsigned int GetNumberOfNodesIndices() override {return (unsigned int)3; }
    
    Triangle<TFloat> GetTriangle() override;
    TFloat GetArea();
    Vector3<TFloat> GetCentroid();
    Vector3<TFloat> GetNormalVector();
    TFloat CalcContributionToVolume(const TFloat *referenceCoords) override;
    void CalcContributionToVolumeJacobian(TFloat *volumeJacobianEntries, TInt *volumeJacobianEntriesIndices) override;
    void GetAreaNormalAndIndices(TFloat *c, TInt *nodesCoordsIndices) override;
    void ApplyPressure(TFloat pressure);
    void CalcPressureJacobian(TFloat pressure);
    void SetBasisAtQuadraturePoint(int i, const Matrix3<TFloat> &basis) override;
    Matrix3<TFloat> *GetBasisAtQuadraturePoint(int i) override;
    
protected:
    std::array<TInt, 3> nodesIndices_;
    std::array<TFloat, 3> weight_ = {{1.0/3.0, 1.0/3.0, 1.0/3.0}};
    Matrix3<TFloat> basisAtQuadraturePoint_;
    void CalcForcesDueToPressure(TFloat pressure, const TInt *nodesCoordsIndices, const TFloat *nodesCoords,
                                 TFloat *forces);
    
private:
    CBElementSurfaceT3(const CBElementSurfaceT3 &);
    void operator=(const CBElementSurfaceT3 &);
    typedef CBElement        Base;
};

#endif // ifndef CB_ELEMENT_SURFACE_T3_H
