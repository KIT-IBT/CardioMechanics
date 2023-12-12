/*
 * File: CBElementSurfaceT6.h
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


#ifndef CB_ELEMENT_SURFACE_T6_H
#define CB_ELEMENT_SURFACE_T6_H

#include "CBElementCavity.h"
#include "Vector3.h"
#include "Triangle.h"
#include <array>

using namespace math_pack;


class CBElementSurfaceT6 : public CBElementCavity
{
public:
    CBElementSurfaceT6(){}
    CBElementSurfaceT6(CBElementSurfaceT6& other);
    
    static CBElement* New();
    CBElement* Clone() override;
    
    virtual ~CBElementSurfaceT6(){}
    void SetNodeIndex(unsigned int i, TInt j) override;
    TInt GetNodeIndex(unsigned int i) override;
    TInt GetNumberOfQuadraturePoints() override {return 1;}
    std::string GetType() override {return(std::string("T6")); }
    unsigned int GetNumberOfNodesIndices() override {return((unsigned int)6); }
    
    Triangle<TFloat> GetTriangle();
    TFloat GetArea();
    Vector3<TFloat> GetCentroid();
    
    void SetBasisAtQuadraturePoint(int i, const Matrix3<TFloat>& basis) override;
    Matrix3<TFloat>* GetBasisAtQuadraturePoint(int i) override;
    
    void ApplyPressure(TFloat pressure) override;
    void CalcPressureJacobian(TFloat pressure) override;
    TFloat CalcContributionToVolume(const TFloat* referenceCoords) override;
    void CalcContributionToVolumeJacobian(TFloat* volumeJacobianEntries, TInt* volumeJacobianEntriesIndices) override;
protected:
    std::array<TInt,6>   nodesIndices_;
    Matrix3<TFloat> basisAtQuadraturePoint_;
    TFloat GetArea(const TFloat* nodesCoords);
    Matrix3<TFloat> GetNormalVectorAtQuadraturePoints(const TFloat* nodesCoords);
    void CalcForcesDueToPressure(TFloat pressure, const TInt* nodesCoordsIndices, const TFloat* nodesCoords, TFloat* forces);
private:
    
    CBElementSurfaceT6(const CBElementSurfaceT6&);
    void operator=(const CBElementSurfaceT6&);
    typedef CBElement        Base;
};

#endif
