/*
 *  CBElementSolidT10.h
 *  CardioMechanics
 *
 *  Created by Thomas Fritz on 22.04.10.
 *  Copyright 2010 Institute for Biomedical Engineering, Karlsruhe Institute of Technology (KIT). All rights reserved.
 *
 */

#ifndef CB_ELEMENT_SOLID_T10_H
#define CB_ELEMENT_SOLID_T10_H

#include "CBElementSolid.h"
#include <array>

class CBElementSolidT10 : public CBElementSolid
{
public:
    CBElementSolidT10() : CBElementSolid() {}
    CBElementSolidT10(CBElementSolidT10& other);
    virtual ~CBElementSolidT10(){}
    static CBElement* New();
    CBElement* Clone();
    
    virtual std::string GetType(){return(std::string("T10"));}
    TInt GetDimension(){return(3);}
    void SetNodeIndex(unsigned int i, TInt j);
    TInt GetNodeIndex(unsigned int i);
    void CheckNodeSorting();
    virtual void UpdateShapeFunctions();
    virtual bool IsElementInverted();
    
    virtual unsigned int GetNumberOfNodesIndices(){return(10);}
    TInt GetNumberOfQuadraturePoints(){return(5);}
    CBStatus CalcNodalForcesJacobian();
    virtual CBStatus CalcNodalForces();
    TFloat GetDeformationEnergy();
    Matrix3<TFloat> GetPK2Stress();
    CBStatus SetNodalForcesToZeroIfElementIsDefect();
    CBStatus SetNodalForcesJacobianToZeroIfElementIsDefect();
    virtual CBStatus CalcNodalForcesWithoutActiveStress();
    virtual CBStatus CalcNodalForcesAndJacobian();
    CBStatus CalcConsistentMassMatrix();
    CBStatus CalcLumpedMassMatrix();
    CBStatus GetDeformationTensor(Matrix3<TFloat>& f);
    CBStatus GetDeformationTensorAtQuadraturePoints(Matrix3<TFloat>* f);
    CBStatus GetCauchyStress(Matrix3<TFloat>& cauchyStress);
    virtual Matrix3<TFloat> GetDeformationTensorWithGlobalBasis(TFloat l1, TFloat l2, TFloat l3, TFloat l4);
    virtual Matrix3<TFloat> GetRightCauchyDeformationTensorWithGlobalBasis(TFloat l1, TFloat l2, TFloat l3, TFloat l4);
    TFloat GetVolume();
    void SetBasisAtQuadraturePoint(int i, const Matrix3<TFloat>& basis);
    virtual CBStatus CalculateLaplacian();
    virtual CBStatus CalculateLaplacianT4();
    Matrix3<TFloat>* GetBasisAtQuadraturePoint(int i);
    TFloat* GetShapeFunctionsDerivatives();
    TFloat* GetT4ShapeFunctionsDerivatives();
    virtual void RotateFiberBy(TFloat phi, TFloat theta);
protected:
    virtual void CalcShapeFunctionDerivatives(TFloat l1, TFloat l2, TFloat l3, TFloat l4, TFloat* dNdX, bool useReferenceNodes = 0);
    void CalcT4ShapeFunctionsDerivatives();
    void CalcShapeFunctionDerivativesAtQuadraturePoints();
    virtual void CalcShapeFunctionDerivativesAtCentroid();
    virtual void CalcDeformationTensorsAtQuadraturePointsWithLocalBasis(const TFloat* nodesCoords, Matrix3<TFloat>* deformationTensors);
    void CalcDeformationTensorsAtQuadraturePointsWithGlobalBasis(const TFloat* nodesCoords, Matrix3<TFloat>* deformationTensors);
    virtual void CalcDeformationTensorsAtCentroidWithLocalBasis(const TFloat* nodesCoords, Matrix3<TFloat>& deformationTensors);
    virtual void CalcDeformationTensorsAtCentroidWithLocalBasisWithT4ShapeFunctions(const TFloat* nodesCoords, Matrix3<TFloat>& deformationTensors);
    void GetNodesCoordsIndices(TInt* nodesCoordsIndices);
    virtual CBStatus CalcNodalForcesHelperFunction(const TFloat* nodesCoords, const bool* boundaryConditions, TFloat* forces);
    
    std::array<TInt, 10> nodesIndices_;
    std::array<TFloat, 150> dNdXW_;  // Derivatives of the shape functions at the 4 quadrature points + center
    std::array<TFloat, 12> dNdXt4_;
    TFloat detJ_ = 0;
    
private:
    typedef CBElement        Base;
    typedef CBElementSolid   Ancestor;
    
    std::array<Matrix3<TFloat>, 5> basisAtQuadraturePoint_;
    
    CBElementSolidT10(const CBElementSolidT10 &) = delete;
    void operator = (const CBElementSolidT10&) = delete;
};

#endif
