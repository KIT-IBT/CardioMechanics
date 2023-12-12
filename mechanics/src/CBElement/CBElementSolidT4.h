/*
 * File: CBElementSolidT4.h
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


#ifndef CB_ELEMENT_SOLID_T4_H
#define CB_ELEMENT_SOLID_T4_H

#include "CBElementSolid.h"
#include <array>

class CBElementSolidT4 : public CBElementSolid {
public:
    CBElementSolidT4() : CBElementSolid() {}
    
    CBElementSolidT4(CBElementSolidT4 &other);
    
    virtual ~CBElementSolidT4() {}
    
    static CBElement *New();
    CBElement *Clone();
    void SetNodeIndex(unsigned int i, TInt j);
    TInt GetNodeIndex(unsigned int i);
    
    unsigned int GetNumberOfNodesIndices() {return (unsigned int)4;}
    
    TInt GetNumberOfQuadraturePoints() {return 1;}
    
    std::string GetType() {return std::string("T4");}
    
    void CheckNodeSorting();
    CBStatus CalcNodalForcesJacobian();
    CBStatus CalcNodalForces();
    CBStatus CalcNodalForcesAndJacobian();
    CBStatus CalcConsistentMassMatrix(); // {std::runtime_error("Error: Function CBElementSolidT4::CalcConsistenMassMatrix() is not implemented yet."); return CBStatus::FAILED;}
    CBStatus CalcLumpedMassMatrix();
    CBStatus CalculateLaplacian();
    virtual bool IsElementInverted();
    
    CBStatus CalcDampingMatrix() {
        std::runtime_error("Error: Function CBElementSolidT4::CalcDampingMatrix() is not implemented yet.");
        return CBStatus::FAILED;
    }
    
    TFloat GetDeformationEnergy();
    Matrix3<TFloat> GetPK2Stress();
    CBStatus GetCauchyStress(Matrix3<TFloat> &cauchyStress);
    virtual CBStatus GetDeformationTensor(Matrix3<TFloat> &f);
    
    void UpdateShapeFunctions() {CalcShapeFunctionsDerivatives();}
    
    TFloat GetVolume();
    void SetBasisAtQuadraturePoint(int i, const Matrix3<TFloat> &basis);
    Matrix3<TFloat> *GetBasisAtQuadraturePoint(int i);
    static std::vector<double> *conds;
    TFloat *GetShapeFunctionsDerivatives();
    
protected:
    void CalcShapeFunctionsDerivatives();
    void CalcShapeFunctionDerivatives(TFloat l1, TFloat l2, TFloat l3, TFloat l4, TFloat *dNdX,
                                      bool useReferenceNodes = false);
    void CalcDeformationTensorWithLocalBasis(const TFloat *nodesCoords, Matrix3<TFloat> &deformationTensor);
    void GetNodesCoordsIndices(TInt *nodesCoordsIndices);
    virtual CBStatus CalcNodalForcesHelperFunction(const TFloat *nodesCoords, const bool *boundaryConditions, TFloat *forces);
    std::array<TInt, 4>   nodesIndices_;
    std::array<TFloat, 12> dNdX_;     // partial derivates of the shape function Ni: âNi/âx Since this element is linear, the derivates are constant over the whole element.
    TFloat detJ_ = 0;
    
    // ek717: needed for the CBAccelerate Plug in
    CBStatus GetDeformationTensorAtQuadraturePoints(Matrix3<TFloat> *f);
    
private:
    typedef CBElement        Base;
    typedef CBElementSolid   Ancestor;
    Matrix3<TFloat> basisAtQuadraturePoint_;
    CBElementSolidT4(const CBElementSolidT4 &);
    void operator=(const CBElementSolidT4 &);
};
#endif // ifndef CB_ELEMENT_SOLID_T4_H
