/*
 *  CBElement.h
 *  CardioMechanics
 *
 *  Created by Thomas Fritz on 23.02.10.
 *  Copyright 2010 Institute for Biomedical Engineering, Karlsruhe Institute of Technology (KIT). All rights reserved.
 *
 */

#ifndef CB_ELEMENT_H
#define CB_ELEMENT_H

#include <type_traits>

#include "Matrix3.h"

#include "CBMaterial.h"
#include "CBStatus.h"
#include "DCType.h"
class CBTensionModel;
class CBElementFactory;

using namespace math_pack;

enum class ElementProperty {
    Undefined = 0,
    Solid   = 1,
    Surface = 2,
    Active  = 4
};

inline ElementProperty operator | (ElementProperty lhs, ElementProperty rhs) {
    return static_cast<ElementProperty>(static_cast<std::underlying_type<ElementProperty>::type>(lhs) | static_cast<std::underlying_type<ElementProperty>::type>(rhs));
}

inline ElementProperty operator & (ElementProperty lhs, ElementProperty rhs) {
    return static_cast<ElementProperty>(static_cast<std::underlying_type<ElementProperty>::type>(lhs) & static_cast<std::underlying_type<ElementProperty>::type>(rhs));
}

inline bool operator && (ElementProperty lhs, ElementProperty rhs) {
    return (static_cast<std::underlying_type<ElementProperty>::type>(lhs) && static_cast<std::underlying_type<ElementProperty>::type>(rhs));
}

inline ElementProperty operator ~(ElementProperty rhs) {
    return static_cast<ElementProperty>(~static_cast<std::underlying_type<ElementProperty>::type>(rhs));
}

inline ElementProperty& operator |= (ElementProperty& lhs, ElementProperty rhs) {
    lhs = lhs | rhs;
    return lhs;
}

inline ElementProperty& operator &= (ElementProperty& lhs, ElementProperty rhs) {
    lhs = lhs & rhs;
    return lhs;
}

class CBElementAdapter;

class CBElement
{
public:
    CBElement() :
    materialIndex_(0),
    properties(ElementProperty::Active),
    adapter_(0),
    material_(0) {}
    
    CBElement(CBElement& other);
    
    virtual ~CBElement(){}
    
    // ------- Pure virtual functions -------
    
    virtual CBElement* Clone() = 0;
    
    virtual std::string GetType() = 0;
    
    virtual void SetNodeIndex(unsigned int i, TInt j) = 0;
    virtual TInt GetNodeIndex(unsigned int i)         = 0;
    virtual unsigned int GetNumberOfNodesIndices() = 0;
    
    virtual TInt GetNumberOfQuadraturePoints()     = 0;
    virtual void SetBasisAtQuadraturePoint(int i, const Matrix3<TFloat>& basis) = 0;
    virtual Matrix3<TFloat>* GetBasisAtQuadraturePoint(int i) = 0;
    
    // --------------------------------------
    
    virtual bool CheckFiberLength();
    
    bool IsActive(){return IsProperty(ElementProperty::Active);}
    void Activate(){SetProperty(ElementProperty::Active);}
    void Deactivate(){DeleteProperty(ElementProperty::Active);}
    
    void PrintNodesIndices();
    
    void SetMaterialIndex(TInt materialIndex){materialIndex_ = materialIndex;}
    TInt GetMaterialIndex() const {return(materialIndex_);}
    void SetMaterial(CBMaterial* material){material_ = material;}
    CBMaterial* GetMaterial(){return(material_);}
    
    void SetTensionModel(CBTensionModel* tension) {tensionModel_ = tension;}
    CBTensionModel* GetTensionModel(){ return(tensionModel_); }
    
    void SetIndex(TInt i){index_ = i;}
    TInt GetIndex() const {return(index_);}
    void SetLocalIndex(TInt i){localIndex_ = i;}
    TInt GetLocalIndex() const {return(localIndex_);}
    
    void SetBasis(const Matrix3<TFloat>& basis){SetBasisAtQuadraturePoint(0, basis);}
    Matrix3<TFloat>* GetBasis(){return(GetBasisAtQuadraturePoint(0));}
    
    void SetParameters(ParameterMap* parameters){parameters_=parameters;}
    
    void SetAdapter(CBElementAdapter* adapter){adapter_ = adapter;}
    CBElementAdapter* GetAdapter(){return(adapter_);}
    
    void SetProperty(ElementProperty property);
    bool IsProperty(ElementProperty property);
    
protected:
    
    void DeleteProperty(ElementProperty property);
    
    TInt index_;
    TInt localIndex_;                           // Used for domain decomposition and multiprocessing
    TInt materialIndex_;
    
    /// Flags that describe the properties of the element, like it is a solid, a surface or
    /// another type of object. See ElemetProperty enum class for more information.
    ElementProperty properties;
    
    std::vector<TInt> neighbors_;
    CBElementAdapter* adapter_;
    ParameterMap* parameters_ = 0;
    CBMaterial*      material_;
    CBTensionModel* tensionModel_ = 0;
    
private:
    //CBElement(const CBElement &) = delete;
    void operator = (const CBElement&) = delete;
};
#endif
