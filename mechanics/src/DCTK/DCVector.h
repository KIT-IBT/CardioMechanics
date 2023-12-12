/*
 * File: DCVector.h
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


#ifndef DC_VECTOR
#define DC_VECTOR

#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>
#include <stdexcept>

#include "DCType.h"

enum DCVectorType
{
    VEC_SEQUENTIAL,
    VEC_PARALLEL
};

class DCCtrl;


class DCVectorBase // Pure class
{
public:
    DCVectorBase(){}
    virtual ~DCVectorBase(){}
    
    DCVectorBase(const DCVectorBase* v);
    
    virtual void SetLocalSize(TInt size) = 0;
    virtual void SetGlobalSize(TInt size) = 0;
    
    virtual TInt GetLocalSize() = 0;
    virtual TInt GetGlobalSize() = 0;
    
    virtual void SetGhostIndices(TInt numGhostIndices, const TInt* ghostIndices) = 0;
    
    virtual void SetType(DCVectorType type)               = 0;
    virtual DCVectorType GetType()                                = 0;
    
    virtual void Build()                                  = 0;
    virtual void Assemble()                               = 0;
    virtual void DuplicateTo(DCVectorBase* v) = 0;
    virtual void CopyValuesTo(DCVectorBase* v) = 0;
    
    virtual void InsertValues(TInt numValues, const TInt* indices, const TFloat* values)=0;
    virtual void InsertValuesLocal(TInt numValues, const TInt* indices, const TFloat* values)=0;
    
    virtual void AddValues(TInt numValues, const TInt* indices, const TFloat* values)=0;
    virtual void AddValuesLocal(TInt numValues, const TInt* indices, const TFloat* values)=0;
    
    virtual void GetValues(TInt numValues, const TInt* indices, TFloat* values)=0;
    virtual void GetValuesLocal(TInt numValues, const TInt* indices, TFloat* values)=0;
    
    virtual void SetLocalToGlobalMapping(TInt* mapping)=0;
    virtual void ConvertLocalToGlobalIndices(TInt numIndices, TInt* indices)=0;
    
    //------------------------ Algorithm -------------------------
    
    
    
    //------------------------------------------------------------
    
protected:
private:
};



class DCVector // Kind of smart pointer to DCVectorBase
{
public:
    DCVector();
    DCVector(const TInt globalSize);
    DCVector(const TInt localSize, const TInt globalSize);
    virtual ~DCVector();
    DCVector(const DCVector& v);
    void operator=(const DCVector& v);
    
    
    
    // ------------------------------------------  Delegation ---------------------------------------
    
    void SetType(DCVectorType type){v_->SetType(type);}
    void SetLocalSize(TInt size){v_->SetLocalSize(size);}
    void SetGlobalSize(TInt size){v_->SetGlobalSize(size);}
    void SetGhostIndices(TInt numGhostIndices,const TInt* ghostIndices){v_->SetGhostIndices(numGhostIndices,ghostIndices);}
    
    TInt GetLocalSize(){return v_->GetLocalSize();}
    TInt GetGlobalSize(){return v_->GetGlobalSize();}
    
    void DuplicateTo(DCVector v){v_->DuplicateTo(v.v_);}
    void DuplicateFrom(DCVector v){v.v_->DuplicateTo(v_);}
    
    void CopyValuesTo(DCVector v){v_->CopyValuesTo(v.v_);}
    void CopyValuesFrom(DCVector v){v.v_->CopyValuesTo(v_);}
    
    
    void Assemble(){v_->Assemble();}
    void Build(){v_->Build();}
    
    
    virtual void InsertValues(TInt numValues, const TInt* indices, const TFloat* values){v_->InsertValues(numValues, indices, values);}
    virtual void InsertValuesLocal(TInt numValues, const TInt* indices, const TFloat* values){v_->InsertValuesLocal(numValues, indices, values);}
    
    virtual void AddValues(TInt numValues, const TInt* indices, const TFloat* values){v_->AddValues(numValues, indices, values);}
    virtual void AddValuesLocal(TInt numValues, const TInt* indices, const TFloat* values){v_->AddValuesLocal(numValues, indices, values);}
    
    virtual void GetValues(TInt numValues, const TInt* indices, TFloat* values){v_->GetValues(numValues, indices, values);}
    virtual void GetValuesLocal(TInt numValues, const TInt* indices, TFloat* values){v_->GetValuesLocal(numValues, indices, values);}
    
    void SetLocalToGlobalMapping(TInt* mapping){v_->SetLocalToGlobalMapping(mapping);}
    void ConvertLocalToGlobalIndices(TInt numIndices, TInt* indices){v_->ConvertLocalToGlobalIndices(numIndices,indices);}
    
    Vec Petsc();
    // ----------------------------------------------------------------------------------------------
    
protected:
private:
    DCVectorBase*  v_ = 0;
};





#endif
