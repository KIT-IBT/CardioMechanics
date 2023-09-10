/*
 *  DCMatrix.h
 *
 *  Created by Thomas Fritz on 29.10.11.
 *  Copyright 2011 Institute for Biomedical Engineering, Karlsruhe Institute of Technology (KIT). All rights reserved.
 *
 */

#ifndef DC_MATRIX
#define DC_MATRIX

#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>
#include <stdexcept>

#include "DCType.h"
enum DCMatrixType
{
    MAT_SEQUENTIAL,
    MAT_SEQUENTIAL_DENSE,
    MAT_SEQUENTIAL_BAIJ,
    MAT_PARALLEL,
    MAT_PARALLEL_DENSE,
    MAT_PARALLEL_BAIJ
};

class DCCtrl;
class DCMatrixPETSc;

class DCMatrixBase // Pure class
{
public:
    DCMatrixBase(){}
    
    virtual ~DCMatrixBase(){}
    
    DCMatrixBase(const DCMatrixBase* v);
    DCMatrixBase(const DCMatrixBase& v);
    
    virtual void SetNumberOfLocalRows(TInt numLocalRows) = 0;
    virtual void SetNumberOfGlobalRows(TInt numLocalRows)             = 0;
    virtual void SetNumberRows(TInt numLocalRows, TInt numGlobalRows) = 0;
    
    virtual void SetNumberOfLocalCols(TInt numLocalCols) = 0;
    virtual void SetNumberOfGlobalCols(TInt numGlobalCols) = 0;
    virtual void SetNumberOfCols(TInt numLocalCols, TInt numGlobalCols) = 0;
    
    virtual void SetNumberOfNonzeros(TInt numNzDiag, TInt numNzOffDiag) = 0;
    virtual void SetNonzeros(TInt* nonZerosDiag, TInt* nonZerosDiagOff) = 0;
    
    virtual void SetType(DCMatrixType type) = 0;
    virtual DCMatrixType GetType()          = 0;
    
    virtual void Build() = 0;
    virtual void Assemble() = 0;
    virtual void DuplicateTo(DCMatrixBase* m)  = 0;
    virtual void CopyValuesTo(DCMatrixBase* m) = 0;
    
    /*  TInt GetLocalSize(){return localSize_;}
     * TInt GetGlobalSize(){}
     *
     *
     * virtual void SetValue(TInt index, TFloat value){}
     * virtual void SetValues(TInt numValues, const TInt* indices, const TFloat* values){}
     * virtual void SetValuesLocal(TInt numValues, const TInt* indices, const TFloat* values){}
     *
     * virtual TFloat GetValue(TInt index){}
     * virtual void GetValues(TInt numValues, const TInt* indices, TFloat* values){}
     * virtual void GetValuesLocal(TInt numValues, const TInt* indices, TFloat* values){}
     *
     * virtual void Zero(){}
     * virtual void AXPY(TInt a, DCMatrix* y){}
     * virtual void AYPX(TInt a, DCMatrix* y){} */
    
protected:
private:
};


class DCMatrix // Kind of smart pointer to DCMatrixBase
{
public:
    DCMatrix();
    virtual ~DCMatrix();
    DCMatrix(const DCMatrix &m);
    void operator = (const DCMatrix &m);
    DCMatrix(DCMatrix && m);
    void operator = (DCMatrix && m);
    
    int a;
    int b;
    
    // ------------------------------------------  Delegation ---------------------------------------
    
    Mat Petsc();
    
    // ------------------------------------------  Delegation ---------------------------------------
    virtual void SetNumberOfLocalRows(TInt numLocalRows){m_->SetNumberOfLocalRows(numLocalRows);}
    virtual void SetNumberOfGlobalRows(TInt numGlobalRows){m_->SetNumberOfGlobalRows(numGlobalRows);}
    virtual void SetNumberRows(TInt numLocalRows, TInt numGlobalRows){m_->SetNumberRows(numLocalRows, numGlobalRows);}
    
    virtual void SetNumberOfLocalCols(TInt numLocalCols){m_->SetNumberOfLocalCols(numLocalCols);}
    virtual void SetNumberOfGlobalCols(TInt numGlobalCols){m_->SetNumberOfGlobalCols(numGlobalCols);}
    virtual void SetNumberOfCols(TInt numLocalCols, TInt numGlobalCols){m_->SetNumberOfCols(numLocalCols, numGlobalCols);}
    
    virtual void SetNumberOfNonzeros(TInt numNzDiag, TInt numNzOffDiag){m_->SetNumberOfNonzeros(numNzDiag, numNzOffDiag);}
    virtual void SetNonzeros(TInt* nonZerosDiag, TInt* nonZerosOffDiag){m_->SetNonzeros(nonZerosDiag, nonZerosOffDiag);}
    
    void SetType(DCMatrixType type){m_->SetType(type);}
    
    void Assemble(){m_->Assemble();}
    void Build(){m_->Build();}
    
    void DuplicateTo(DCMatrix m){m_->DuplicateTo(m.m_);}
    void DuplicateFrom(DCMatrix m){m.m_->DuplicateTo(m_);}
    
    void CopyValuesTo(DCMatrix m){m_->CopyValuesTo(m.m_);}
    void CopyValuesFrom(DCMatrix m){m.m_->CopyValuesTo(m_);}
    
    // ----------------------------------------------------------------------------------------------
    
protected:
private:
    DCMatrixBase* m_ = 0;
};
#endif
