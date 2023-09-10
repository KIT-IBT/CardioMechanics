//
//  DCMatrixPETSc.cpp
//  CardioMechanics
//
//  Created by Thomas Fritz on 30.10.11.
//  Copyright (c) 2011 Institute for Biomedical Engineering, Karlsruhe Institute of Technology (KIT). All rights reserved.
//

#include <iostream>
#include <petscksp.h>

#include "DCMatrixPETSc.h"
#include "DCVectorPETSc.h"
#include "DCCtrl.h"

void DCMatrixPETSc::SetNumberOfLocalRows(PetscInt numLocalRows)
{
    numLocalRows_ = numLocalRows;
}

void DCMatrixPETSc::SetNumberOfGlobalRows(PetscInt numGlobalRows)
{
    numGlobalRows_ = numGlobalRows;
}

void DCMatrixPETSc::SetNumberRows(PetscInt numLocalRows, PetscInt numGlobalRows)
{
    numLocalRows_ = numLocalRows;
    numGlobalRows_ = numGlobalRows;
}

void DCMatrixPETSc::SetNumberOfLocalCols(PetscInt numLocalCols)
{
    numLocalCols_ = numLocalCols;
}

void DCMatrixPETSc::SetNumberOfGlobalCols(PetscInt numGlobalCols)
{
    numGlobalCols_ = numGlobalCols;
}

void DCMatrixPETSc::SetNumberOfCols(PetscInt numLocalCols, PetscInt numGlobalCols)
{
    numLocalCols_ = numLocalCols;
    numGlobalCols_ = numGlobalCols;
}

void DCMatrixPETSc::SetNumberOfNonzeros(PetscInt numNzDiag, PetscInt numNzOffDiag)
{
    numNzDiag_ = numNzDiag;
    numNzOffDiag_ = numNzOffDiag;
}

void DCMatrixPETSc::SetNonzeros(PetscInt* nonZerosDiag, PetscInt* nonZerosOffDiag)
{
    nonZerosDiag_ = nonZerosDiag;
    nonZerosOffDiag_ = nonZerosOffDiag;
}


void DCMatrixPETSc::Assemble()
{
    MatAssemblyBegin(mat_,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(mat_,MAT_FINAL_ASSEMBLY);
}



void DCMatrixPETSc::SetType(DCMatrixType type)
{
    type_ = type;
}

DCMatrixType DCMatrixPETSc::GetType()
{
    return type_;
}


DCMatrixPETSc::~DCMatrixPETSc()
{
    if(mat_)
        MatDestroy(&mat_);
}

void DCMatrixPETSc::DuplicateTo(DCMatrixBase* v)
{
    if(v)
        delete v;
    
    v = new DCMatrixPETSc;
    *reinterpret_cast<DCMatrixPETSc*>(v) = *this;
    
    MatDuplicate(mat_, MAT_DO_NOT_COPY_VALUES, &(reinterpret_cast<DCMatrixPETSc*>(v))->mat_);
}

void DCMatrixPETSc::DuplicateFrom(DCMatrixBase* v)
{
    if(v)
        delete v;
    v = new DCMatrixPETSc;
    *reinterpret_cast<DCMatrixPETSc*>(v) = *this;
    
    MatDuplicate(mat_, MAT_DO_NOT_COPY_VALUES, &(reinterpret_cast<DCMatrixPETSc*>(v))->mat_);
}


void DCMatrixPETSc::CopyValuesTo(DCMatrixBase* v)
{
    MatCopy(mat_, reinterpret_cast<DCMatrixPETSc*>(v)->mat_,DIFFERENT_NONZERO_PATTERN);
}


void DCMatrixPETSc::CopyValuesFrom(DCMatrixBase* v)
{
    MatCopy(mat_, reinterpret_cast<DCMatrixPETSc*>(v)->mat_,DIFFERENT_NONZERO_PATTERN);
}


void DCMatrixPETSc::Build()
{
    if(!hasBeenBuilt_)
    {
        if(DCCtrl::GetNumberOfProcesses() == 1)
        {
            if(numGlobalRows_ != 0)
                MatCreateSeqAIJ(PETSC_COMM_WORLD,numGlobalRows_,numGlobalCols_,numNzDiag_,nonZerosDiag_,&mat_);
            else
                MatCreateSeqAIJ(PETSC_COMM_WORLD,numLocalRows_,numLocalCols_,numNzDiag_,nonZerosDiag_,&mat_);
        }
        else
        {   
            MatCreateAIJ(PETSC_COMM_WORLD,numLocalRows_,numLocalCols_,numGlobalRows_,numGlobalCols_,numNzDiag_,nonZerosDiag_,numNzOffDiag_,nonZerosOffDiag_,&mat_);
        }
        
        hasBeenBuilt_ = true;
    }
    else
        throw std::runtime_error("void DCMatrixPETSc::Build(): Matrix can't be built twice");
}

