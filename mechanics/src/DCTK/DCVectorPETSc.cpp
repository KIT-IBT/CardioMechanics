/*
 * File: DCVectorPETSc.cpp
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


#include <iostream>

#include <petscksp.h>

#include "DCCtrl.h"
#include "DCVectorPETSc.h"
#include "DCVector.h"

void DCVectorPETSc::Assemble()
{
    if (localForm_) 
    {
        VecDestroy(&localForm_);
        localForm_ = 0;
    }
    
    VecAssemblyBegin(vec_);
    VecAssemblyEnd(vec_);
    
}


void DCVectorPETSc::SetLocalSize(PetscInt size)
{
    localSize_ = size;
}

void DCVectorPETSc::SetGlobalSize(PetscInt size)
{
    globalSize_ = size;
}

PetscInt DCVectorPETSc::GetLocalSize()
{
    return localSize_;
}

PetscInt DCVectorPETSc::GetGlobalSize()
{
    return globalSize_;
}

void DCVectorPETSc::SetGhostIndices(PetscInt numGhostIndices, const PetscInt* ghostIndices)
{
    numGhostIndices_ = numGhostIndices;
    ghostIndices_ = ghostIndices;
    isGhost_ = true; 
}

void DCVectorPETSc::SetType(DCVectorType type)
{
    type_ = type;
}

DCVectorType DCVectorPETSc::GetType()
{
    return type_;
}


DCVectorPETSc::~DCVectorPETSc()
{
    if(vec_)
        VecDestroy(&vec_);
}

void DCVectorPETSc::DuplicateTo(DCVectorBase* v)
{
    
    if(v)
        delete v;
    
    v = new DCVectorPETSc;
    *reinterpret_cast<DCVectorPETSc*>(v) = *this;
    VecDuplicate(vec_, &(reinterpret_cast<DCVectorPETSc*>(v))->vec_);
}

void DCVectorPETSc::CopyValuesTo(DCVectorBase* v)
{
    VecCopy(vec_, reinterpret_cast<DCVectorPETSc*>(v)->vec_);
}

void DCVectorPETSc::Build()
{
    if(DCCtrl::GetNumberOfProcesses() == 1)
    {
        if(globalSize_ != 0)
            VecCreateSeq(PETSC_COMM_SELF,globalSize_,&vec_);
        else
            VecCreateSeq(PETSC_COMM_SELF,localSize_,&vec_);
    }
    else
    {   
        if(numGhostIndices_ == 0)
            VecCreateMPI(PETSC_COMM_WORLD,localSize_,globalSize_,&vec_);
        else
            VecCreateGhost(PETSC_COMM_WORLD,localSize_,globalSize_,numGhostIndices_,ghostIndices_,&vec_);
    }
}

void DCVectorPETSc::InsertValues(PetscInt numValues, const PetscInt* indices, const PetscScalar* values)
{
    VecSetValues(vec_, numValues, indices, values, INSERT_VALUES);
}

void DCVectorPETSc::InsertValuesLocal(PetscInt numValues, const PetscInt* indices, const PetscScalar* values)
{
    VecSetValuesLocal(vec_, numValues, indices, values, INSERT_VALUES);    
}

void DCVectorPETSc::AddValues(PetscInt numValues, const PetscInt* indices, const PetscScalar* values)
{
    VecSetValues(vec_, numValues, indices, values, ADD_VALUES);
}

void DCVectorPETSc::AddValuesLocal(PetscInt numValues, const PetscInt* indices, const PetscScalar* values)
{
    VecSetValuesLocal(vec_, numValues, indices, values, ADD_VALUES);    
}

void DCVectorPETSc::GetValues(PetscInt numValues, const PetscInt* indices, PetscScalar* values)
{
    if(isGhost_)
        throw std::runtime_error("void DCVectorPETSc::GetValues(PetscInt numValues, const PetscInt* indices, PetscScalar* values): Vector is ghost vector, use GetValuesLocal()");
    else
        VecGetValues(vec_, numValues, indices, values);
}

void DCVectorPETSc::GetValuesLocal(PetscInt numValues, const PetscInt* indices, PetscScalar* values)
{
    if(isGhost_)
    {
        if(!localForm_)
        {
            VecGhostUpdateBegin(vec_, INSERT_VALUES, SCATTER_FORWARD);
            VecGhostUpdateEnd(vec_, INSERT_VALUES, SCATTER_FORWARD);
            VecGhostGetLocalForm(vec_, &localForm_);
        }
        VecGetValues(localForm_, numValues, indices, values);
    }
    throw std::runtime_error("void DCVectorPETSc::GetValuesLocal(PetscInt numValues, const PetscInt* indices, PetscScalar* values): Use only for ghost vectors, If you have local indices, use void DCVectorPETSc::ConvertToGlobalIndices(PetscInt numIndices, TInt* indices)");
}

void DCVectorPETSc::ConvertLocalToGlobalIndices(PetscInt numIndices, PetscInt* indices)
{
    ISLocalToGlobalMapping mapping;
    VecGetLocalToGlobalMapping(vec_, &mapping);
    ISLocalToGlobalMappingApply(mapping,numIndices,indices,indices);
}

void DCVectorPETSc::SetLocalToGlobalMapping(PetscInt* mapping)
{
    if(isGhost_)
        throw std::runtime_error("void DCVectorPETSc::void DCVectorPETSc::SetLocalToGlobalMapping(PetscInt* mapping): Vector is ghost vector, mapping is generated automatically");
    else
    {
        ISLocalToGlobalMapping localToGlobalMapping;
        ISLocalToGlobalMappingCreate(PETSC_COMM_WORLD,1,localSize_,mapping,PETSC_COPY_VALUES,&localToGlobalMapping);
        VecSetLocalToGlobalMapping(vec_, localToGlobalMapping);
    }
}



