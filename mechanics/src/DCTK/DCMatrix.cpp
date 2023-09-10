//
//  DCMatrix.cpp
//  CardioMechanics_unstable
//
//  Created by Thomas Fritz on 10.12.13.
//
//

#include "DCMatrix.h"
#include "DCMatrixPETSc.h"
#include "DCCtrl.h"

DCMatrix::DCMatrix()
{
    m_ = DCCtrl::NewMatrix();
}

DCMatrix::~DCMatrix()
{
    delete m_;
}

DCMatrix::DCMatrix(const DCMatrix &m)
{
    m.m_->DuplicateTo(m_);
    m.m_->CopyValuesTo(m_);
}

void DCMatrix::operator = (const DCMatrix &m)
{
    m.m_->DuplicateTo(m_);
    m.m_->CopyValuesTo(m_);
}

DCMatrix::DCMatrix(DCMatrix && m)
{
    m_ = m.m_;
    m.m_ = 0;
}

void DCMatrix::operator = (DCMatrix && m)
{
    m_ = m.m_;
    m.m_ = 0;
}

Mat DCMatrix::Petsc()
{
    DCMatrixPETSc* p = dynamic_cast<DCMatrixPETSc*>(m_);
    if(p)
        return p->GetMat();
    else
        throw std::runtime_error("Mat DCMatrix::Petsc(): DCCtrl has not been initialized with Petsc");
}
