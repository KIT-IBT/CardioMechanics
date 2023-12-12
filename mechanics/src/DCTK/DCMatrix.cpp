/*
 * File: DCMatrix.cpp
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
