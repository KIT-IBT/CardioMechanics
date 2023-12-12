/*
 * File: DCCtrlPETSc.cpp
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


#include "DCCtrl.h"
#include "DCCtrlPETSc.h"
#include "DCVectorPETSc.h"
#include "DCMatrixPETSc.h"
#include<sstream>

int DCCtrlPETSc::InitImpl(int argc, char** argv)
{
    PetscErrorCode ierr;
    
    ierr = PetscInitialize(&argc, &argv, (char*)0, (char*)0);
    CHKERRQ(ierr);
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &mpiRank_);
    CHKERRQ(ierr);
    ierr = MPI_Comm_size(PETSC_COMM_WORLD, &mpiSize_);
    CHKERRQ(ierr);
    return(0);
}


unsigned int DCCtrlPETSc::GetProcessIDImpl()
{
    return(mpiRank_);
}


unsigned int DCCtrlPETSc::GetNumberOfProcessesImpl()
{
    return(mpiSize_);
}

bool DCCtrlPETSc::IsProcessZeroImpl()
{
    if(mpiRank_ == 0)
        return true;
    else
        return false;
}


bool DCCtrlPETSc::IsParallelImpl()
{
    if(mpiSize_ > 1)
        return true;
    else
        return false;
}


DCVectorBase* DCCtrlPETSc::NewVectorImpl()
{
    return(new DCVectorPETSc);
}


DCMatrixBase* DCCtrlPETSc::NewMatrixImpl()
{
    return(new DCMatrixPETSc);
}


void DCCtrlPETSc::GatherToZeroImpl(std::vector<int>& a)
{
    if(IsParallel())
    {
        if(IsProcessZero())
        {
            std::vector<int> sizePerProcess;
            
            for(int i=1; i < GetNumberOfProcesses(); i++)
            {
                int n=0;
                MPI_Status stat;
                MPI_Recv(&n, 1, MPI_INT, i, i, PETSC_COMM_WORLD, &stat);
                sizePerProcess.push_back(n);
            }
            
            int sizeGlobal = 0;
            for(auto i:sizePerProcess)
                sizeGlobal+=i;
            
            
            int* recv = new int[sizeGlobal];
            int* p = recv;
            for(int i=1; i < GetNumberOfProcesses(); i++)
            {
                MPI_Status stat;
                MPI_Recv(p, sizePerProcess[i-1], MPI_INT, i, i, PETSC_COMM_WORLD, &stat);
                p += sizePerProcess[i-1];
            }
            
            a.insert(a.end(),recv,recv+sizeGlobal);
        }
        else
        {
            int n=a.size();
            
            for(int i=1; i < DCCtrl::GetNumberOfProcesses(); i++)
                if(DCCtrl::GetProcessID()==i)
                    MPI_Send(&n, 1, MPI_INT, 0, i, PETSC_COMM_WORLD);
            
            for(int i=1; i < DCCtrl::GetNumberOfProcesses(); i++)
                if(DCCtrl::GetProcessID()==i)
                    MPI_Send(&a[0],n,MPI_INT,0,i,PETSC_COMM_WORLD);
        }
    }
    
    // else do nothing;
    
}

void DCCtrlPETSc::GatherToZeroImpl(std::vector<long>& a)
{
    
    if(DCCtrl::IsParallel())
    {
        if(DCCtrl::IsProcessZero())
        {
            std::vector<long> totalSize;
            totalSize.push_back(0); // Sum of all
            totalSize.push_back(0); // from 0
            
            for(int i=1; i < DCCtrl::GetNumberOfProcesses(); i++)
            {
                int n=0;
                MPI_Status stat;
                MPI_Recv(&n, 1, MPI_INT, i, i, PETSC_COMM_WORLD, &stat);
                totalSize.push_back(n);
            }
            
            for(auto i:totalSize)
                totalSize.at(0)+=i;
            for (int i=2; i<DCCtrl::GetNumberOfProcesses(); i++)
                totalSize.at(i)+=totalSize.at(i-1);
            
            long* recv = new long[totalSize.at(0)];
            
            for(int i=1; i < DCCtrl::GetNumberOfProcesses(); i++)
            {
                MPI_Status stat;
                MPI_Recv(&recv[totalSize.at(i)], totalSize.at(i+1), MPI_LONG, i, i, PETSC_COMM_WORLD, &stat);
            }
            
            a.insert(a.end(),recv,recv+totalSize.at(0));
        }
        else
        {
            int n=a.size();
            
            for(int i=1; i < DCCtrl::GetNumberOfProcesses(); i++)
                if(DCCtrl::GetProcessID()==i)
                    MPI_Send(&n, 1, MPI_INT, 0, i, PETSC_COMM_WORLD);
            
            for(int i=1; i < DCCtrl::GetNumberOfProcesses(); i++)
                if(DCCtrl::GetProcessID()==i)
                    MPI_Send(&a[0],n,MPI_LONG,0,i,PETSC_COMM_WORLD);
        }
    }
    
    // else do nothing;
    
}


void DCCtrlPETSc::GatherToZeroImpl(std::vector<float>& a)
{
    
    if(DCCtrl::IsProcessZero())
    {
        std::vector<int> totalSize;
        totalSize.push_back(0); // Sum of all
        totalSize.push_back(0); // from 0
        
        for(int i=1; i < DCCtrl::GetNumberOfProcesses(); i++)
        {
            int n=0;
            MPI_Status stat;
            MPI_Recv(&n, 1, MPI_INT, i, i, PETSC_COMM_WORLD, &stat);
            totalSize.push_back(n);
        }
        
        for(auto i:totalSize)
            totalSize.at(0)+=i;
        for (int i=2; i<DCCtrl::GetNumberOfProcesses(); i++)
            totalSize.at(i)+=totalSize.at(i-1);
        
        float* recv = new float[totalSize.at(0)];
        
        for(int i=1; i < DCCtrl::GetNumberOfProcesses(); i++)
        {
            MPI_Status stat;
            MPI_Recv(&recv[totalSize.at(i)], totalSize.at(i+1), MPI_FLOAT, i, i, PETSC_COMM_WORLD, &stat);
        }
        
        a.insert(a.end(),recv,recv+totalSize.at(0));
    }
    else
    {
        int n=a.size();
        
        for(int i=1; i < DCCtrl::GetNumberOfProcesses(); i++)
            if(DCCtrl::GetProcessID()==i)
                MPI_Send(&n, 1, MPI_INT, 0, i, PETSC_COMM_WORLD);
        
        for(int i=1; i < GetNumberOfProcesses(); i++)
            if(GetProcessID()==i)
                MPI_Send(&a[0],n,MPI_FLOAT,0,i,PETSC_COMM_WORLD);
    }
    
}


void DCCtrlPETSc::GatherToZeroImpl(std::vector<double>& a)
{
    
    if(IsProcessZero())
    {
        std::vector<int> totalSize;
        totalSize.push_back(0); // Sum of all
        totalSize.push_back(0); // from 0
        
        for(int i=1; i < GetNumberOfProcesses(); i++)
        {
            int n=0;
            MPI_Status stat;
            MPI_Recv(&n, 1, MPI_INT, i, i, PETSC_COMM_WORLD, &stat);
            totalSize.push_back(n);
        }
        
        for(auto i:totalSize)
            totalSize.at(0)+=i;
        for (int i=2; i<DCCtrl::GetNumberOfProcesses(); i++)
            totalSize.at(i)+=totalSize.at(i-1);
        
        double* recv = new double[totalSize.at(0)];
        
        for(int i=1; i < GetNumberOfProcesses(); i++)
        {
            MPI_Status stat;
            MPI_Recv(&recv[totalSize.at(i)], totalSize.at(i+1), MPI_DOUBLE, i, i, PETSC_COMM_WORLD, &stat);
        }
        
        a.insert(a.end(),recv,recv+totalSize.at(0));
    }
    else
    {
        int n=a.size();
        
        for(int i=1; i < GetNumberOfProcesses(); i++)
            if(GetProcessID()==i)
                MPI_Send(&n, 1, MPI_INT, 0, i, PETSC_COMM_WORLD);
        
        for(int i=1; i < GetNumberOfProcesses(); i++)
            if(GetProcessID()==i)
                MPI_Send(&a[0],n,MPI_DOUBLE,0,i,PETSC_COMM_WORLD);
    }
    
}

void DCCtrlPETSc::GatherToZeroImpl(std::vector<std::string>& a)
{
    
    if(IsParallel())
    {
        
        if(IsProcessZero())
        {
            for(int i=1; i < GetNumberOfProcesses(); i++)
            {
                MPI_Status stat;
                int n = 0;
                MPI_Recv(&n, 1, MPI_INT, i, i, PETSC_COMM_WORLD, &stat);
                
                for(int j=0; j < n; j++)
                {
                    int c = 0;
                    MPI_Recv(&c,1,MPI_LONG, i, i, PETSC_COMM_WORLD, &stat);
                    char* s = new char[c];
                    MPI_Recv(s,c+1,MPI_CHAR, i, i, PETSC_COMM_WORLD, &stat);
                    a.push_back(std::string(s));
                    delete s;
                }
            }
            
        }
        else
        {
            for(int i=1; i < GetNumberOfProcesses(); i++)
                if(GetProcessID()==i)
                {
                    int n=a.size();
                    MPI_Send(&n,1,MPI_LONG, 0, i, PETSC_COMM_WORLD);
                    
                    for(int j=0; j < n; j++)
                    {
                        int c = a.at(j).length();
                        MPI_Send(&c,1,MPI_LONG, 0, i, PETSC_COMM_WORLD);
                        MPI_Send((void*)(a.at(j).c_str()),c+1,MPI_CHAR, 0, i, PETSC_COMM_WORLD);
                    }
                }
            
        }
    }
    
    // else do nothing;
}

void DCCtrlPETSc::CopyFromZeroToAllImpl(std::vector<int>& a)
{
    
    int n=0;
    
    if(IsProcessZero())
        n = a.size();
    
    MPI_Bcast(&n, 1, MPI_INT, 0, PETSC_COMM_WORLD);
    
    int* t = new int[n];
    
    if(IsProcessZero())
        memcpy(t, &a[0], n*sizeof(int));
    
    MPI_Bcast(t, n, MPI_INT, 0, PETSC_COMM_WORLD);
    a = std::vector<int>(t,t+n);
    
}

void DCCtrlPETSc::CopyFromZeroToAllImpl(std::vector<long>& a)
{
    
    int n=0;
    
    if(IsProcessZero())
        n = a.size();
    
    MPI_Bcast(&n, 1, MPI_INT, 0, PETSC_COMM_WORLD);
    
    long* t = new long[n];
    
    if(IsProcessZero())
        memcpy(t, &a[0], n*sizeof(long));
    
    MPI_Bcast(t, n, MPI_LONG, 0, PETSC_COMM_WORLD);
    a = std::vector<long>(t,t+n);
    
}

void DCCtrlPETSc::CopyFromZeroToAllImpl(std::vector<float>& a)
{
    
    int n=0;
    
    if(IsProcessZero())
        n = a.size();
    
    MPI_Bcast(&n, 1, MPI_INT, 0, PETSC_COMM_WORLD);
    
    float* t = new float[n];
    
    if(IsProcessZero())
        memcpy(t, &a[0], n*sizeof(float));
    
    MPI_Bcast(t, n, MPI_FLOAT, 0, PETSC_COMM_WORLD);
    a = std::vector<float>(t,t+n);
    
}


void DCCtrlPETSc::CopyFromZeroToAllImpl(std::vector<double>& a)
{
    
    int n=0;
    
    if(IsProcessZero())
        n = a.size();
    
    MPI_Bcast(&n, 1, MPI_INT, 0, PETSC_COMM_WORLD);
    
    double* t = new double[n];
    
    if(IsProcessZero())
        memcpy(t, &a[0], n*sizeof(double));
    
    MPI_Bcast(t, n, MPI_DOUBLE, 0, PETSC_COMM_WORLD);
    a = std::vector<double>(t,t+n);
    
}

/// copies a value from all processes to zero process and computes the average value under given weights, e.g. numLocalElements
void DCCtrlPETSc::WeightedAverageImpl(double& localAverage, double& localWeight, double& globalAverage)
{
    if(DCCtrl::IsParallel())
    {
        std::vector<double> t;
        t.push_back(localAverage);
        t.push_back(localWeight);
        Petsc::GatherToZero(t);
        globalAverage = 0;
        
        if(DCCtrl::IsProcessZero())
        {
            double globalWeightsSum = 0;
            for(int i = 0; i < t.size(); i += 2)
            {
                globalAverage      += t.at(i)*t.at(i+1);
                globalWeightsSum += t.at(i+1);
            }
            if (globalWeightsSum!=0)
                globalAverage /= globalWeightsSum;
        }
    }
    else
        globalAverage = localAverage;
}

void DCCtrlPETSc::CreateVector(PetscInt n,PetscInt N, Vec* v)
{
    if(IsParallel())
        VecCreateMPI(PETSC_COMM_WORLD, n, N, v);
    else
        VecCreateSeq(PETSC_COMM_WORLD, n, v);
}

void DCCtrlPETSc::CreateSeqVector(PetscInt n,Vec* v)
{
    VecCreateSeq(PETSC_COMM_WORLD, n, v);
}

void DCCtrlPETSc::CreateMatrix(PetscInt n, PetscInt m,PetscInt N, PetscInt M, PetscInt numNonZeros, Mat* v)
{
    if(DCCtrlPETSc::IsParallel())
    {
        MatCreateAIJ(PETSC_COMM_WORLD, n, m, N, M,numNonZeros, PETSC_NULL, numNonZeros, PETSC_NULL, v);
    }
    else
    {
        MatCreateSeqAIJ(PETSC_COMM_WORLD, N, M,numNonZeros, PETSC_NULL, v);
    }
}

void DCCtrlPETSc::CreateSeqMatrix(PetscInt N, PetscInt M, PetscInt numNonZeros, Mat* v)
{
    MatCreateSeqAIJ(PETSC_COMM_WORLD, N, M,numNonZeros, PETSC_NULL, v);
}

std::string DCCtrlPETSc::SNESReasonToString(int i)
{
    if (snesReasonStr_.count(i)==1)
        return snesReasonStr_.at(i);
    else
    {
        std::stringstream ss;
        std::string str;
        ss << i;
        ss << "_unknown";
        ss >> str;
        return str;
    }
}

std::string DCCtrlPETSc::KSPReasonToString(int i)
{
    if (kspReasonStr_.count(i)==1)
        return kspReasonStr_.at(i);
    else
    {
        std::stringstream ss;
        std::string str;
        ss << i;
        ss << "_unknown";
        ss >> str;
        return str;
    }
}

std::map<int,std::string> DCCtrlPETSc::snesReasonStr_ = {
    {2,"SNES_CONVERGED_FNORM_ABS"},
    {3,"SNES_CONVERGED_FNORM_RELATIVE"},
    {4,"SNES_CONVERGED_SNORM_RELATIVE"},
    {5,"SNES_CONVERGED_ITS"},
    {7,"SNES_CONVERGED_TR_DELTA"},
    {-1,"SNES_DIVERGED_FUNCTION_DOMAIN"},
    {-2,"SNES_DIVERGED_FUNCTION_COUNT"},
    {-3,"SNES_DIVERGED_LINEAR_SOLVE"},
    {-4,"SNES_DIVERGED_FNORM_NAN"},
    {-5,"SNES_DIVERGED_MAX_IT"},
    {-6,"SNES_DIVERGED_LINE_SEARCH"},
    {-7,"SNES_DIVERGED_INNER"},
    {-8,"SNES_DIVERGED_LOCAL_MIN"},
    {0,"SNES_CONVERGED_ITERATING"}
};

std::map<int,std::string> DCCtrlPETSc::kspReasonStr_ = {
    {1,"KSP_CONVERGED_RTOL_NORMAL"},
    {9,"KSP_CONVERGED_ATOL_NORMAL"},
    {2,"KSP_CONVERGED_RTOL"},
    {3,"KSP_CONVERGED_ATOL"},
    {4,"KSP_CONVERGED_ITS"},
    {5,"KSP_CONVERGED_CG_NEG_CURVE"},
    {6,"KSP_CONVERGED_CG_CONSTRAINED"},
    {7,"KSP_CONVERGED_STEP_LENGTH"},
    {8,"KSP_CONVERGED_HAPPY_BREAKDOWN"},
    {-2,"KSP_DIVERGED_NULL"},
    {-3,"KSP_DIVERGED_ITS"},
    {-4,"KSP_DIVERGED_DTOL"},
    {-5,"KSP_DIVERGED_BREAKDOWN"},
    {-6,"KSP_DIVERGED_BREAKDOWN_BICG"},
    {-7,"KSP_DIVERGED_NONSYMMETRIC"},
    {-8,"KSP_DIVERGED_INDEFINITE_PC"},
    {-9,"KSP_DIVERGED_NANORINF"},
    {-10,"KSP_DIVERGED_INDEFINITE_MAT"},
    {-11,"KSP_DIVERGED_PCSETUP_FAILED"},
    {0,"KSP_CONVERGED_ITERATING"}
};

void GatherMessagesOnProcessZero()
{
    
}
