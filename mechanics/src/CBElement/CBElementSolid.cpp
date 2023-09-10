//
//  CBElementSolid.cpp
//  CardioMechanics_unstable
//
//  Created by Thomas Fritz on 13.05.13.
//
//

extern "C" double dgesvd_(const char*,const char*,int*,int*,double*,int*,double*,double*,int*,double*,int*,double*,int*,int*);
// #endif

#include "CBElementSolid.h"

CBElementSolid::CBElementSolid(CBElementSolid& other) : CBElement(other){
    initialVolume_ = other.initialVolume_;
    isDefect_ = other.isDefect_;
    targetStrain = other.targetStrain;
}

void CBElementSolid::RepairDeformationTensorIfInverted(Matrix3<TFloat>& deformationTensor)
{
    if(deformationTensor.Det() < 0)
    {
        int m = 3, n = 3, lda = 3, ldu = 3, ldvt = 3, info, lwork;
        double wkopt;
        double* work;
        /* Local arrays */
        double s[3], u[9], vt[9];
        
        double a[9] =
        {   deformationTensor(0,0),deformationTensor(1,0),deformationTensor(2,0),
            deformationTensor(0,1),deformationTensor(1,1),deformationTensor(2,1),
            deformationTensor(0,2),deformationTensor(1,2),deformationTensor(2,2)
        };
        
        lwork = -1;
        
        dgesvd_( "All", "All", &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, &wkopt, &lwork,
                &info );
        lwork = (int)wkopt;
        work = (double*)malloc( lwork*sizeof(double) );
        /* Compute SVD */
        dgesvd_( "All", "All", &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, work, &lwork,
                &info );
        /* Check for convergence */
        if( info > 0 ) {
            printf( "The algorithm computing SVD failed to converge.\n" );
            exit( 1 );
        }
        
        
        int mi=0;
        
        for(int i=1; i< 3; i++)
            if(s[i] < s[mi])
                mi = i;
        
        s[mi] *= -1;
        
        
        Matrix3<TFloat> em( s[0],   0,  0,
                           0, s[1],  0,
                           0,    0,s[2]);
        Matrix3<TFloat> um(u);
        Matrix3<TFloat> vm(vt);
        
        um.Transpose();
        vm.Transpose();
        
        free( (void*)work );
        
        deformationTensor = um * em * vm;
    }
}

bool CBElementSolid::CheckFiberLength() {
    for(int i=0; i < this->GetNumberOfQuadraturePoints(); i++) {
        for (int j=0; j<3; j++) {
            TFloat fibernorm = this->GetBasisAtQuadraturePoint(i)->GetCol(j).Norm();
            if (fibernorm < 0.99 || fibernorm > 1.01) {
                return true;
            }
        }
    }
    return false;
}
