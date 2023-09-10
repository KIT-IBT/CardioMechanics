//
//  Quaternion.h
//
//  Created by Lukas Baron on 23.06.2017.
//  Copyright 2009 Institute for Biomedical Engineering, Karlsruhe Institute of Technology (KIT). All rights reserved.
//

#pragma once

#include "Matrix3.h"

namespace math_pack
{

std::vector<TFloat> RotationMatrixToQuaternion(Matrix3<TFloat> m);
Matrix3<TFloat> QuaternionToRotationMatrix(std::vector<TFloat> q);
TFloat QuaternionScp(std::vector<TFloat> u, std::vector<TFloat> v);


/// converts a rotation matrix to a quaternion (given as vector)
/// source: http://en.wikipedia.org/wiki/Rotation_matrix#Quaternion
std::vector<TFloat> RotationMatrixToQuaternion(Matrix3<TFloat> m)
{
    TFloat Qxx = m(0,0);
    TFloat Qyy = m(1,1);
    TFloat Qzz = m(2,2);
    TFloat Qzy = m(2,1);
    TFloat Qxz = m(0,2);
    TFloat Qyx = m(1,0);
    TFloat Qyz = m(1,2);
    TFloat Qzx = m(2,0);
    TFloat Qxy = m(0,1);
    
    TFloat t = Qxx+Qyy+Qzz; // (trace of Q)
    TFloat r = sqrt(1+t);
    TFloat w = 0.5*r;
    TFloat x = std::copysign(0.5*sqrt(1+Qxx-Qyy-Qzz), Qzy-Qyz);
    TFloat y = std::copysign(0.5*sqrt(1-Qxx+Qyy-Qzz), Qxz-Qzx);
    TFloat z = std::copysign(0.5*sqrt(1-Qxx-Qyy+Qzz), Qyx-Qxy);
    
    std::vector<TFloat> v;
    
    v.push_back(w);
    v.push_back(x);
    v.push_back(y);
    v.push_back(z);
    return(v);
}

/// converts a quaternion (given as vector) to a rotation matrix
/// source: http://en.wikipedia.org/wiki/Rotation_matrix#Quaternion
Matrix3<TFloat> QuaternionToRotationMatrix(std::vector<TFloat> q)
{
    Matrix3<TFloat> m;
    TFloat          w = q[0];
    TFloat          x = q[1];
    TFloat          y = q[2];
    TFloat          z = q[3];
    TFloat          n = w * w + x * x + y * y + z * z;
    TFloat          s = 0;
    
    if(n != 0)
        s = 2./n;
    
    TFloat wx = s * w * x;
    TFloat wy = s * w * y;
    TFloat wz = s * w * z;
    TFloat xx = s * x * x;
    TFloat xy = s * x * y;
    TFloat xz = s * x * z;
    TFloat yy = s * y * y;
    TFloat yz = s * y * z;
    TFloat zz = s * z * z;
    
    m(0,0) = 1 - (yy + zz);
    m(1,0) =      xy + wz;
    m(2,0) =      xz - wy;
    
    m(0,1) =      xy - wz;
    m(1,1) = 1 - (xx + zz);
    m(2,1) =      yz + wx;
    
    m(0,2) =       xz + wy;
    m(1,2) =       yz - wx;
    m(2,2) =  1 - (xx + yy);
    return(m);
}



/// Scalarproduct between two vectors of arbitrary size
TFloat QuaternionScp(std::vector<TFloat> u, std::vector<TFloat> v)
{
    TFloat sum = 0;
    if((u.size() != 4) || (v.size() != 4))
        throw std::runtime_error("CBGenerateFiberOrientation::QuaternionScp() : One or both quaternions do not have the size 4!");
    if(u.size() == v.size())
        for(int i = 0; i < u.size(); i++)
            sum += u[i]*v[i];
    else
        throw std::runtime_error("CBGenerateFiberOrientation::QuaternionScp() : Vectors do not have the same size!");
    
    return(sum);
}

} // namespace math_pack
