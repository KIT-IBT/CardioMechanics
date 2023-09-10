/*
 *  MatrixN.h
 *  This class is only made for small (max. N = 256), quadratic (NxN) Matrices!
 *
 *  Created by Thomas Fritz on 24.08.09, extended by Steffen Schuler
 *  Copyright 2009 Institute for Biomedical Engineering, Karlsruhe Institute of Technology (KIT). All rights reserved.
 *
 */
#ifndef MATHPACK_MATRIXN_H
#define MATHPACK_MATRIXN_H

#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>

#include "VectorN.h"

namespace math_pack
{
template<typename T> class VectorN;

template<typename T>
class MatrixN
{
public:
    friend class VectorN<T>;
    
    MatrixN()
    {
        N_ = 0;
        NN_ = 0;
        e_.resize(NN_);
    }
    
    MatrixN(ushort N)
    {
        N_ = N;
        NN_ = N_*N_;
        e_.resize(NN_);
        std::fill(e_.begin(), e_.end(), 0);
    }
    
    MatrixN(ushort N, T val)
    {
        N_ = N;
        NN_ = N_*N_;
        e_.resize(NN_);
        std::fill(e_.begin(), e_.end(), val);
    }
    
    MatrixN(ushort N, std::initializer_list<T> l)
    {
#ifndef NDEBUG
        if(l.size() != N*N) throw std::runtime_error("MatrixN<T>::MatrixN(ushort N, std::initializer_list l) initializer list has to exactly contain N*N elements ");
#endif
        N_ = N;
        NN_ = N_*N_;
        e_.resize(NN_);
        std::copy(l.begin(), l.end(), e_.begin());
    }
    
    MatrixN(const MatrixN<T> &m) = default;
    
    MatrixN(ushort N, const T* array)
    {
        N_ = N;
        NN_ = N_*N_;
        e_.resize(NN_);
        std::copy(array, array+NN_, e_.begin());
    }
    
    inline void Resize(ushort N)
    {
        if(N < N_)
        {
            e_.erase(e_.end()-(N_-N)*N_, e_.end());
            for(ushort i = N; i > 0; i--)
                e_.erase(e_.begin()+(i-1)*N_+N, e_.begin()+i*N_);
        }
        else if(N > N_)
        {
            for(ushort i = N_; i > 0; i--)
                e_.insert(e_.begin()+i*N_, N-N_, 0);
            e_.insert(e_.end(), (N-N_)*N, 0);
        }
        else
            return;
        
        N_ = N;
        NN_ = N_*N_;
    }
    
    inline ushort GetDimension() const
    {
        return(N_);
    }
    
    inline T& operator()(ushort i, ushort j)
    {
#ifndef NDEBUG
        if(i >= N_ || j >= N_)
            throw std::runtime_error("MatrixN<T>::operator()(ushort i,ushort j) i or j out of range");
#endif
        return((T&)e_[N_*i+j]);
    }
    
    inline const T& operator()(ushort i, ushort j) const
    {
#ifndef NDEBUG
        if(i >= N_ || j >= N_)
            throw std::runtime_error("MatrixN<T>::operator()(ushort i,ushort j) i or j out of range");
#endif
        return((T&)e_[N_*i+j]);
    }
    
    inline T& operator()(uint i)
    {
#ifndef NDEBUG
        if(i >= NN_)
            throw std::runtime_error("MatrixN<T>::operator()(ushort i) i out of range");
#endif
        return((T&)e_[i]);
    }
    
    inline const T& operator()(uint i) const
    {
#ifndef NDEBUG
        if(i >= NN_)
            throw std::runtime_error("MatrixN<T>::operator()(ushort i) i out of range");
#endif
        return((T&)e_[i]);
    }
    
    inline T& Set(ushort i, ushort j)
    {
#ifndef NDEBUG
        if(i >= N_ || j >= N_)
            throw std::runtime_error("MatrixN<T>::Set(ushort i,ushort j) i or j out of range");
#endif
        return((T&)e_[N_*i+j]);
    }
    
    inline T& Set(uint i)
    {
#ifndef NDEBUG
        if(i >= NN_)
            throw std::runtime_error("MatrixN<T>::Set(ushort i) i out of range");
#endif
        return((T&)e_[i]);
    }
    
    inline T Get(uint i) const
    {
#ifndef NDEBUG
        if(i >= NN_)
            throw std::runtime_error("MatrixN<T>::Get(ushort i) i out of range");
#endif
        return(e_[i]);
    }
    
    inline T Get(ushort i, ushort j) const
    {
#ifndef NDEBUG
        if(i >= N_ || j >= N_)
            throw std::runtime_error("MatrixN<T>::Get(ushort i,ushort j) i or j out of range");
#endif
        return(e_[N_*i+j]);
    }
    
    inline T* GetArray()
    {
        return(e_.data());
    }
    
    inline const T* GetArray() const
    {
        return(e_.data());
    }
    
    inline std::vector<T> GetStdVector()
    {
        return(e_);
    }
    
    inline void SetArray(const T* const array)
    {
        std::copy(array, array+NN_, e_.begin());
    }
    
    inline void SetStdVector(const std::vector<T>* const vector)
    {
        ushort N = sqrt(vector->size());
        if(N*N != vector->size())
            throw std::runtime_error("SetStdVector(const std::vector<T>* const vector) Size of vector is not of the form n^2");
        
        N_ = N;
        NN_ = N_*N_;
        e_ = *vector;
    }
    
    
    // ------------------- Operator -------------------------------
    
    inline bool operator==(const MatrixN<T> &a) const;
    inline MatrixN<T> operator*(const MatrixN<T> &a) const;
    inline MatrixN<T> operator+(const MatrixN<T> &a) const;
    inline MatrixN<T> operator-(const MatrixN<T> &a) const;
    
    inline void operator*=(const MatrixN<T> &a);
    inline void operator+=(const MatrixN<T> &a);
    inline void operator-=(const MatrixN<T> &a);
    
    inline MatrixN<T> operator*(const T &a) const;
    inline MatrixN<T> operator/(const T &a) const;
    
    inline VectorN<T> operator*(const VectorN<T> &b) const;
    
    inline void operator*=(const T &a);
    inline void operator/=(const T &a);
    
    // -------------------------------------------------------------
    
    inline VectorN<T> GetRow(ushort i) const;
    
    inline VectorN<T> GetCol(ushort i) const;
    inline void SetRow(ushort i, VectorN<T> b);
    inline void SetCol(ushort i, VectorN<T> b);
    
    inline MatrixN<T> GetSubMatrix(ushort row, ushort col) const;
    
    inline bool Invert();
    
    inline MatrixN<T> GetInverse() const;
    inline void Transpose();
    
    inline MatrixN<T> GetTranspose() const;
    
    inline void SetToZero();
    inline void SetToIdentityMatrix();
    inline T Det() const;
    
    inline T Trace() const
    {
        T trace = 0;
        for(ushort i = 0; i < N_; i++)
            trace += e_[i*(N_+1)];
        return(trace);
    }
    
    inline T Invariant1() const { return(Trace()); }
    inline T Invariant2() const { return(0.5 * (Trace()*Trace() - ((*this) * (*this)).Trace())); }
    inline T Invariant3() const { return(Det()); }
    
    friend std::ostream& operator<<(std::ostream& str,const MatrixN<T> &a)
    {
        str << std::setprecision(5) << std::fixed;
        for(uint i = 0; i < a.N_; i++)
        {
            for(uint j = 0; j < a.N_; j++)
                str << std::setw(10) << a.e_[a.N_*i+j] << " ";
            str << std::endl;
        }
        return(str);
    }
    
    void Print() const;
    
    template<class C>
    MatrixN<C> GetConvertedMatrix() const
    {
        MatrixN<C> mat(N_);
        for(uint i = 0; i < NN_; i++)
            mat(i) = static_cast<C>(e_[i]);
        return(mat);
    }
    
protected:
    
private:
    inline T RecursiveDet(const MatrixN<T> &mat) const;
    
    ushort N_;
    uint NN_;
    std::vector<T> e_;
};

template<typename T>
inline bool MatrixN<T>::operator==(const MatrixN<T> &a) const
{
    if(a.GetDimension()!= N_)
        throw std::runtime_error("MatrixN<T>::operator==(const MatrixN<T> &a) Matrices have different dimensions");
    
    for(uint i = 0; i < NN_; i++)
        if(e_[i] != a.e_[i])
            return(false);
    
    return(true);
}

template<typename T>
inline MatrixN<T> MatrixN<T>::operator+(const MatrixN<T> &b) const
{
    if(b.GetDimension()!= N_)
        throw std::runtime_error("MatrixN<T>::operator+(const MatrixN<T> &b) Matrices have different dimensions");
    
    MatrixN<T> mat(N_);
    for(uint i = 0; i < NN_; i++)
        mat(i) = e_[i] + b.e_[i];
    return(mat);
}

template<typename T>
inline MatrixN<T> MatrixN<T>::operator-(const MatrixN<T> &b) const
{
    if(b.GetDimension()!= N_)
        throw std::runtime_error("MatrixN<T>::operator-(const MatrixN<T> &b) Matrices have different dimensions");
    
    MatrixN<T> mat(N_);
    for(uint i = 0; i < NN_; i++)
        mat(i) = e_[i] - b.e_[i];
    return(mat);
}

template<typename T>
inline MatrixN<T> MatrixN<T>::operator*(const MatrixN<T> &b) const
{
    if(b.GetDimension()!= N_)
        throw std::runtime_error("MatrixN<T>::operator*(const MatrixN<T> &b) Matrices have different dimensions");
    
    MatrixN<T> mat(N_, (T)0);
    for(ushort i = 0; i < N_; i++)
        for(ushort j = 0; j < N_; j++)
            for(ushort k = 0; k < N_; k++)
                mat(i,j) += Get(i,k) * b(k,j);
    return(mat);
}

template<typename T>
inline VectorN<T> MatrixN<T>::operator*(const VectorN<T> &b) const
{
    if(b.GetDimension()!= N_)
        throw std::runtime_error("MatrixN<T>::operator*(const VectorN<T> &b) Vector does not have the same dimension as the matrix");
    
    VectorN<T> vec(N_, (T)0);
    for(ushort i = 0; i < N_; i++)
        for(ushort j = 0; j < N_; j++)
            vec(i) += Get(i,j) * b.e_[j];
    return(vec);
}

template<typename T>
inline void MatrixN<T>::operator+=(const MatrixN<T> &b)
{
    if(b.GetDimension()!= N_)
        throw std::runtime_error("MatrixN<T>::operator+=(const MatrixN<T> &b) Matrices have different dimensions");
    
    for(uint i = 0; i < NN_; i++)
        e_[i] +=  b.e_[i];
}

template<typename T>
inline void MatrixN<T>::operator-=(const MatrixN<T> &b)
{
    if(b.GetDimension()!= N_)
        throw std::runtime_error("MatrixN<T>::operator-=(const MatrixN<T> &b) Matrices have different dimensions");
    
    for(uint i = 0; i < NN_; i++)
        e_[i] -=  b.e_[i];
}

template<typename T>
inline void MatrixN<T>::operator*=(const MatrixN<T> &b)
{
    *this = (*this) * b;
}

template<typename T>
inline MatrixN<T> MatrixN<T>::operator*(const T &b) const
{
    MatrixN<T> mat(N_);
    for(uint i = 0; i < NN_; i++)
        mat(i) = e_[i] * b;
    return(mat);
}

template<typename T>
inline MatrixN<T> MatrixN<T>::operator/(const T &b) const
{
    MatrixN<T> mat(N_);
    for(uint i = 0; i < NN_; i++)
        mat(i) = e_[i] / b;
    return(mat);
}

template<typename T>
inline MatrixN<T> operator*(const T &b, MatrixN<T> a)
{
    return(a*b);
}

template<typename T>
inline MatrixN<T> operator/(const T &b, MatrixN<T> a)
{
    return(a/b);
}

template<typename T>
inline void MatrixN<T>::operator *= (const T &b)
{
    for(uint i = 0; i < NN_; i++)
        e_[i] *= b;
}

template<typename T>
inline void MatrixN<T>::operator /= (const T &b)
{
    for(uint i = 0; i < NN_; i++)
        e_[i] /= b;
}

template<typename T>
inline VectorN<T> MatrixN<T>::GetRow(ushort i) const
{
#ifndef NDEBUG
    if(i >= N_)
        throw std::runtime_error("MatrixN::GetRow(ushort i) i out of range");
#endif
    VectorN<T> vec(N_);
    for(ushort j = 0; j < N_; j++)
        vec(j) = e_[N_*i+j];
    return(vec);
}

template<typename T>
inline VectorN<T> MatrixN<T>::GetCol(ushort j) const
{
#ifndef NDEBUG
    if(j >= N_)
        throw std::runtime_error("MatrixN<T>::GetCol(ushort j) j out of range");
#endif
    VectorN<T> vec(N_);
    for(ushort i = 0; i < N_; i++)
        vec(i) = e_[N_*i+j];
    return(vec);
}

template<typename T>
inline void MatrixN<T>::SetRow(ushort i, VectorN<T> b)
{
#ifndef NDEBUG
    if(i >= N_)
        throw std::runtime_error("MatrixN<T>::SetRow(ushort i) i out of range");
    if(b.GetDimension() != N_)
        throw std::runtime_error("MatrixN<T>::SetRow(ushort i, VectorN<T> &b) Vector has not the same dimension as matrix");
#endif
    for(ushort j = 0; j < N_; j++)
        e_[N_*i+j] = b.e_[j];
}

template<typename T>
inline void MatrixN<T>::SetCol(ushort j, VectorN<T> b)
{
#ifndef NDEBUG
    if(j >= N_)
        throw std::runtime_error("MatrixN<T>::SetCol(ushort j) j out of range");
    if(b.GetDimension() != N_)
        throw std::runtime_error("MatrixN<T>::SetCol(ushort i, VectorN<T> &b) Vector has not the same dimension as matrix");
#endif
    for(ushort i = 0; i < N_; i++)
        e_[N_*i+j] =b.e_[i];
}

template<typename T>
inline T MatrixN<T>::Det() const
{
    switch(N_)
    {
        case 1:
            return(e_[0]);
            
        case 2:
            return(e_[0]*e_[3] - e_[1]*e_[2]);
            
        case 3:
            return(e_[0]*e_[4]*e_[8] + e_[3]*e_[7]*e_[2] + e_[6]*e_[1]*e_[5] - e_[2]*e_[4]*e_[6] - e_[5]*e_[7]*e_[0] - e_[8]*e_[3]*e_[1]);
            
        case 4:
            T d[4];
            d[0] =  e_[5]*e_[10]*e_[15] - e_[5]*e_[11]*e_[14] - e_[9]*e_[6]*e_[15] + e_[9]*e_[7]*e_[14] + e_[13]*e_[6]*e_[11] - e_[13]*e_[7]*e_[10];
            d[1] = -e_[4]*e_[10]*e_[15] + e_[4]*e_[11]*e_[14] + e_[8]*e_[6]*e_[15] - e_[8]*e_[7]*e_[14] - e_[12]*e_[6]*e_[11] + e_[12]*e_[7]*e_[10];
            d[2] =  e_[4]*e_[ 9]*e_[15] - e_[4]*e_[11]*e_[13] - e_[8]*e_[5]*e_[15] + e_[8]*e_[7]*e_[13] + e_[12]*e_[5]*e_[11] - e_[12]*e_[7]*e_[ 9];
            d[3] = -e_[4]*e_[ 9]*e_[14] + e_[4]*e_[10]*e_[13] + e_[8]*e_[5]*e_[14] - e_[8]*e_[6]*e_[13] - e_[12]*e_[5]*e_[10] + e_[12]*e_[6]*e_[ 9];
            return(e_[0]*d[0] + e_[1]*d[1] + e_[2]*d[2] + e_[3]*d[3]);
            
        default:
            return RecursiveDet(*this);
    }
}

template<typename T>
inline bool MatrixN<T>::Invert()
{
    switch(N_)
    {
        case 1:
        {
            if(e_[0] == 0)
            {
                e_[0] = INFINITY;
                return(false);
            }
            e_[0] = 1 / e_[0];
            return(true);
        }
        case 2:
        {
            T det = e_[0]*e_[3] - e_[1]*e_[2];
            if(det == 0)
            {
                std::fill(e_.begin(), e_.end(), INFINITY);
                return(false);
            }
            T i[4];
            std::copy(e_.begin(), e_.end(), i);
            e_[0] =  i[3] / det;
            e_[1] = -i[1] / det;
            e_[2] = -i[2] / det;
            e_[3] =  i[0] / det;
            return(true);
        }
        case 3:
        {
            T i[9];
            i[0] = -e_[5]*e_[7] + e_[4]*e_[8];
            i[1] =  e_[2]*e_[7] - e_[1]*e_[8];
            i[2] = -e_[2]*e_[4] + e_[1]*e_[5];
            i[3] =  e_[5]*e_[6] - e_[3]*e_[8];
            i[4] = -e_[2]*e_[6] + e_[0]*e_[8];
            i[5] =  e_[2]*e_[3] - e_[0]*e_[5];
            i[6] = -e_[4]*e_[6] + e_[3]*e_[7];
            i[7] =  e_[1]*e_[6] - e_[0]*e_[7];
            i[8] = -e_[1]*e_[3] + e_[0]*e_[4];
            
            T det = e_[0]*i[0] + e_[1]*i[3] + e_[2]*i[6];
            
            if(det == 0)
            {
                std::fill(e_.begin(), e_.end(), INFINITY);
                return(false);
            }
            for(int j = 0; j < 9; j++)
                e_[j] = i[j] / det;
            return(true);
        }
        case 4:
        {
            T i[16];
            i[0]  =  e_[5]*e_[10]*e_[15] - e_[5]*e_[11]*e_[14] - e_[9]*e_[6]*e_[15] + e_[9]*e_[7]*e_[14] + e_[13]*e_[6]*e_[11] - e_[13]*e_[7]*e_[10];
            i[1]  = -e_[1]*e_[10]*e_[15] + e_[1]*e_[11]*e_[14] + e_[9]*e_[2]*e_[15] - e_[9]*e_[3]*e_[14] - e_[13]*e_[2]*e_[11] + e_[13]*e_[3]*e_[10];
            i[2]  =  e_[1]*e_[ 6]*e_[15] - e_[1]*e_[ 7]*e_[14] - e_[5]*e_[2]*e_[15] + e_[5]*e_[3]*e_[14] + e_[13]*e_[2]*e_[ 7] - e_[13]*e_[3]*e_[ 6];
            i[3]  = -e_[1]*e_[ 6]*e_[11] + e_[1]*e_[ 7]*e_[10] + e_[5]*e_[2]*e_[11] - e_[5]*e_[3]*e_[10] - e_[ 9]*e_[2]*e_[ 7] + e_[ 9]*e_[3]*e_[ 6];
            i[4]  = -e_[4]*e_[10]*e_[15] + e_[4]*e_[11]*e_[14] + e_[8]*e_[6]*e_[15] - e_[8]*e_[7]*e_[14] - e_[12]*e_[6]*e_[11] + e_[12]*e_[7]*e_[10];
            i[5]  =  e_[0]*e_[10]*e_[15] - e_[0]*e_[11]*e_[14] - e_[8]*e_[2]*e_[15] + e_[8]*e_[3]*e_[14] + e_[12]*e_[2]*e_[11] - e_[12]*e_[3]*e_[10];
            i[6]  = -e_[0]*e_[ 6]*e_[15] + e_[0]*e_[ 7]*e_[14] + e_[4]*e_[2]*e_[15] - e_[4]*e_[3]*e_[14] - e_[12]*e_[2]*e_[ 7] + e_[12]*e_[3]*e_[ 6];
            i[7]  =  e_[0]*e_[ 6]*e_[11] - e_[0]*e_[ 7]*e_[10] - e_[4]*e_[2]*e_[11] + e_[4]*e_[3]*e_[10] + e_[ 8]*e_[2]*e_[ 7] - e_[ 8]*e_[3]*e_[ 6];
            i[8]  =  e_[4]*e_[ 9]*e_[15] - e_[4]*e_[11]*e_[13] - e_[8]*e_[5]*e_[15] + e_[8]*e_[7]*e_[13] + e_[12]*e_[5]*e_[11] - e_[12]*e_[7]*e_[ 9];
            i[9]  = -e_[0]*e_[ 9]*e_[15] + e_[0]*e_[11]*e_[13] + e_[8]*e_[1]*e_[15] - e_[8]*e_[3]*e_[13] - e_[12]*e_[1]*e_[11] + e_[12]*e_[3]*e_[ 9];
            i[10] =  e_[0]*e_[ 5]*e_[15] - e_[0]*e_[ 7]*e_[13] - e_[4]*e_[1]*e_[15] + e_[4]*e_[3]*e_[13] + e_[12]*e_[1]*e_[ 7] - e_[12]*e_[3]*e_[ 5];
            i[11] = -e_[0]*e_[ 5]*e_[11] + e_[0]*e_[ 7]*e_[ 9] + e_[4]*e_[1]*e_[11] - e_[4]*e_[3]*e_[ 9] - e_[ 8]*e_[1]*e_[ 7] + e_[ 8]*e_[3]*e_[ 5];
            i[12] = -e_[4]*e_[ 9]*e_[14] + e_[4]*e_[10]*e_[13] + e_[8]*e_[5]*e_[14] - e_[8]*e_[6]*e_[13] - e_[12]*e_[5]*e_[10] + e_[12]*e_[6]*e_[ 9];
            i[13] =  e_[0]*e_[ 9]*e_[14] - e_[0]*e_[10]*e_[13] - e_[8]*e_[1]*e_[14] + e_[8]*e_[2]*e_[13] + e_[12]*e_[1]*e_[10] - e_[12]*e_[2]*e_[ 9];
            i[14] = -e_[0]*e_[ 5]*e_[14] + e_[0]*e_[ 6]*e_[13] + e_[4]*e_[1]*e_[14] - e_[4]*e_[2]*e_[13] - e_[12]*e_[1]*e_[ 6] + e_[12]*e_[2]*e_[ 5];
            i[15] =  e_[0]*e_[ 5]*e_[10] - e_[0]*e_[ 6]*e_[ 9] - e_[4]*e_[1]*e_[10] + e_[4]*e_[2]*e_[ 9] + e_[ 8]*e_[1]*e_[ 6] - e_[ 8]*e_[2]*e_[ 5];
            
            T det = e_[0]*i[0] + e_[1]*i[4] + e_[2]*i[8] + e_[3]*i[12];
            
            if(det == 0)
            {
                std::fill(e_.begin(), e_.end(), INFINITY);
                return(false);
            }
            for(int j = 0; j < 16; j++)
                e_[j] = i[j] / det;
            return(true);
        }
        default:
        {
            T det = RecursiveDet(*this);
            
            if(det == 0)
            {
                std::fill(e_.begin(), e_.end(), INFINITY);
                return(false);
            }
            MatrixN<T> mat = *this;
            for(ushort i = 0; i < N_; i++)
            {
                for(ushort j = 0; j < N_; j++)
                {
                    MatrixN<T> subMat = mat.GetSubMatrix(j,i);
                    T cofactor = ((i+j)%2 == 1 ? -1.0 : 1.0) * subMat.Det();
                    Set(i,j) = cofactor / det;
                }
            }
            return(true);
        }
    }
}

template<typename T>
MatrixN<T> MatrixN<T>::GetSubMatrix(ushort row, ushort col) const
{
    // extract submatrix needed to calculate the minor of element (row,col)
    MatrixN<T> subMat(N_-1);
    ushort rowCount = 0;
    for(ushort i = 0; i < N_; i++)
    {
        if(i != row)
        {
            ushort colCount = 0;
            for(ushort j = 0; j < N_; j++)
            {
                if(j != col)
                {
                    subMat(rowCount,colCount) = Get(i,j);
                    colCount++;
                }
            }
            rowCount++;
        }
    }
    return subMat;
}

template<typename T>
T MatrixN<T>::RecursiveDet(const MatrixN<T> &mat) const
{
    // calculate the determinant recursively using Laplace's formula
    T det = 0;
    for(ushort j = 0; j < N_; j++ )
    {
        MatrixN<T> subMat = mat.GetSubMatrix(0,j);
        T cofactor = (j%2 == 1 ? -1.0 : 1.0) * subMat.Det();
        det += mat(0,j) * cofactor;
    }
    return det;
}

template<typename T>
inline MatrixN<T> MatrixN<T>::GetInverse() const
{
    MatrixN<T> mat = *this;
    mat.Invert();
    return(mat);
}

template<typename T>
inline void MatrixN<T>::Transpose()
{
    *this = this->GetTranspose();
}

template<typename T>
inline MatrixN<T> MatrixN<T>::GetTranspose() const
{
    MatrixN<T> mat(N_);
    for(ushort i = 0; i < N_; i++)
        for(ushort j = 0; j < N_; j++)
            mat(i,j) = Get(j,i);
    return(mat);
}

template<typename T>
inline void MatrixN<T>::SetToZero()
{
    std::fill(e_.begin(), e_.end(), 0);
}

template<typename T>
inline void MatrixN<T>::SetToIdentityMatrix()
{
    this->SetToZero();
    for(ushort i = 0; i < N_; i++)
        e_[i*(N_+1)] = 1;
}

template<typename T>
void MatrixN<T>::Print() const
{
    std::cout << *this;
}

template<typename T>
inline MatrixN<T> DyadicProduct(const VectorN<T>& a, const VectorN<T>& b)
{
    ushort N = a.GetDimension();
    if(N != b.GetDimension())
        throw std::runtime_error("MatrixN<T> DyadicProduct(const VectorN<T>& a, const VectorN<T>& b) a and b have different dimensions, but I can only handle quadratic matrices");
    
    MatrixN<T> mat(N);
    for(ushort i = 0; i < N; i++)
        for(ushort j = 0; j < N; j++)
            mat(i,j) = a(i) * b(j);
    return(mat);
}

template<typename T>
MatrixN<T> Identity(ushort N_)
{
    MatrixN<T> eye(N_);
    eye.SetToIdentityMatrix();
    return(eye);
}
}

#endif
