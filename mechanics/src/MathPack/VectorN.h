/*
 * File: VectorN.h
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

#ifndef MATHPACK_VECTORN_H
#define MATHPACK_VECTORN_H

#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>
#include <initializer_list>
#include <cassert>
#include "MatrixN.h"

namespace math_pack
{
template<typename T> class MatrixN;

template<typename T>
class VectorN
{
    
public:
    friend class MatrixN<T>;
    
    VectorN()
    {
        N_ = 0;
    }
    
    VectorN(ushort N)
    {
        N_ = N;
        e_.resize(N_);
        std::fill(e_.begin(), e_.end(), 0);
    }
    
    VectorN(ushort N, T val)
    {
        N_ = N;
        e_.resize(N_);
        std::fill(e_.begin(), e_.end(), val);
    }
    
    VectorN(ushort N, const T* const array)
    {
        N_ = N;
        e_.resize(N_);
        std::copy(array, array+N_, e_.begin());
    }
    
    VectorN(ushort N, std::initializer_list<T> l)
    {
#ifndef NDEBUG
        if(l.size() != N) throw std::runtime_error("VectorN<T>::VectorN(unsigned short N, std::initializer_list l) initializer list has to exactly contain N elements ");
#endif
        N_ = N;
        e_.resize(N_);
        std::copy(l.begin(), l.end(), e_.begin());
    }
    
    inline void Resize(ushort N)
    {
        N_ = N;
        e_.resize(N_);
    }
    
    inline ushort GetDimension() const
    {
        return(N_);
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
    
    inline T& operator()(unsigned short i)
    {
#ifndef NDEBUG
        if(i >= N_) throw std::runtime_error("VectorN<T>::operator()(unsigned short i) i out of range");
#endif
        return((T&)e_[i]);
    }
    
    inline const T& operator()(unsigned short i) const
    {
#ifndef NDEBUG
        if(i >= N_) throw std::runtime_error("VectorN<T>::operator()(unsigned short i) i out of range");
#endif
        return((T&)e_[i]);
    }
    
    inline T& Set(unsigned short i)
    {
#ifndef NDEBUG
        if(i >= N_) throw std::runtime_error("VectorN<T>::Set(unsigned short i) i out of range");
#endif
        return((T&)e_[i]);
    }
    
    inline T Get(unsigned short i) const
    {
#ifndef NDEBUG
        if(i >= N_) throw std::runtime_error("VectorN<T>::Get(unsigned short i) i out of range");
#endif
        return(e_[i]);
    }
    
    inline void SetArray(const T* const array)
    {
        std::copy(array, array+N_, e_.begin());
    }
    
    inline void SetStdVector(const std::vector<T>* const vector)
    {
        e_ = *vector;
        N_ = e_.size();
    }
    
    
    // ------------------- Operator ----------------------
    
    inline bool operator==(const VectorN<T>& a) const;
    
    inline T operator*(const VectorN<T>& a) const;
    inline VectorN<T> operator+(const T& a) const;
    inline VectorN<T> operator-(const T& a) const;
    inline VectorN<T> operator+(const VectorN<T>& a) const;
    inline VectorN<T> operator-(const VectorN<T>& a) const;
    
    inline VectorN<T> operator*(const T& b) const;
    inline VectorN<T> operator/(const T& b) const;
    
    inline VectorN<T> operator*(const MatrixN<T> &b);
    
    inline void operator*=(const T& b);
    inline void operator/=(const T& b);
    inline void operator+=(const T& b);
    inline void operator-=(const T& b);
    inline void operator+=(const VectorN<T>& a);
    inline void operator-=(const VectorN<T>& a);
    
    // ----------------------------------------------------
    
    friend std::ostream& operator<<(std::ostream& str, const VectorN<T>& a)
    {
        for(ushort i = 0; i < a.N_; i++)
            str << a.e_[i] << "\n";
        return(str);
    }
    
    void Print() const;
    inline T Norm() const;
    inline bool Normalize();
    
    std::ostream& Print(std::ostream& outputStream) const;
    
    template<class C> VectorN<C> GetConvertedVector()
    {
        VectorN<C> vec(N_);
        for(ushort i = 0; i < N_; i++)
            vec(i) = static_cast<C>(e_[i]);
        return(vec);
    }
    
protected:
    
private:
    ushort N_;
    std::vector<T> e_;
};

template<typename T>
inline bool VectorN<T>::operator==(const VectorN<T>& b) const
{
    if(b.GetDimension()!= N_)
        throw std::runtime_error("VectorN<T>::operator==(const VectorN<T>& b) Vectors have different dimensions");
    
    for(ushort i = 0; i < N_; i++)
        if(e_[i] !=  b.e_[i])
            return(false);
    return(true);
}

template<typename T>
inline VectorN<T> VectorN<T>::operator+(const T& b) const
{
    VectorN<T> vec(N_);
    for(ushort i = 0; i < N_; i++)
        vec(i) = e_[i] + b;
    return(vec);
}

template<typename T>
inline VectorN<T> VectorN<T>::operator-(const T& b) const
{
    VectorN<T> vec(N_);
    for(ushort i = 0; i < N_; i++)
        vec(i) = e_[i] - b;
    return(vec);
}

template<typename T>
inline VectorN<T> VectorN<T>::operator+(const VectorN<T>& b) const
{
    if(b.GetDimension()!= N_)
        throw std::runtime_error("VectorN<T>::operator+(const VectorN<T>& b) Vectors have different dimensions");
    
    VectorN<T> vec(N_);
    for(ushort i = 0; i < N_; i++)
        vec(i) = e_[i] + b.e_[i];
    return(vec);
}

template<typename T>
inline VectorN<T> VectorN<T>::operator-(const VectorN<T>& b) const
{
    if(b.GetDimension()!= N_)
        throw std::runtime_error("VectorN<T>::operator-(const VectorN<T>& b) Vectors have different dimensions");
    
    VectorN<T> vec(N_);
    for(ushort i = 0; i < N_; i++)
        vec(i) = e_[i] - b.e_[i];
    return(vec);
}

template<typename T>
inline T VectorN<T>::operator*(const VectorN<T>& b) const
{
    if(b.GetDimension()!= N_)
        throw std::runtime_error("VectorN<T>::operator+(const VectorN<T>& b) Vectors have different dimensions");
    
    T productSum = 0;
    for(ushort i = 0; i < N_; i++)
        productSum += e_[i] * b.e_[i];
    return(productSum);
}

template<typename T>
inline VectorN<T> VectorN<T>::operator*(const T& b) const
{
    VectorN<T> vec(N_);
    for(ushort i = 0; i < N_; i++)
        vec(i) = e_[i] * b;
    return(vec);
}

template<typename T>
inline VectorN<T> VectorN<T>::operator*(const MatrixN<T> &b)
{
    if(b.GetDimension()!= N_)
        throw std::runtime_error("VectorN<T> VectorN<T>::operator*(const MatrixN<T> &b) Matrix has not the same dimension as vector");
    
    // calculates the product of the transposed vector with a matrix --> row vector
    VectorN<T> vec(N_, (T)0);
    for(ushort i = 0; i < N_; i++)
        for(ushort j = 0; j < N_; j++)
            vec(j) += e_[i] * b(i,j);
    return(vec);
}

template<typename T>
inline VectorN<T> operator*(const T& b, const VectorN<T>& a)
{
    return(a*b);
}

template<typename T>
inline VectorN<T> VectorN<T>::operator/(const T& b) const
{
    VectorN<T> vec(N_);
    for(ushort i = 0; i < N_; i++)
        vec(i) = e_[i] / b;
    return(vec);
}

template<typename T>
inline VectorN<T> operator/(const T& b, const VectorN<T>& a)
{
    return(a/b);
}

template<typename T>
inline void VectorN<T>::operator*=(const T& b)
{
    for(ushort i = 0; i < N_; i++)
        e_[i] *=  b;
}

template<typename T>
inline void VectorN<T>::operator/=(const T& b)
{
    for(ushort i = 0; i < N_; i++)
        e_[i] /=  b;
}

template<typename T>
inline void VectorN<T>::operator+=(const T& b)
{
    for(ushort i = 0; i < N_; i++)
        e_[i] +=  b;
}

template<typename T>
inline void VectorN<T>::operator-=(const T& b)
{
    for(ushort i = 0; i < N_; i++)
        e_[i] -=  b;
}

template<typename T>
inline void VectorN<T>::operator+=(const VectorN<T>& b)
{
    if(b.GetDimension()!= N_)
        throw std::runtime_error("VectorN<T>::operator+=(const VectorN<T>& b) Vectors have different dimensions");
    
    for(ushort i = 0; i < N_; i++)
        e_[i] +=  b.e_[i];
}

template<typename T>
inline void VectorN<T>::operator-=(const VectorN<T>& b)
{
    if(b.GetDimension()!= N_)
        throw std::runtime_error("VectorN<T>::operator+=(const VectorN<T>& b) Vectors have different dimensions");
    
    for(ushort i = 0; i < N_; i++)
        e_[i] -=  b.e_[i];
}

template<typename T>
inline T VectorN<T>::Norm() const
{
    return(sqrt((*this) * (*this)));
}

template<typename T>
inline bool VectorN<T>::Normalize()
{
    T norm = Norm();
    if(norm == T(0.0))
    {
#ifndef NDEBUG
        throw std::runtime_error("VectorN<T><T>::Normalize() Division by Zero");
#endif
        return(false);
    }
    
    for(ushort i = 0; i < N_; i++)
        e_[i] /= norm;
    
    return(true);
}

template<typename T>
inline std::ostream & VectorN<T>::Print(std::ostream& outputStream) const
{
    outputStream << std::scientific;
    for(ushort i = 0; i < N_; i++)
        outputStream << e_[i] << " ";
    return(outputStream);
}

template<typename T>
void VectorN<T>::Print() const
{
    std::cout << *this;
}
}
#endif
