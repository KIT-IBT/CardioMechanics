/*
 * File: Matrix4.h
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

#ifndef MATHPACK_MATRIX4_H
#define MATHPACK_MATRIX4_H

#include <cmath>
#include <iomanip>
#include <iostream>

#include "Vector4.h"

// forward declaration of Matrix4

namespace math_pack
{
template<typename T> class Vector4;

template<typename T>
class Matrix4

{
public:
    friend class Vector4<T>;
    
    Matrix4(T e0, T e1, T e2, T e3, T e4, T e5, T e6, T e7, T e8, T e9, T e10, T e11, T e12, T e13, T e14, T e15)
    {
        elements_[0]  = e0;
        elements_[1]  = e1;
        elements_[2]  = e2;
        elements_[3]  = e3;
        elements_[4]  = e4;
        elements_[5]  = e5;
        elements_[6]  = e6;
        elements_[7]  = e7;
        elements_[8]  = e8;
        elements_[9]  = e9;
        elements_[10] = e10;
        elements_[11] = e11;
        elements_[12] = e12;
        elements_[13] = e13;
        elements_[14] = e14;
        elements_[15] = e15;
    }
    
    Matrix4(const Vector4<T> &a, const Vector4<T> &b, const Vector4<T> &c, const Vector4<T> &d)
    {
        elements_[0]  = a(0);
        elements_[1]  = b(0);
        elements_[2]  = c(0);
        elements_[3]  = d(0);
        
        elements_[4]  = a(1);
        elements_[5]  = b(1);
        elements_[6]  = c(1);
        elements_[7]  = d(1);
        
        elements_[8]  = a(2);
        elements_[9]  = b(2);
        elements_[10] = c(2);
        elements_[11] = d(2);
        
        elements_[12] = a(3);
        elements_[13] = b(3);
        elements_[14] = c(3);
        elements_[15] = d(3);
    }
    
    Matrix4(){}
    Matrix4(const Matrix4<T> &m) = default;
    
    Matrix4(const T* array)
    {
        memcpy(elements_,array,16*sizeof(T));
    }
    
    Matrix4(std::initializer_list<T> l)
    {
#ifndef NDEBUG
        if(l.size() != 16)
            throw std::runtime_error("Matrix4<T>::Matrix4(std::initializer_list l) initializer list has to exactly contain 9 elements ");
#endif
        T* e = elements_;
        
        for(auto i:l)
        {
            *e = i;
            ++e;
        }
    }
    
    inline T& operator()(unsigned short i, unsigned short j)
    {
#ifndef NDEBUG
        if((i > 3) || (j > 3))
            throw std::runtime_error("Matrix4<T>::GetElement(unsigned short i,unsigned short j) i or j out of range");
#endif
        return((T&)elements_[4*i+j]);
    }
    
    inline const T& operator()(unsigned short i, unsigned short j) const
    {
#ifndef NDEBUG
        if((i > 3) || (j > 3))
            throw std::runtime_error("Matrix4<T>::GetElement(unsigned short i,unsigned short j) i or j out of range");
#endif
        return((T&)elements_[4*i+j]);
    }
    
    
    inline T& operator()(unsigned short i)
    {
#ifndef NDEBUG
        if((i < 0) || (i > 15))
            throw std::runtime_error("Matrix4<T>::GetElement(unsigned short i) i out of range");
#endif
        return((T&)elements_[i]);
    }
    
    inline const T& operator()(unsigned short i) const
    {
#ifndef NDEBUG
        if((i < 0) || (i > 15))
            throw std::runtime_error("Matrix4<T>::GetElement(unsigned short i) i out of range");
#endif
        return((T&)elements_[i]);
    }
    
    
    inline T& Set(unsigned short i, unsigned short j)
    {
#ifndef NDEBUG
        if((i > 3) || (j > 3))
            throw std::runtime_error("Matrix4<T>::GetElement(unsigned short i,unsigned short j) i or j out of range");
#endif
        return((T&)elements_[4*i+j]);
    }
    
    
    inline T& Set(unsigned short i)
    {
#ifndef NDEBUG
        if((i < 0) || (i > 15))
            throw std::runtime_error("Matrix4<T>::GetElement(unsigned short i) i out of range");
#endif
        return((T&)elements_[i]);
    }
    
    
    inline T Get(unsigned short i) const
    {
#ifndef NDEBUG
        if((i < 0) || (i > 15))
            throw std::runtime_error("Matrix4<T>::GetElement(unsigned short i) i out of range");
#endif
        return(elements_[i]);
    }
    
    
    inline T Get(unsigned short i, unsigned short j) const
    {
#ifndef NDEBUG
        if((i > 3) || (j > 3))
            throw std::runtime_error("Matrix4<T>::GetElement(unsigned short i,unsigned short j) i or j out of range");
#endif
        return(elements_[4*i+j]);
    }
    
    
    inline T* GetArray(){return(elements_);}
    inline const T* GetArray() const {return(elements_);}
    
    inline void SetArray(T* array)
    {
        elements_[0]  = array[0];
        elements_[1]  = array[1];
        elements_[2]  = array[2];
        elements_[3]  = array[3];
        elements_[4]  = array[4];
        elements_[5]  = array[5];
        elements_[6]  = array[6];
        elements_[7]  = array[7];
        elements_[8]  = array[8];
        elements_[9]  = array[9];
        elements_[10] = array[10];
        elements_[11] = array[11];
        elements_[12] = array[12];
        elements_[13] = array[13];
        elements_[14] = array[14];
        elements_[15] = array[15];
    }
    
    
    // ------------------- Operator -------------------------------
    
    inline bool operator == (const Matrix4<T> &a) const;
    inline Matrix4<T> operator*  (const Matrix4<T> &a) const;
    inline Matrix4<T> operator+  (const Matrix4<T> &a) const;
    inline Matrix4<T> operator-  (const Matrix4<T> &a) const;
    
    inline void operator *= (const Matrix4<T> &a);
    inline void operator += (const Matrix4<T> &a);
    inline void operator -= (const Matrix4<T> &a);
    
    inline Matrix4<T> operator*  (const T &a) const;
    inline Matrix4<T> operator/  (const T &a) const;
    
    inline Vector4<T> operator*  (const Vector4<T> &b) const;
    
    inline void operator *= (const T &a);
    inline void operator /= (const T &a);
    
    // -------------------------------------------------------------
    
    inline Vector4<T> GetRow(unsigned short i) const;
    
    inline Vector4<T> GetCol(unsigned short i) const;
    inline void SetRow(unsigned short i, Vector4<T> b);
    inline void SetCol(unsigned short i, Vector4<T> b);
    
    inline bool Invert();
    
    inline Matrix4<T> GetInverse() const;
    inline void Transpose();
    
    inline Matrix4<T> GetTranspose() const;
    
    inline void SetToZero();
    inline void SetToIdentityMatrix();
    inline T Det() const;
    
    inline T Trace() const {return(elements_[0] + elements_[4] + elements_[8] + elements_[15]);}
    inline T Invariant1() const {return(Trace());}
    inline T Invariant2() const {return(0.5 * (Trace()*Trace() - (*this * *this).Trace()));}
    inline T Invariant3() const {return(Det());}
    
    static Matrix4<T> Identity()
    {
        return Matrix4<T>(1,0,0,0,
                          0,1,0,0,
                          0,0,1,0,
                          0,0,0,1);
    }
    
    friend std::ostream& operator<<(std::ostream& str,const Matrix4<T> &a)
    {
        str << a.elements_[0]  << " " << a.elements_[1]  << " " << a.elements_[2]  << " " << a.elements_[3]  << "\n"
        << a.elements_[4]  << " " << a.elements_[5]  << " " << a.elements_[6]  << " " << a.elements_[7]  << "\n"
        << a.elements_[8]  << " " << a.elements_[9]  << " " << a.elements_[10] << " " << a.elements_[11] << "\n"
        << a.elements_[12] << " " << a.elements_[13] << " " << a.elements_[14] << " " << a.elements_[15] << "\n";
        return(str);
    }
    
    void Print() const;
    
    template<class C>
    Matrix4<C> GetConvertedMatrix() const
    {
        Matrix4<C> mat;
        
        for(unsigned int i = 0; i < 4; i++)
            for(unsigned int j = 0; j < 4; j++)
                mat(i,j) = static_cast<C>(Get(i,j));
        
        return(mat);
    }
    
    
protected:
private:
    T elements_[16];
};

template<typename T>
inline bool Matrix4<T>::operator == (const Matrix4<T> &a) const
{
    for(int i = 0; i < 16; i++)
        if(elements_[i] != a.elements_[i])
            return(false);
    
    return(true);
}

template<typename T>
inline Matrix4<T> Matrix4<T>::operator+(const Matrix4<T> &b) const
{
    return(Matrix4<T>(elements_[0]  + b.elements_[0],
                      elements_[1]  + b.elements_[1],
                      elements_[2]  + b.elements_[2],
                      elements_[3]  + b.elements_[3],
                      elements_[4]  + b.elements_[4],
                      elements_[5]  + b.elements_[5],
                      elements_[6]  + b.elements_[6],
                      elements_[7]  + b.elements_[7],
                      elements_[8]  + b.elements_[8],
                      elements_[9]  + b.elements_[9],
                      elements_[10] + b.elements_[10],
                      elements_[11] + b.elements_[11],
                      elements_[12] + b.elements_[12],
                      elements_[13] + b.elements_[13],
                      elements_[14] + b.elements_[14],
                      elements_[15] + b.elements_[15]
                      ));
}

template<typename T>
inline Matrix4<T> Matrix4<T>::operator-(const Matrix4<T> &b) const
{
    return(Matrix4<T>(elements_[0]  - b.elements_[0],
                      elements_[1]  - b.elements_[1],
                      elements_[2]  - b.elements_[2],
                      elements_[3]  - b.elements_[3],
                      elements_[4]  - b.elements_[4],
                      elements_[5]  - b.elements_[5],
                      elements_[6]  - b.elements_[6],
                      elements_[7]  - b.elements_[7],
                      elements_[8]  - b.elements_[8],
                      elements_[9]  - b.elements_[9],
                      elements_[10] - b.elements_[10],
                      elements_[11] - b.elements_[11],
                      elements_[12] - b.elements_[12],
                      elements_[13] - b.elements_[13],
                      elements_[14] - b.elements_[14],
                      elements_[15] - b.elements_[15]
                      ));
}

template<typename T>
inline Matrix4<T> Matrix4<T>::operator*(const Matrix4<T> &b) const
{
    //                        0  1  2  3
    //                        4  5  6  7
    //                        8  9 10 11
    //                       12 13 14 15
    //         0  1  2  3
    //         4  5  6  7
    //         8  9 10 11
    //        12 13 14 15
    
    return(Matrix4<T>(elements_[0]   * b.elements_[0] + elements_[1]  * b.elements_[4] + elements_[2] * b.elements_[8]   + elements_[3]  * b.elements_[12],
                      elements_[0]   * b.elements_[1] + elements_[1]  * b.elements_[5] + elements_[2] * b.elements_[9]   + elements_[3]  * b.elements_[13],
                      elements_[0]   * b.elements_[2] + elements_[1]  * b.elements_[6] + elements_[2] * b.elements_[10]  + elements_[3]  * b.elements_[14],
                      elements_[0]   * b.elements_[3] + elements_[1]  * b.elements_[7] + elements_[2] * b.elements_[11]  + elements_[3]  * b.elements_[15],
                      
                      elements_[4]   * b.elements_[0] + elements_[5]  * b.elements_[4] + elements_[6] * b.elements_[8]   + elements_[7]  * b.elements_[12],
                      elements_[4]   * b.elements_[1] + elements_[5]  * b.elements_[5] + elements_[6] * b.elements_[9]   + elements_[7]  * b.elements_[13],
                      elements_[4]   * b.elements_[2] + elements_[5]  * b.elements_[6] + elements_[6] * b.elements_[10]  + elements_[7]  * b.elements_[14],
                      elements_[4]   * b.elements_[3] + elements_[5]  * b.elements_[7] + elements_[6] * b.elements_[11]  + elements_[7]  * b.elements_[15],
                      
                      elements_[8]   * b.elements_[0] + elements_[9]  * b.elements_[4] + elements_[10] * b.elements_[8]  + elements_[11] * b.elements_[12],
                      elements_[8]   * b.elements_[1] + elements_[9]  * b.elements_[5] + elements_[10] * b.elements_[9]  + elements_[11] * b.elements_[13],
                      elements_[8]   * b.elements_[2] + elements_[9]  * b.elements_[6] + elements_[10] * b.elements_[10] + elements_[11] * b.elements_[14],
                      elements_[8]   * b.elements_[3] + elements_[9]  * b.elements_[7] + elements_[10] * b.elements_[11] + elements_[11] * b.elements_[15],
                      
                      elements_[12]  * b.elements_[0] + elements_[13] * b.elements_[4] + elements_[14] * b.elements_[8]  + elements_[15] * b.elements_[12],
                      elements_[12]  * b.elements_[1] + elements_[13] * b.elements_[5] + elements_[14] * b.elements_[9]  + elements_[15] * b.elements_[13],
                      elements_[12]  * b.elements_[2] + elements_[13] * b.elements_[6] + elements_[14] * b.elements_[10] + elements_[15] * b.elements_[14],
                      elements_[12]  * b.elements_[3] + elements_[13] * b.elements_[7] + elements_[14] * b.elements_[11] + elements_[15] * b.elements_[15]
                      ));
}

template<typename T>
inline Vector4<T>  Matrix4<T>::operator*(const Vector4<T> &b) const
{
    return(Vector4<T>(elements_[0]  * b.elements_[0] + elements_[1]  * b.elements_[1] + elements_[2]  * b.elements_[2] + elements_[3]  * b.elements_[3],
                      elements_[4]  * b.elements_[0] + elements_[5]  * b.elements_[1] + elements_[6]  * b.elements_[2] + elements_[7]  * b.elements_[3],
                      elements_[8]  * b.elements_[0] + elements_[9]  * b.elements_[1] + elements_[10] * b.elements_[2] + elements_[11] * b.elements_[3],
                      elements_[12] * b.elements_[0] + elements_[13] * b.elements_[1] + elements_[14] * b.elements_[2] + elements_[15] * b.elements_[3]
                      ));
}

template<typename T>
inline void Matrix4<T>::operator += (const Matrix4<T> &b)
{
    elements_[0] +=  b.elements_[0];
    elements_[1] +=  b.elements_[1];
    elements_[2] +=  b.elements_[2];
    elements_[3] +=  b.elements_[3];
    
    elements_[4] +=  b.elements_[4];
    elements_[5] +=  b.elements_[5];
    elements_[6] +=  b.elements_[6];
    elements_[7] +=  b.elements_[7];
    
    elements_[8]  +=  b.elements_[8];
    elements_[9]  +=  b.elements_[9];
    elements_[10] +=  b.elements_[10];
    elements_[11] +=  b.elements_[11];
    
    elements_[12] +=  b.elements_[12];
    elements_[13] +=  b.elements_[13];
    elements_[14] +=  b.elements_[14];
    elements_[15] +=  b.elements_[15];
}

template<typename T>
inline void Matrix4<T>::operator -= (const Matrix4<T> &b)
{
    elements_[0] -=  b.elements_[0];
    elements_[1] -=  b.elements_[1];
    elements_[2] -=  b.elements_[2];
    
    elements_[3] -=  b.elements_[3];
    elements_[4] -=  b.elements_[4];
    elements_[5] -=  b.elements_[5];
    
    elements_[6] -=  b.elements_[6];
    elements_[7] -=  b.elements_[7];
    elements_[8] -=  b.elements_[8];
}

template<typename T>
inline void Matrix4<T>::operator *= (const Matrix4<T> &b)
{
    *this *= b;
}

template<typename T>
inline Matrix4<T> Matrix4<T>::operator*(const T &b) const
{
    return(Matrix4<T>(elements_[0]  * b,
                      elements_[1]  * b,
                      elements_[2]  * b,
                      elements_[3]  * b,
                      
                      elements_[4]  * b,
                      elements_[5]  * b,
                      elements_[6]  * b,
                      elements_[7]  * b,
                      
                      elements_[8]  * b,
                      elements_[9]  * b,
                      elements_[10] * b,
                      elements_[11] * b,
                      
                      elements_[12] * b,
                      elements_[13] * b,
                      elements_[14] * b,
                      elements_[15] * b
                      ));
}

template<typename T>
inline Matrix4<T> Matrix4<T>::operator/(const T &b) const
{
    return(Matrix4<T>(elements_[0]  / b,
                      elements_[1]  / b,
                      elements_[2]  / b,
                      elements_[3]  / b,
                      
                      elements_[4]  / b,
                      elements_[5]  / b,
                      elements_[6]  / b,
                      elements_[7]  / b,
                      
                      elements_[8]  / b,
                      elements_[9]  / b,
                      elements_[10] / b,
                      elements_[11] / b,
                      
                      elements_[12] / b,
                      elements_[13] / b,
                      elements_[14] / b,
                      elements_[15] / b
                      ));
}

template<typename T>
inline Matrix4<T> operator*(const T &b, Matrix4<T> a)
{
    return(a*b);
}

template<typename T>
inline Matrix4<T> operator/(const T &b, Matrix4<T> a)
{
    return(a/b);
}

template<typename T>
inline void Matrix4<T>::operator *= (const T &b)
{
    elements_[0] *= b;
    elements_[1] *= b;
    elements_[2] *= b;
    elements_[3] *= b;
    
    elements_[4] *= b;
    elements_[5] *= b;
    elements_[6] *= b;
    elements_[7] *= b;
    
    elements_[8]  *= b;
    elements_[9]  *= b;
    elements_[10] *= b;
    elements_[11] *= b;
    
    elements_[12] *= b;
    elements_[13] *= b;
    elements_[14] *= b;
    elements_[15] *= b;
}

template<typename T>
inline void Matrix4<T>::operator /= (const T &b)
{
    elements_[0] /= b;
    elements_[1] /= b;
    elements_[2] /= b;
    elements_[3] /= b;
    
    elements_[4] /= b;
    elements_[5] /= b;
    elements_[6] /= b;
    elements_[7] /= b;
    
    elements_[8]  /= b;
    elements_[9]  /= b;
    elements_[10] /= b;
    elements_[11] /= b;
    
    elements_[12] /= b;
    elements_[13] /= b;
    elements_[14] /= b;
    elements_[15] /= b;
}

template<typename T>
inline Vector4<T> Matrix4<T>::GetRow(unsigned short i) const
{
#ifndef NDEBUG
    if(i > 3)
        throw std::runtime_error("Matrix4::GetRow(unsigned short i) i out of range");
#endif
    return(Vector4<T>(elements_[4*i], elements_[4*i+1], elements_[4*i+2],elements_[4*i+3]));
}


template<typename T>
inline Vector4<T> Matrix4<T>::GetCol(unsigned short i) const
{
#ifndef NDEBUG
    if(i > 3)
        throw std::runtime_error("Matrix4<T>::GetCol(unsigned short i) i out of range");
#endif
    return(Vector4<T>(elements_[i], elements_[i+4], elements_[i+8],elements_[i+12]));
}


template<typename T>
inline void Matrix4<T>::SetRow(unsigned short i, Vector4<T> b)
{
#ifndef NDEBUG
    if((i < 0) || (i > 3))
        throw std::runtime_error("Matrix4<T>::SetRow(unsigned short i) i out of range");
#endif
    elements_[4*i]   = b.elements_[0];
    elements_[4*i+1] = b.elements_[1];
    elements_[4*i+2] = b.elements_[2];
    elements_[4*i+3] = b.elements_[3];
}


template<typename T>
inline void Matrix4<T>::SetCol(unsigned short i, Vector4<T> b)
{
#ifndef NDEBUG
    if((i < 0) || (i > 3))
        throw std::runtime_error("Matrix4<T>::SetCol(unsigned short i) i out of range");
#endif
    elements_[i]    = b.elements_[0];
    elements_[i+4]  = b.elements_[1];
    elements_[i+8]  = b.elements_[2];
    elements_[i+12] = b.elements_[3];
}


template<typename T>
inline T Matrix4<T>::Det() const
{
    T d[4];
    
    d[0] =  elements_[5] * elements_[10] * elements_[15] - elements_[5] * elements_[11] * elements_[14] - elements_[9] * elements_[6] * elements_[15] + elements_[9] * elements_[7] * elements_[14] + elements_[13] * elements_[6] * elements_[11] - elements_[13] * elements_[7] * elements_[10];
    d[1] = -elements_[4] * elements_[10] * elements_[15] + elements_[4] * elements_[11] * elements_[14] + elements_[8] * elements_[6] * elements_[15] - elements_[8] * elements_[7] * elements_[14] - elements_[12] * elements_[6] * elements_[11] + elements_[12] * elements_[7] * elements_[10];
    d[2] =  elements_[4] * elements_[9]  * elements_[15] - elements_[4] * elements_[11] * elements_[13] - elements_[8] * elements_[5] * elements_[15] + elements_[8] * elements_[7] * elements_[13] + elements_[12] * elements_[5] * elements_[11] - elements_[12] * elements_[7] * elements_[9];
    d[3] = -elements_[4] * elements_[9]  * elements_[14] + elements_[4] * elements_[10] * elements_[13] + elements_[8] * elements_[5] * elements_[14] - elements_[8] * elements_[6] * elements_[13] - elements_[12] * elements_[5] * elements_[10] + elements_[12] * elements_[6] * elements_[9];
    
    return(elements_[0] * d[0] + elements_[1] * d[1] + elements_[2] * d[2] + elements_[3] * d[3]);
}


template<typename T>
inline bool Matrix4<T>::Invert()
{
    T i[16];
    
    i[0]  =  elements_[5] * elements_[10] * elements_[15] - elements_[5] * elements_[11] * elements_[14] - elements_[9] * elements_[6] * elements_[15] + elements_[9] * elements_[7] * elements_[14] + elements_[13] * elements_[6] * elements_[11] - elements_[13] * elements_[7] * elements_[10];
    i[1]  = -elements_[1] * elements_[10] * elements_[15] + elements_[1] * elements_[11] * elements_[14] + elements_[9] * elements_[2] * elements_[15] - elements_[9] * elements_[3] * elements_[14] - elements_[13] * elements_[2] * elements_[11] + elements_[13] * elements_[3] * elements_[10];
    i[2]  =  elements_[1] * elements_[6]  * elements_[15] - elements_[1] * elements_[7]  * elements_[14] - elements_[5] * elements_[2] * elements_[15] + elements_[5] * elements_[3] * elements_[14] + elements_[13] * elements_[2] * elements_[7]  - elements_[13] * elements_[3] * elements_[6];
    i[3]  = -elements_[1] * elements_[6]  * elements_[11] + elements_[1] * elements_[7]  * elements_[10] + elements_[5] * elements_[2] * elements_[11] - elements_[5] * elements_[3] * elements_[10] - elements_[9]  * elements_[2] * elements_[7]  + elements_[9]  * elements_[3] * elements_[6];
    i[4]  = -elements_[4] * elements_[10] * elements_[15] + elements_[4] * elements_[11] * elements_[14] + elements_[8] * elements_[6] * elements_[15] - elements_[8] * elements_[7] * elements_[14] - elements_[12] * elements_[6] * elements_[11] + elements_[12] * elements_[7] * elements_[10];
    i[5]  =  elements_[0] * elements_[10] * elements_[15] - elements_[0] * elements_[11] * elements_[14] - elements_[8] * elements_[2] * elements_[15] + elements_[8] * elements_[3] * elements_[14] + elements_[12] * elements_[2] * elements_[11] - elements_[12] * elements_[3] * elements_[10];
    i[6]  = -elements_[0] * elements_[6]  * elements_[15] + elements_[0] * elements_[7]  * elements_[14] + elements_[4] * elements_[2] * elements_[15] - elements_[4] * elements_[3] * elements_[14] - elements_[12] * elements_[2] * elements_[7]  + elements_[12] * elements_[3] * elements_[6];
    i[7]  =  elements_[0] * elements_[6]  * elements_[11] - elements_[0] * elements_[7]  * elements_[10] - elements_[4] * elements_[2] * elements_[11] + elements_[4] * elements_[3] * elements_[10] + elements_[8]  * elements_[2] * elements_[7]  - elements_[8]  * elements_[3] * elements_[6];
    i[8]  =  elements_[4] * elements_[9]  * elements_[15] - elements_[4] * elements_[11] * elements_[13] - elements_[8] * elements_[5] * elements_[15] + elements_[8] * elements_[7] * elements_[13] + elements_[12] * elements_[5] * elements_[11] - elements_[12] * elements_[7] * elements_[9];
    i[9]  = -elements_[0] * elements_[9]  * elements_[15] + elements_[0] * elements_[11] * elements_[13] + elements_[8] * elements_[1] * elements_[15] - elements_[8] * elements_[3] * elements_[13] - elements_[12] * elements_[1] * elements_[11] + elements_[12] * elements_[3] * elements_[9];
    i[10] =  elements_[0] * elements_[5]  * elements_[15] - elements_[0] * elements_[7]  * elements_[13] - elements_[4] * elements_[1] * elements_[15] + elements_[4] * elements_[3] * elements_[13] + elements_[12] * elements_[1] * elements_[7]  - elements_[12] * elements_[3] * elements_[5];
    i[11] = -elements_[0] * elements_[5]  * elements_[11] + elements_[0] * elements_[7]  * elements_[9]  + elements_[4] * elements_[1] * elements_[11] - elements_[4] * elements_[3] * elements_[9]  - elements_[8]  * elements_[1] * elements_[7]  + elements_[8]  * elements_[3] * elements_[5];
    i[12] = -elements_[4] * elements_[9]  * elements_[14] + elements_[4] * elements_[10] * elements_[13] + elements_[8] * elements_[5] * elements_[14] - elements_[8] * elements_[6] * elements_[13] - elements_[12] * elements_[5] * elements_[10] + elements_[12] * elements_[6] * elements_[9];
    i[13] =  elements_[0] * elements_[9]  * elements_[14] - elements_[0] * elements_[10] * elements_[13] - elements_[8] * elements_[1] * elements_[14] + elements_[8] * elements_[2] * elements_[13] + elements_[12] * elements_[1] * elements_[10] - elements_[12] * elements_[2] * elements_[9];
    i[14] = -elements_[0] * elements_[5]  * elements_[14] + elements_[0] * elements_[6]  * elements_[13] + elements_[4] * elements_[1] * elements_[14] - elements_[4] * elements_[2] * elements_[13] - elements_[12] * elements_[1] * elements_[6]  + elements_[12] * elements_[2] * elements_[5];
    i[15] =  elements_[0] * elements_[5]  * elements_[10] - elements_[0] * elements_[6]  * elements_[9]  - elements_[4] * elements_[1] * elements_[10] + elements_[4] * elements_[2] * elements_[9]  + elements_[8]  * elements_[1] * elements_[6]  - elements_[8]  * elements_[2] * elements_[5];
    
    T det = elements_[0] * i[0] + elements_[1] * i[4] + elements_[2] * i[8] + elements_[3] * i[12];
    
    if(det == 0)
    {
        *this = Matrix4<T>(INFINITY,INFINITY,INFINITY,INFINITY,INFINITY,INFINITY,INFINITY,INFINITY,INFINITY,INFINITY,INFINITY,INFINITY,INFINITY,INFINITY,INFINITY,INFINITY);
        return(false);
    }
    
    for(int j = 0; j < 16; j++)
        elements_[j] = i[j]/det;
    
    return(true);
}


template<typename T>
inline Matrix4<T> Matrix4<T>::GetInverse() const
{
    Matrix4<T> m = *this;
    m.Invert();
    return(m);
}


template<typename T>
inline void Matrix4<T>::Transpose()
{
    T temp;
    
    // 0  1  2  3
    // 4  5  6  7
    // 8  9  10 11
    // 12 13 14 15
    
    temp         = elements_[4];
    elements_[4] = elements_[1];
    elements_[1] = temp;
    
    temp         = elements_[8];
    elements_[8] = elements_[2];
    elements_[2] = temp;
    
    temp          = elements_[12];
    elements_[12] = elements_[3];
    elements_[3]  = temp;
    
    temp         = elements_[9];
    elements_[9] = elements_[6];
    elements_[6] = temp;
    
    temp          = elements_[13];
    elements_[13] = elements_[7];
    elements_[7]  = temp;
    
    temp          = elements_[14];
    elements_[14] = elements_[11];
    elements_[11] = temp;
}


template<typename T>
inline Matrix4<T> Matrix4<T>::GetTranspose() const
{
    return(Matrix4<double>(elements_[0], elements_[4], elements_[8], elements_[12],
                           elements_[1], elements_[5], elements_[9], elements_[13],
                           elements_[2], elements_[6], elements_[10],elements_[14],
                           elements_[3], elements_[7], elements_[11],elements_[15]));
}


template<typename T>
inline void Matrix4<T>::SetToZero()
{
    elements_[0]  = 0; elements_[1]  = 0; elements_[2]  = 0; elements_[3]  = 0;
    elements_[4]  = 0; elements_[5]  = 0; elements_[6]  = 0; elements_[7]  = 0;
    elements_[8]  = 0; elements_[9]  = 0; elements_[10] = 0; elements_[11] = 0;
    elements_[12] = 0; elements_[13] = 0; elements_[14] = 0; elements_[15] = 0;
}


template<typename T>
inline void Matrix4<T>::SetToIdentityMatrix()
{
    elements_[0]  = 1; elements_[1]  = 0; elements_[2]  = 0; elements_[3]  = 0;
    elements_[4]  = 0; elements_[5]  = 1; elements_[6]  = 0; elements_[7]  = 0;
    elements_[8]  = 0; elements_[9]  = 0; elements_[10] = 1; elements_[11] = 0;
    elements_[12] = 0; elements_[13] = 0; elements_[14] = 0; elements_[15] = 1;
}


template<typename T>
inline Matrix4<T> DyadicProduct(const Vector4<T>& a,const Vector4<T>& b)
{
    return(Matrix4<T>(a(0)*b(0), a(0)*b(1), a(0)*b(2), a(0)*b(3),
                      a(1)*b(0), a(1)*b(1), a(1)*b(2), a(1)*b(3),
                      a(2)*b(0), a(2)*b(1), a(2)*b(2), a(2)*b(3),
                      a(3)*b(0), a(3)*b(1), a(3)*b(3), a(3)*b(3)
                      ));
}


template<typename T>
void Matrix4<T>::Print() const
{
    std::cout << *this;
}
}
#endif
