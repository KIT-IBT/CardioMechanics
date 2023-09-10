/*
 *  Vector4.h
 *
 *
 *  Created by Thomas Fritz on 24.08.09.
 *  Copyright 2009 Institute for Biomedical Engineering, Karlsruhe Institute of Technology (KIT). All rights reserved.
 *
 */
#ifndef MATHPACK_VECTOR4_H
#define MATHPACK_VECTOR4_H

#include <cmath>
#include <iomanip>
#include <iostream>
#include <initializer_list>
#include <cassert>
#include <algorithm>
#include "Matrix4.h"


// forward declaration of Matrix3

namespace math_pack {
template<typename T> class Matrix4;

template<typename T>
class Vector4 {
public:
    friend class Matrix4<T>;
    Vector4(T x, T y, T z, T w) {
        elements_[0] = x;
        elements_[1] = y;
        elements_[2] = z;
        elements_[3] = w;
    }
    
    Vector4(const T *const array) {
        elements_[0] = array[0];
        elements_[1] = array[1];
        elements_[2] = array[2];
        elements_[3] = array[3];
    }
    
    Vector4(std::initializer_list<T> l) {
#ifndef NDEBUG
        if (l.size() != 4) throw std::runtime_error(
                                                    "Vector4<T>::Vector4(std::initializer_list l) initializer list has to exactly contain 3 elements ");
#endif // ifndef NDEBUG
        T *e = elements_;
        for (auto i : l) {
            *e = i;
            ++e;
        }
    }
    
    Vector4() {}
    
    inline T *GetArray() {
        return elements_;
    }
    
    inline const T *GetArray() const {
        return elements_;
    }
    
    inline T & operator()(unsigned short i) {
#ifndef NDEBUG
        if (i > 3) throw std::runtime_error("Vector4<T>::GetElement(unsigned short i) i out of range");
#endif // ifndef NDEBUG
        return (T &)elements_[i];
    }
    
    inline const T & operator()(unsigned short i) const {
#ifndef NDEBUG
        if (i > 3) throw std::runtime_error("Vector4<T>::GetElement(unsigned short i) i out of range");
#endif // ifndef NDEBUG
        return (T &)elements_[i];
    }
    
    inline T & Set(unsigned short i) {
#ifndef NDEBUG
        if (i > 3) throw std::runtime_error("Vector4<T>::GetElement(unsigned short i) i out of range");
#endif // ifndef NDEBUG
        return (T &)elements_[i];
    }
    
    inline T Get(unsigned short i) const {
#ifndef NDEBUG
        if (i > 3) throw std::runtime_error("Vector4<T>::GetElement(unsigned short i) i out of range");
#endif // ifndef NDEBUG
        return elements_[i];
    }
    
    inline T & X() {
        return (T &)elements_[0];
    }
    
    inline T & Y() {
        return (T &)elements_[1];
    }
    
    inline T & Z() {
        return (T &)elements_[2];
    }
    
    inline T & W() {
        return (T &)elements_[3];
    }
    
    inline const T & X() const {
        return (T &)elements_[0];
    }
    
    inline const T & Y() const {
        return (T &)elements_[1];
    }
    
    inline const T & Z() const {
        return (T &)elements_[2];
    }
    
    inline const T & W() const {
        return (T &)elements_[3];
    }
    
    inline void SetArray(const T *const array) {
        elements_[0] = array[0];
        elements_[1] = array[1];
        elements_[2] = array[2];
        elements_[3] = array[3];
    }
    
    // ------------------- Operator ----------------------
    
    inline bool operator==(const Vector4<T> &a) const;
    
    inline T operator*(const Vector4<T> &a) const;
    inline Vector4<T> operator+(const Vector4<T> &a) const;
    inline Vector4<T> operator-(const Vector4<T> &a) const;
    
    inline Vector4<T> operator*(const T &b) const;
    inline Vector4<T> operator/(const T &b) const;
    
    inline void operator*=(const T &b);
    inline void operator/=(const T &b);
    inline void operator+=(const Vector4<T> &a);
    inline void operator-=(const Vector4<T> &a);
    
    // ----------------------------------------------------
    
    friend std::ostream & operator<<(std::ostream &str, const Vector4<T> &a) {
        str << a.elements_[0] << "\n" << a.elements_[1] << "\n" << a.elements_[2] << "\n" << a.elements_[4] << "\n";
        return str;
    }
    
    void Print() const;
    inline T Norm() const;
    inline T Max() const;
    inline T Min() const;
    inline bool Normalize();
    std::ostream & Print(std::ostream &outputStream) const;
    
    template<class C> Vector4<C> GetConvertedVector() {
        Vector4<C> vec((C)elements_[0], (C)elements_[1], (C)elements_[2], (C)elements_[3]);
        
        return vec;
    }
    
protected:
private:
    T elements_[4];
};

template<typename T>
inline bool Vector4<T>::operator==(const Vector4<T> &b) const {
    if ((elements_[0] ==  b.elements_[0]) && (elements_[1] ==  b.elements_[1]) && (elements_[2] ==  b.elements_[2]) &&
        (b.elements_[3] == elements_[3]))
        return true;
    else
        return false;
}

template<typename T>
inline Vector4<T> Vector4<T>::operator+(const Vector4<T> &b) const {
    return Vector4<T>(elements_[0] + b.elements_[0], elements_[1] + b.elements_[1], elements_[2] + b.elements_[2],
                      elements_[3] + b.elements_[3]);
}

template<typename T>
inline Vector4<T> Vector4<T>::operator-(const Vector4<T> &b) const {
    return Vector4<T>(elements_[0] - b.elements_[0], elements_[1] - b.elements_[1], elements_[2] - b.elements_[2],
                      elements_[3] - b.elements_[3]);
}

template<typename T>
inline T Vector4<T>::operator*(const Vector4<T> &b) const {
    return elements_[0]*b.elements_[0] + elements_[1]*b.elements_[1] + elements_[2]*b.elements_[2] + elements_[3]*
    b.elements_[3];
}

template<typename T>
inline Vector4<T> Vector4<T>::operator*(const T &b) const {
    return Vector4<T>(elements_[0]*b, elements_[1]*b, elements_[2]*b, elements_[3]*b);
}

template<typename T>
inline Vector4<T> operator*(const T &b, const Vector4<T> &a) {
    return a*b;
}

template<typename T>
inline Vector4<T> Vector4<T>::operator/(const T &b) const {
    return Vector4<T>(elements_[0]/b, elements_[1]/b, elements_[2]/b, elements_[3]/b);
}

template<typename T>
inline Vector4<T> operator/(const T &b, const Vector4<T> &a) {
    return a/b;
}

template<typename T>
inline void Vector4<T>::operator+=(const Vector4<T> &b) {
    elements_[0] +=  b.elements_[0];
    elements_[1] +=  b.elements_[1];
    elements_[2] +=  b.elements_[2];
    elements_[3] +=  b.elements_[3];
}

template<typename T>
inline void Vector4<T>::operator-=(const Vector4<T> &b) {
    elements_[0] -=  b.elements_[0];
    elements_[1] -=  b.elements_[1];
    elements_[2] -=  b.elements_[2];
    elements_[3] -=  b.elements_[3];
}

template<typename T>
inline T Vector4<T>::Norm() const {
    return sqrt(
                elements_[0]*elements_[0] + elements_[1]*elements_[1] + elements_[2]*elements_[2]+elements_[3]*elements_[3]);
}

template<typename T>
inline T Vector4<T>::Max() const {
    return *std::max_element(elements_, elements_ + 4);
}

template<typename T>
inline T Vector4<T>::Min() const {
    return *std::min_element(elements_, elements_ + 4);
}

template<typename T>
inline bool Vector4<T>::Normalize() {
    T norm = Norm();
    
    if (norm == T(0.0)) {
#ifndef NDEBUG
        throw std::runtime_error("Vector4<T><T>::Normalize() Division by Zero");
#endif // ifndef NDEBUG
        return false;
    }
    
    elements_[0] /= norm;
    elements_[1] /= norm;
    elements_[2] /= norm;
    elements_[3] /= norm;
    return true;
}

template<typename T>
inline std::ostream & Vector4<T>::Print(std::ostream &outputStream) const {
    return outputStream << std::scientific << elements_[0] << " " << elements_[1] << " " << elements_[2] << elements_[3];
}

template<typename T>
void Vector4<T>::Print() const {
    std::cout << *this;
}
}
#endif // ifndef MATHPACK_VECTOR4_H
