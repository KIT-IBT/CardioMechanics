//
//  MathPack.h
//  CardioMechanics
//
//  Created by Thomas Fritz on 12.12.11.
//  Copyright (c) 2011 Institute for Biomedical Engineering, Karlsruhe Institute of Technology (KIT). All rights reserved.
//

#ifndef MATHPACK_TRIANGLE_H
#define MATHPACK_TRIANGLE_H

#include "Vector3.h"

namespace math_pack
{

template<typename T> class Vector3;

template<typename T>
class Triangle
{
    
public:
    Triangle() 
    {
        v_[0] = Vector3<T>(0,0,0);
        v_[1] = Vector3<T>(0,0,0);
        v_[2] = Vector3<T>(0,0,0);
        n_ = Vector3<T>(0,0,0);
        d_ = 0;
    }
    
    Triangle(Vector3<T> a, Vector3<T> b, Vector3<T> c)
    {
        v_[0] = a;
        v_[1] = b;
        v_[2] = c;
        n_ = CrossProduct((v_[1]-v_[0]),(v_[2]-v_[0]));
        n_.Normalize();
        d_ = n_*v_[0];
    }
    Triangle(T* p)
    {
        v_[0] = Vector3<T>(p);
        p+=3;
        v_[1] = Vector3<T>(p);
        p+=3;
        v_[2] = Vector3<T>(p);
        n_ = CrossProduct((v_[1]-v_[0]),(v_[2]-v_[0]));
        n_.Normalize();
        d_ = n_*v_[0];
    }
    
    void SetNodes(Vector3<T> a, Vector3<T> b, Vector3<T> c)
    {
        v_[0] = a;
        v_[1] = b;
        v_[2] = c;
        n_ = CrossProduct((v_[1]-v_[0]),(v_[2]-v_[0]));
        n_.Normalize();
        d_ = n_*v_[0];
    }
    void SetNodes(T* p)
    {
        v_[0] = Vector3<T>(p);
        p+=3;
        v_[1] = Vector3<T>(p);
        p+=3;
        v_[2] = Vector3<T>(p);
        n_ = CrossProduct((v_[1]-v_[0]),(v_[2]-v_[0]));
        n_.Normalize();
        d_ = n_*v_[0];
    }
    
    void SetNode(int i, const Vector3<T>& p)
    {
        v_[i] = p;
        n_ = CrossProduct((v_[1]-v_[0]),(v_[2]-v_[0]));
        n_.Normalize();
        d_ = n_*v_[0];
    }
    
    Vector3<T> GetNode(int i) const
    {
        return v_[i];
    }
    
    Vector3<T> GetNormalVector() const {return n_;}
    Vector3<T> GetCentroid() const {return ((v_[0] + v_[1] + v_[2])/3.0);}
    Vector3<T> GetGaussPoint3(int i) const
    {
        switch(i)
        {
            case 0:
                return v_[0]*(2.0/3.0) + v_[1]*(1.0/6.0) + v_[2]*(1.0/6.0);
                break;
            case 1:
                return v_[0]*(1.0/6.0) + v_[1]*(2.0/3.0) + v_[2]*(1.0/6.0);
                break;
            case 2:
                return v_[0]*(1.0/6.0) + v_[1]*(1.0/6.0) + v_[2]*(2.0/3.0);
                break;
            default:
                throw std::runtime_error("Triangle<T>::GetGaussPoint3(int i) i is out of range (0-2)");
        }
    }
    
    void GetPlaneEquation(Vector3<T>* p,T* d){*p=n_;*d=d_;} const
    
    T GetDistanceTo(const Vector3<T>& p) const {return fabs( n_*p - d_);}
    Vector3<T> GetFoot(const Vector3<T>& p) const  {return p + (n_*p - d_) * n_;}
    
    Vector3<T> CalcIntersectionPoint(const Vector3<T>& p, const Vector3<T>& q) const
    {
        T t = ((v_[0]-p)*n_) / (q*n_);
        return p + (t * q);
    }
    
    bool IsPointWithinTriangle(const Vector3<T>& p) const
    {
        
        T a = fabs(CrossProduct((v_[0]-p),(v_[1]-p)).Norm()) + fabs(CrossProduct( (v_[1]-p),(v_[2]-p)).Norm()) + fabs(CrossProduct( (v_[0]-p),(v_[2]-p)).Norm());
        T b = 2*GetArea();
        if(fabs(a-b) < 1e-6)
            return true;
        else
            return false;
    }
    inline T GetSignedArea() const  {return 0.5 * CrossProduct((v_[1]-v_[0]),(v_[2]-v_[0])).Norm();}
    T GetArea() const  {return fabs(GetSignedArea());}
    T GetSubArea(const Vector3<T>& p,int i) const
    {
        Vector3<TFloat> f = GetFoot(p);
        
        switch(i)
        {
            case 0 : return Triangle<TFloat>(f,v_[1],v_[2]).GetArea();
                break;
            case 1 : return Triangle<TFloat>(v_[0],f,v_[2]).GetArea();
                break;
            case 2 : return Triangle<TFloat>(v_[0],v_[1],f).GetArea();
                break;
            default:
                return 0;
        }
        
    }
private:
    
    Vector3<T> v_[3];
    Vector3<T> n_;
    T d_;
};
}





#endif
