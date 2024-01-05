/*
 * File: CBTensionModel.h
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


#pragma once

// #include "CBFileManager.h"
// #include "CBElement.h"

#include "DCType.h" // defines TFloat
#include "Matrix3.h"
#include "CBElementSolid.h"
#include "CBDataFromFile.h"
#include "CBFileManager.h"
#include "math.h"
#include <functional>

/// base class for classes implementing functions to calculate tensions / read tension from files
class CBTensionModel {
public:
    CBTensionModel() {}
    
    /// tension is the force in fiber direction (a scalar value)
    virtual double CalcActiveTension(const math_pack::Matrix3<double> &deformation, const double time) = 0;
    
    /// active stress is the three-dimensional PK2 stress matrix resulting from tension
    virtual math_pack::Matrix3<double> CalcActiveStress(const math_pack::Matrix3<double> &deformation,
                                                        const double time) {
        double activeStress = CalcActiveTension(deformation, time);
        
        double I4_f = deformation.GetCol(0) * deformation.GetCol(0);
        double I4_s = deformation.GetCol(1) * deformation.GetCol(1);
        double I4_n = deformation.GetCol(2) * deformation.GetCol(2);
        
        return activeStress *
        Matrix3<TFloat> {1./sqrt(I4_f), 0, 0,  0, 1./sqrt(I4_s), 0,  0, 0, 1./sqrt(I4_n)} *stressCoefficients_;

        // For benchmark problems/examples to be correct
        // return activeStress * Matrix3<TFloat> {1, 0, 0, 0, 0, 0, 0, 0, 0};
    }
    
    virtual CBStatus SetActiveTensionAtQuadraturePoint(int indexQP, TFloat activeTension) {
        return CBStatus::SUCCESS;
    }
    
    virtual void SetfibreRatio(Vector3<TFloat> ffRatio) {}
    
    virtual void SetStressCoefficients(Matrix3<TFloat> stressCoeff) {
        stressCoefficients_ = stressCoeff;
    }
    
    /// helper function: perform one explicit euler step
    static double ExplicitEulerStep(double u, double dt, std::function<double(double)> f) {
        double k1 = f(u);
        
        return u + dt*k1;
    }
    
    /// helper function: performs one 4th order Runge-Kutta step (TODO: something similar exists in CBCircModel.hpp ::Integrate() ? )
    static double RungeKutta4Step(double u, double dt, std::function<double(double)> f) {
        double k1 = f(u);
        double k2 = f(u+dt/2.*k1);
        double k3 = f(u+dt/2.*k2);
        double k4 = f(u + dt*k3);
        
        return u + dt*(k1 + 2.*k2 + 2.*k3 + k4)/6.;
    }
    
    /// Method to properly set the activation time attribute of the element
    void SetActivationTime(TFloat activationTime) {
        activationTime_ = activationTime;
    }
    
    /// delays the activation time by the given factor
    void AddToActivationTime(TFloat activationTime) {
        activationTime_ += activationTime;
    }
    
    /// Method to retrieve final beginning of the activation
    virtual TFloat GetActivationTime() {
        return activationTime_;
    }
    
protected:
    TFloat Tmax_ = 0;
    
private:
    /// contains local activation time and lat-offset (from lat-reader)
    TFloat activationTime_ = 0.0;
    Matrix3<TFloat> stressCoefficients_ = {1, 0, 0, 0, 0, 0, 0, 0, 0};
};

/// A simple class that produces no tension at all. Needed as default case.
class CBNoTension : public CBTensionModel {
public:
    virtual double CalcActiveTension(const math_pack::Matrix3<double> &deformation, const double time) override {
        return 0;
    }
    
    virtual math_pack::Matrix3<double> CalcActiveStress(const math_pack::Matrix3<double> &deformation,
                                                        const double time) override {
        return Matrix3<double>();
    }
    
    double GetActivationTime() override {
        return 0;
    }
};

/// Provides tensions from a tension.list file. The actual data is outsourced and managed by a FileManager object that lives within the solver class. The tension between stored data gets interpolated from CBDataFromFile* object.
class CBFileTension : public CBTensionModel {
protected:
    CBDataFromFile *data_;
    const CBElementSolid *e_;
    
public:
    CBFileTension(const CBFileManager *fileManager, CBElementSolid *e) {
        Tmax_ = e->GetMaterial()->GetProperties()->tensionMax_;
        e_ = e;
        assert(e_ != nullptr);
        data_ = fileManager->GetDataFromFileObject(e);
        assert(data_ != nullptr);
    }
    
    ~CBFileTension() {}
    
    virtual double CalcActiveTension(const math_pack::Matrix3<double> &deformation, const double time) override {
        return Tmax_ * data_->Get(time, e_->GetIndex());
    }
};
