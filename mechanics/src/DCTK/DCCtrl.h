/*
 * File: DCCtrl.h
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


#ifndef DC_CTRL_H
#define DC_CTRL_H

#include <cmath>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <stdexcept>

#include "DCType.h"

class DCVectorBase;
class DCMatrixBase;
class DCVector;
class DCMatrix;


class DCCtrl
{
public:
    
    friend class DCMatrix;
    friend class DCVector;
    
    virtual ~DCCtrl(){}
    
    static int Init(int argc, char** argv);
    static bool IsReady(){if(instance_)return true;else return false;}
    
    static int Finalize()
    {
        DCCtrl::instance_->logFile_.close();
        delete instance_;
        return 0;
    }
    
    
    static unsigned int GetProcessID();
    static unsigned int GetNumberOfProcesses();
    
    static bool IsProcessZero();
    static bool IsParallel();
    
    template <class T>
    static void GatherToZero(std::vector<T>& a);
    template <class T>
    static void CopyFromZeroToAll(std::vector<T>& a);
    template <class T>
    static void WeightedAverage(T localAverage, T localWeight, T& globalAverage);
    
    // Specializations
    
    TFloat FloatType; // Just used to decltype the template paramters
    TInt   IntType;   // Just used to decltype the template paramters
    
    class Print
    {
    public:
        template <class T>
        friend Print& operator<<(Print& me, T s)
        {
            if(!DCCtrl::IsReady())
                throw std::runtime_error("Print& Messages::operator<<(Print& me,T str): DCCtrl has not been intialized");
            std::stringstream ss;
            ss << s;
            me.str_ += ss.str();
            me.Flush();
            return me;
        }
        typedef std::ostream& (*ManipPtr)(std::ostream&);
        friend Print& operator<<(Print& me,ManipPtr s)
        {
            if(!DCCtrl::IsReady())
                throw std::runtime_error("Print& Print::operator<<(Print& me,T str): DCCtrl has not been intialized");
            
            if(s == static_cast<ManipPtr>(std::endl))
            {
                me.str_ += "\n";
            }
            if(s == static_cast<ManipPtr>(std::flush))
                me.Flush();
            
            return me;
        }
        void Flush()
        {
            if(DCCtrl::IsProcessZero())
            {
                if( (this == &DCCtrl::print || this == &DCCtrl::cprint) ||
                   ( (this == &DCCtrl::verbose || this == &DCCtrl::cverbose) && (DCCtrl::instance_->isVerbose_) ) ||
                   ( (this == &DCCtrl::debug || this == &DCCtrl::cdebug) && (DCCtrl::instance_->isDebug_)) )
                {
                    std::cout << str_ << std::flush;
                    if(DCCtrl::instance_->logFile_.good())
                    {
                        std::size_t i = 0;
                        while(i != std::string::npos)
                        {
                            i = str_.find("\xd");
                            if(i != std::string::npos)
                            {
                                str_.erase(i, 1);
                                str_.insert(i,"\n");
                            }
                        }
                        DCCtrl::instance_->logFile_ << log_ << str_;
                        DCCtrl::instance_->logFile_.flush();
                        log_.clear();
                    }
                    else
                        log_ += str_;
                    
                }
            }
            str_.clear();
        }
    protected:
        std::string str_;
        std::string log_;
    };
    
    class CPrint : public Print
    {
        template<class T>
        friend CPrint& operator<<(CPrint& me, T s)
        {
            if(!DCCtrl::IsReady())
                throw std::runtime_error("Print& Messages::operator<<(Print& me,T str): DCCtrl has not been intialized");
            std::stringstream ss;
            ss << s;
            me.str_ += ss.str();
            return me;
        }
        
        friend CPrint& operator<<(CPrint& me,ManipPtr s)
        {
            
            if(!DCCtrl::IsReady())
                throw std::runtime_error("Print& Print::operator<<(Print& me,T str): DCCtrl has not been intialized");
            
            if(s == static_cast<ManipPtr>(std::endl))
                me.str_ += "\n";
            
            if(s == static_cast<ManipPtr>(std::flush))
            {
                std::vector<std::string> m = {"Process " + std::to_string(DCCtrl::GetProcessID()) + ": " + me.str_};
                DCCtrl::GatherToZero(m);
                me.str_.clear();
                for(auto s: m)
                    me.str_ += (s + "\n");
                me.Flush();
            }
            
            return me;
        }
        
    };
    
    static void SetLogFile(std::string filename);
    
    static Print print;   // only process zero prints these messages
    static Print verbose;
    static Print debug;
    
    static CPrint cprint;   // only process zero prints these messages
    static CPrint cverbose;
    static CPrint cdebug;
    
    bool isVerbose_ = false;
    bool isDebug_ = false;
    
    void static VerboseOn();
    void static DebugOn();
    void static VerboseOff();
    void static DebugOff();
    bool static GetVerbose();
    bool static GetDebug();
    
protected:
    
    DCCtrl(){};
    
private:
    
    static DCVectorBase* NewVector();
    static DCMatrixBase* NewMatrix();
    
    virtual DCVectorBase* NewVectorImpl() = 0;
    virtual DCMatrixBase* NewMatrixImpl() = 0;
    
    virtual int InitImpl(int argc, char** argv) = 0;
    
    virtual unsigned int GetProcessIDImpl() = 0;
    virtual unsigned int GetNumberOfProcessesImpl() = 0;
    
    virtual bool IsProcessZeroImpl() = 0;
    virtual bool IsParallelImpl() = 0;
    
    virtual void GatherToZeroImpl(std::vector<int>& a)    = 0;
    virtual void GatherToZeroImpl(std::vector<long>& a)   = 0;
    virtual void GatherToZeroImpl(std::vector<float>& a)  = 0;
    virtual void GatherToZeroImpl(std::vector<double>& a) = 0;
    virtual void GatherToZeroImpl(std::vector<std::string>& a) = 0;
    
    virtual void CopyFromZeroToAllImpl(std::vector<int>& a)    = 0;
    virtual void CopyFromZeroToAllImpl(std::vector<long>& a)   = 0;
    virtual void CopyFromZeroToAllImpl(std::vector<float>& a)  = 0;
    virtual void CopyFromZeroToAllImpl(std::vector<double>& a) = 0;
    
    virtual void WeightedAverageImpl(double& localAverage, double& localWeight, double& globalAverage) = 0;
    
    std::ofstream logFile_;
    static DCCtrl* instance_;
    std::vector<long> startTime_;
};

// ----- Gather to Zero -----

template <class T>
void DCCtrl::GatherToZero(std::vector<T>& a)
{
    if(instance_)
        instance_->GatherToZeroImpl(a);
    else
        throw std::runtime_error("class DCCtrl::GatherToZero(std::vector<T>& a): DCCtrl is not yet inizialized");
}

// ----- ----- ----- ----- ----- -----


// ----- Copy from zero to all -----

template <class T>
void DCCtrl::CopyFromZeroToAll(std::vector<T>& a)
{
    if(instance_)
        instance_->CopyFromZeroToAllImpl(a);
    else
        throw std::runtime_error(" class DCCtrl::CopyFromZeroToAll(std::vector<T>& a): DCCtrl is not yet inizialized");
}

// ----- ----- ----- ----- ----- -----
// ----- -- Weighted average --- -----

template <class T>
void DCCtrl::WeightedAverage(T localAverage, T localWeight, T& globalAverage)
{
    if(instance_)
        instance_->WeightedAverageImpl(localAverage, localWeight, globalAverage);
    else
        throw std::runtime_error(" class void DCCtrl::WeightedAverage(T& localAverage, T& localWeight, T& globalAverage): DCCtrl is not yet inizialized");
}



#include "DCVector.h"
#include "DCMatrix.h"
#include "DCCtrlPETSc.h"

#endif
