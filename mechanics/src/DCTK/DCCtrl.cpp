//
//  DCCtrl.cpp
//  CardioMechanics
//
//  Created by Thomas Fritz on 30.10.11.
//  Copyright (c) 2011 Institute for Biomedical Engineering, Karlsruhe Institute of Technology (KIT). All rights reserved.
//


#include <iostream>
#include <vector>
#include <stdexcept>
#include <chrono>
#include "filesystem.h"

#include "DCCtrlPETSc.h"
#include "DCVector.h"



DCCtrl* DCCtrl::instance_ = 0;
DCCtrl::Print DCCtrl::print;   // only process zero prints these messages
DCCtrl::CPrint DCCtrl::cprint;  // messages which are collected from all processes

DCCtrl::Print DCCtrl::verbose;
DCCtrl::CPrint DCCtrl::cverbose;

DCCtrl::Print DCCtrl::debug;
DCCtrl::CPrint DCCtrl::cdebug;

int DCCtrl::Init(int argc, char** argv)
{
    if(instance_ == 0)
        instance_ = new DCCtrlPETSc;
    instance_->InitImpl(argc,argv);
    instance_->startTime_ = {std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now().time_since_epoch()).count()};
    GatherToZero(instance_->startTime_);
    return 0;
}

unsigned int DCCtrl::GetProcessID()
{
    if(instance_)
        return instance_->GetProcessIDImpl();
    else
        throw std::runtime_error("unsigned int DCCtrl::GetProcessID(): DCCtrl is not yet inizialized");
}


unsigned int DCCtrl::GetNumberOfProcesses()
{
    if(instance_)
        return instance_->GetNumberOfProcessesImpl();
    else
        throw std::runtime_error("unsigned int DCCtrl::GetNumberOfProcesses(): DCCtrl is not yet inizialized");
}


bool DCCtrl::IsProcessZero()
{
    if(instance_)
        return instance_->IsProcessZeroImpl();
    else
        throw std::runtime_error("bool DCCtrl::IsProcessZero(): DCCtrl is not yet inizialized");
}



bool DCCtrl::IsParallel()
{
    if(instance_)
        return instance_->IsParallelImpl();
    else
        throw std::runtime_error("bool DCCtrl::IsParallel(): DCCtrl is not yet inizialized");
}

DCVectorBase* DCCtrl::NewVector()
{
    if(instance_)
        return instance_->NewVectorImpl();
    else
        throw std::runtime_error(" class DCCtrl::NewVector(): DCCtrl is not yet inizialized");
    
}


DCMatrixBase* DCCtrl::NewMatrix()
{
    if(instance_)
        return instance_->NewMatrixImpl();
    else
        throw std::runtime_error("DCMatrixBase* DCCtrl::NewMatrix(): DCCtrl is not yet inizialized");
    
}

void DCCtrl::VerboseOn()
{
    if(instance_)
        instance_->isVerbose_=true;
    else
        throw std::runtime_error("void DCCtrl::VerboseOn(): DCCtrl is not yet inizialized");
    
}

void DCCtrl::DebugOn()
{
    if(instance_)
    {
        instance_->isDebug_=true;
        instance_->isVerbose_=true;
    }
    else
        throw std::runtime_error("void DCCtrl::DebugOn(): DCCtrl is not yet inizialized");
    
}

void DCCtrl::VerboseOff()
{
    if(instance_)
    {
        instance_->isVerbose_=false;
        instance_->isDebug_=false;
    }
    else
        throw std::runtime_error("void DCCtrl::VerboseOff(): DCCtrl is not yet inizialized");
    
}

void DCCtrl::DebugOff()
{
    if(instance_)
        instance_->isDebug_=false;
    else
        throw std::runtime_error("void DCCtrl::DebugOff(): DCCtrl is not yet inizialized");
    
}

bool DCCtrl::GetVerbose()
{
    if (instance_)
        return instance_->isVerbose_;
    else
        throw std::runtime_error("void DCCtrl::IsVerbose(): DCCtrl is not yet initialized");
}

bool DCCtrl::GetDebug()
{
    if (instance_)
        return instance_->isDebug_;
    else
        throw std::runtime_error("void DCCtrl::IsDebug(): DCCtrl is not yet initialized");
}

void DCCtrl::SetLogFile(std::string filename)
{
    if (IsProcessZero()) { // only the first processor needs access
        if(instance_)
        {
            if(filename != "")
            {
                size_t pos = filename.find_last_of("/\\");
                if(pos != std::string::npos)
                {
                    std::string dir        = filename.substr(0, pos);
                    frizzle::filesystem::CreateDirectory(dir);
                }
                instance_->logFile_.open(filename.c_str());
            }
        }
        else
            throw std::runtime_error("void DCCtrl::SetLogFile(std::string filename): DCCtrl is not yet inizialized");
    }
}


