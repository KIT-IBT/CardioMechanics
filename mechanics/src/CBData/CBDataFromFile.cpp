/*
 *  CBDataFromFile.cpp
 *  CardioMechanics
 *
 *  Created by Thomas Fritz on 14.09.11.
 *  Copyright 2011 Institute for Biomedical Engineering, Karlsruhe Institute of Technology (KIT). All rights reserved.
 *
 */

#include "CBDataFromFile.h"
#include <algorithm> // needed for std::min()
#include <stdlib.h> // needed for realpath()
#include "filesystem.h"


CBDataFromFile::CBDataFromFile() {
    dataBegin_.clear();
    dataEnd_.clear();
    timeDataBegin_ = 0;
    timeDataEnd_   = 0;
    isDataSetLoaded_ = false;
}

CBDataFromFile::CBDataFromFile(ParameterMap *parameters, std::string parameterKey) {
    CBDataFromFile();
    
    startTime_ = parameters->Get<TFloat>(parameterKey + ".StartTime", 0.0);
    period_ = parameters->Get<TFloat>(parameterKey + ".Period", 1.0);
    
    /// File at <Filename> is structured as followed and later stored in fileList_:
    /// # TIME  # PATH to file containing the tension values for each element at TIME
    ///
    /// The individual files at PATH containing the element tension have to be binary with the following structure
    /// int32     nrElements
    /// double  tension[0]
    /// double  tension[1]
    /// double  tension[...]
    /// double  tension[n]
    std::string filename = parameters->Get<std::string>(parameterKey + ".Filename");
    
    Init(filename);
}

void CBDataFromFile::Init(std::string filename) {
    LoadFileList(filename);
    LoadDataSet(0);
}

void CBDataFromFile::LoadFileList(std::string filename) {
    if (sizeof(double) != 8)
        throw std::runtime_error("UUUggh, At my days, sizeof(double) used to be 8");
    
    std::ifstream file(filename.c_str());
    if (!file.good()) {
        throw std::runtime_error(
                                 "CBData::LoadFileList(std::string filename): File list: " + filename + " does not exist !!!");
    }
    
    std::string::size_type pos = filename.rfind('/');
    std::string path = filename.substr(0, pos+1);
    
    std::string str;
    while (getline(file, str)) {
        if (str.empty())
            continue;
        
        std::stringstream ss(str);
        
        std::size_t fc = str.find_first_not_of(" \n\t");
        if (str.at(fc) == '#')
            continue;
        std::pair<TFloat, std::string> dataFile;
        
        if (ss.good()) {
            ss >> dataFile.first;
            
            // use the full path information, since tension.list mostly contains paths relative to itself
            std::string name;
            ss >> name;
            dataFile.second = name;
        } else {
            throw std::runtime_error("CBData::LoadFileList(std::string filename): File list: " + filename + " is corrupt.");
        }
        fileList_.push_back(dataFile);
    }
    file.close();
    
    // Check if file list is correct
    for (std::vector<std::pair<TFloat, std::string>>::iterator it = fileList_.begin(); it != fileList_.end(); it++) {
        file.open(it->second.c_str());
        if (!file.good()) {
            throw std::runtime_error(
                                     "CBData::LoadFileList(std::string filename): File list: " + filename + " is corrupt. File " + it->second +
                                     " does not exist !!!");
        }
        file.close();
    }
    isDataSetLoaded_ = true;
    timeDataBegin_ = 0;
    timeDataEnd_   = 0;
} // CBDataFromFile::LoadFileList

bool CBDataFromFile::LoadDataSet(TFloat time) {
    int32_t numValues;
    
    // change interval bounds and reload corresponding data only if time is not
    // within the currently loaded interval
    if ((time < timeDataBegin_) || (time >= timeDataEnd_) || (dataBegin_.size() == 0) || (dataEnd_.size() == 0)) {
        std::vector<std::pair<TFloat, std::string>>::iterator it = fileList_.begin();
        while (time >= it->first && it != fileList_.end()) {
            it++;
        }
        
        if (it != fileList_.end()) {
            std::ifstream file;
            
            file.open(it->second.c_str(), std::ios::in | std::ios::binary);
            file.read((char *)(&numValues), sizeof(int32_t));
            dataEnd_.clear();
            dataEnd_.reserve(numValues);
            double *values = new double[numValues];
            file.read((char *)values, sizeof(double) * numValues);
            if (file.fail())
                throw std::runtime_error("CBDataFromFile::LoadDataSet(): File " + it->second + " is corrupt !");
            for (int i = 0; i < numValues; i++)
                dataEnd_.push_back(values[i]);
            file.close();
            timeDataEnd_ = it->first;
            
            if (it != fileList_.begin()) // time < first element list -> load first element time twice
                it--;
            int32_t dummy;
            
            file.open(it->second.c_str(), std::ios::in | std::ios::binary);
            file.read((char *)(&dummy), sizeof(int32_t));
            if (dummy != numValues)
                throw std::runtime_error("CBDataFromFile::LoadDataSet(): Files don't fit to each other !");
            dataBegin_.clear();
            dataBegin_.reserve(numValues);
            file.read((char *)values, sizeof(double) * numValues);
            if (file.fail())
                throw std::runtime_error("CBDataFromFile::LoadDataSet(): File " + it->second + " is corrupt !");
            for (int i = 0; i < numValues; i++) {
                dataBegin_.push_back(values[i]);
            }
            delete[] values;
            file.close();
            timeDataBegin_ = it->first;
            if (it != fileList_.begin()) // time < first time in list -> set interval bounds to avoid reloading
                timeDataBegin_ = std::min(time, it->first);
        } else {
            dataBegin_.clear();
            dataEnd_.clear();
        }
    }
    return true;
} // CBDataFromFile::LoadDataSet

TFloat CBDataFromFile::Get(TFloat time, TInt index) {
    TFloat t = time - startTime_;
    
    if ((t < 0) || (index >= dataBegin_.size())) {
        return 0.0;
    } else {
        t = fmod(t, period_);
        
        if ((t <= timeDataBegin_) || (t > timeDataEnd_) )
            LoadDataSet(t);
        if ((dataEnd_.size() == 0) || (dataBegin_.size() == 0)) {
            return 0;
        } else if (timeDataBegin_ == timeDataEnd_) {
            return dataBegin_[index];
        }
        
        /// do linear interpolation for time t if not in fileList_
        return dataBegin_[index] + (dataEnd_[index] - dataBegin_[index]) / (timeDataEnd_ - timeDataBegin_) *
        (t - timeDataBegin_);
    }
}

std::vector<TFloat> CBDataFromFile::GetTimePoints() {
    std::vector<TFloat> timePoints;
    
    for (std::vector<std::pair<TFloat, std::string>>::iterator it = fileList_.begin(); it != fileList_.end(); it++) {
        timePoints.push_back(it->first);
    }
    return timePoints;
}
