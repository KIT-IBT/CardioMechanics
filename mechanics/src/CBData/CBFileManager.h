/*
 *  CBFileManager.h
 *  CardioMechanics
 *
 *  Created by Lukas Baron on Thu Apr 27 2017.
 *  Copyright 2017 Institute for Biomedical Engineering, Karlsruhe Institute of Technology (KIT). All rights reserved.
 *
 */

#pragma once

#include "CBTiming.h"
#include "ParameterMap.h"
#include "CBDataFromFile.h"
#include "CBElementSolid.h"

/// manage multiple tension file objects of type CBDataFromFile, and takes care
/// that each file is opened just once.
/// There should exist only one opened file object per material, / the file
/// object gets accessed via Get(element), it identifies the needed tension file
/// directly from a passed element instance / elements material

class CBFileManager {
protected:
    ParameterMap* parameters_;
    CBTiming* timing_;
    std::map<TInt, CBDataFromFile*> tensionFiles_;
    std::vector<std::pair<int, double>> activationTimes_;
    float latOffset_;
    bool processLAT_;
public:
    CBFileManager() { }
    ~CBFileManager();
    
    /// initialize with passed object instances, and create one DataFromFile object for each materials that needs tension from file
    void Init(ParameterMap* parameters, CBTiming* timing);
    
    /// returns the corresponding DataFromFile object for a certain element (detected by it's material index)
    CBDataFromFile* GetDataFromFileObject(const CBElementSolid* ele) const;
    
    /// returns a scalar value corresponding to one element, similar to the syntax of CBData's Get() function
    /// better rename to GetDatum()?
    TFloat Get(TFloat time, CBElementSolid* ele) const;
    
    /// Returns the activationTimes_ vector to access it
    std::vector<std::pair<int, double>> GetActivationTimesVector() {
        return activationTimes_;
    }
    
    /// Returns the latOffset_ variable to access it
    float GetLatOffset() {
        return latOffset_;
    }
    
    /// read the cell activation times for each cell from a file
    int ReadActivationTimesFromFile(std::string filePath);
    
    /// Helper-method to help split strings
    std::vector<std::string> SplitString(std::string input);
    
    /// Return the processLAT value. If it's true a LAT Map was given and should be processed
    bool GetProcessLAT() {
        return processLAT_;
    }
};

