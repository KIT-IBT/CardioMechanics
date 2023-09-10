/*
 *  CBFileManager.cpp
 *  CardioMechanics
 *
 *  Created by Lukas Baron on Thu Apr 28 2017.
 *  Copyright 2017 Institute for Biomedical Engineering, Karlsruhe Institute of Technology (KIT). All rights reserved.
 *
 */

#include "CBFileManager.h"
#include "math.h"
#include <limits>

CBFileManager::~CBFileManager() {
    for (auto m : tensionFiles_)
        delete m.second;
    tensionFiles_.clear();
}


void CBFileManager::Init(ParameterMap* parameters, CBTiming* timing) {
    parameters_ = parameters;
    timing_ = timing;
    
    // look in parameterMap for materials that need tension from files, and load/initialize the respective files
    std::vector<std::string> matStrings = parameters->GetChildNodes("Materials");
    for (auto &mat : matStrings) {
        if (parameters->Get<std::string>(mat + ".TensionModel", "None") == "File") {
            TInt matIndex = -1; // TODO: find correct material, this could be somewhere in the elements where they get their GetMaterialIndex() entry, or where SetMaterialIndex() gets called, or in CBConstitutiveModel
            //TODO: Is this working?
            std::vector<std::string> tokens;
            ParameterMap::tokenize(mat, tokens, "_");
            assert(tokens.size()==2);
            // TODO: this does not work in the "Default" case, i.e. when no specific material number is given
            if (tokens.at(1) != "Default") {
                matIndex = std::stoi(tokens.at(1));
            }
            else {
                std::cout << "Warning from CBFileManager::Init : Default material found, do not know what to set as material index for the tension model." << std::endl;
                throw std::runtime_error("CBFileManager::Init : TensionModel in Mat_Default is currently not supported.");
            }
            CBDataFromFile* file = new CBDataFromFile(parameters, mat + ".TensionFromFile");
            tensionFiles_.insert(std::pair<int, CBDataFromFile*>(matIndex, file));
            // insert this file for all other materials, unless different tension model is set?
        }
    }
    
    //----------------------- Local Activation Time Reader -----------------------
    
    //read preloading offset time from parameterMap. This offset is later added to the activation time to ensure proper function
    //Set the offset to zero in case it isn't given in the xml, otherwise set it to the value written in the xml
    
    latOffset_ = parameters->Get<TFloat>("Mesh.LAT.Offset", 0.0);
    
    //get path to the LAT(local activation times)-file in parameterMap
    std::string filePathToLAT;
    int retCode = 0;
    processLAT_ = false;
    
    try {
        filePathToLAT = parameters->Get<std::string>("Mesh.LAT.FilePath", "null");
        //if it was successfull process the given file
        if (filePathToLAT.compare("null") != 0) {
            //compare returns 0 on equality, we want anything but 0 here
            processLAT_ = true;
            
            //read the activation times from the given file path and process them
            retCode = ReadActivationTimesFromFile(filePathToLAT);
        } else {
            //if there was no lat map given, we end up in this section here. Nothing to do
        }
    } catch (const std::runtime_error& e) {
        //no lat map was given
        //nothing to do
    }
    
    //if there was no lat map given, this section gets skipped
    if (processLAT_ && retCode != 1) {
        //handle error regarding parsing the activation times from file
        throw std::runtime_error("CBFileManager::Init : Reading the activation times from the file \"" + filePathToLAT + "\" was unsuccessful.");
    } else if (processLAT_ && (activationTimes_.size() == 0)) {
        std::cout << "There were no activation times successfully added." << std::endl;
    }
}


CBDataFromFile* CBFileManager::GetDataFromFileObject(const CBElementSolid* ele) const {
    return tensionFiles_.at(ele->GetMaterialIndex());
}


TFloat CBFileManager::Get(TFloat time, CBElementSolid* ele) const {
    TInt index = ele->GetIndex();
    CBDataFromFile* f = GetDataFromFileObject(ele);
    return f->Get(time, index);
}

/// Method to read the activation times for each cell from a (text) file.
/// This Method should store the activation times in the FileManager,
/// so you have the necessary data centralized
/// WARNING: This Method is higly unsafe to multi-threaded usage. Since
/// you only need to read activation times once, this should be no
/// problem if it's done by only one thread.
/* @author Armin Mueller
 * @param filePath
 *		The path to the file (.txt, .csv) in which the activation times are stored.
 *		(The format of the data in this file should be: id activationTimeInSeconds. For example:
 *		0 0.000
 *		1 0.005
 *		2 0.010)
 * @return int
 *		The methods returning integer symbolizes the success or failure of the method.
 *		On success it returns 1, on failure it returns 0.
 */
int CBFileManager::ReadActivationTimesFromFile(std::string filePath) {
    std::string line;
    std::vector<std::string> temp;
    std::pair<int, double> retValue;
    std::string actTime;
    int lineCtr = 0;
    bool reachedEOF = false;
    //activationTimes_ is defined globally
    
    std::ifstream ifs;
    ifs.open(filePath.c_str(), std::ifstream::in);	//open file to read
    
    //check if file is ok
    if (!ifs.good()) {
        throw std::runtime_error("CBFileManager::Init : Failed to read LAT-file.");
        //method failed, close the file, return with error code 0
        ifs.close();
        return 0;
    } else {
        while (ifs.good()) {
            getline(ifs, line);	//read next line from file, save in: line
            lineCtr++;
            
            if (ifs.eof() && line.empty()) {
                //reached end of file, which is a empty line. Return with success code 1
                ifs.close();
                return 1;
            } else if (ifs.eof() && !line.empty()) {
                //reached end of file in last line of document. Set reachedEOF to true,
                //so this last line gets processed and after that close ifs
                reachedEOF = true;
            }
            
            if (lineCtr == 1) {
                //the first line of the file contains header information. So skip this line or
                //maybe use the information for error prevention (e.g. too few/much data points given, ...)
                continue;
            }
            
            //process new dataset
            temp = SplitString(line);
            
            if (temp.size() != 2) {
                //if it doesn't hold 2 elements, something went wrong. Return with error code 0.
                ifs.close();
                return 0;
            }
            
            retValue.first = std::stoi(temp.at(0));		//set the id of the activation time pair
            
            actTime = temp.at(1);
            
            //handle "NaN"-values for activation times
            if (actTime.compare("NaN") == 0) {
                //compare returns 0 for equality
                retValue.second = 0.0;
            } else {
                retValue.second = std::stod(actTime);
            }
            
            activationTimes_.push_back(retValue);	//add the new dataset
            
            //check whether it was the last line and it contained also the EOF
            if (reachedEOF) {
                //EOF and last data item were in the same line. So close ifs and return 1 for success.
                ifs.close();
                return 1;
            }
        }//end while
    }//end if
    
    //if the method reaches this point something went wrong, exit with error code 0
    ifs.close();
    return 0;
}

/// Method to split an string at whitespaces.
/* @author Armin Mueller
 * @param input
 *		The string which should be splitted at the whitespaces
 * @return std::vector<std::string>
 *		The method returns an vector which holds in each line
 *		one word of the input string
 */
std::vector<std::string> CBFileManager::SplitString(std::string input) {
    std::string buffer;
    std::stringstream ss(input); //insert the string into a stream
    std::vector<std::string> tokens; //create vector to hold the words
    
    while (ss >> buffer) {
        tokens.push_back(buffer);
    }
    
    return tokens;
}
