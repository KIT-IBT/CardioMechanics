/*
 *  ParameterMap.h
 *  CardioMechanics
 *
 *  Created by Thomas Fritz on 20.04.10.
 *  Copyright 2010 Institute for Biomedical Engineering, Karlsruhe Institute of Technology (KIT). All rights reserved.
 *
 */

#ifndef PARAMETER_MAP_H
#define PARAMETER_MAP_H

#include <map>
#include <iterator>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <set>
#include <sstream>
#include <stdexcept>
#include "stringtools.h"

class ParameterMap
{
public:
    ParameterMap(){}
    ~ParameterMap(){}
    
    std::vector<std::string> GetChildNodes(std::string parentNode);
    
    template<class T> T Get(std::string key); template<class T> std::vector<T> GetArray(std::string key);
    template<class T> T Get(std::string key, T defaultValue);
    template<class T> std::vector<T> GetArray(std::string key, std::vector<T> defaultArray);
    template<class T> void Set(std::string key, T value);
    template<class T> void SetVector(std::string key, std::vector<T> value);
    
    bool IsAvailable(std::string key);
    void Set(std::string key, std::string value);
    void ReadXml(std::string Filename);
    void WriteXml(std::string Filename);
    void Print();
    void PrintRequestedParameters();
    void PrintIgnoredParameter();
    
    static void tokenize(const std::string& str, std::vector<std::string>& tokens, const std::string& delimiters);
    
    
    std::vector<std::string> GetParameters();
    std::vector<std::string> GetRequestedParameters();
    std::vector<std::string> GetIgnoredParameters();
    
    // C++ Iterators
    
    typedef typename std::map<std::string, std::string>::iterator         iterator;
    typedef typename std::map<std::string, std::string>::const_iterator   const_iterator;
    
    iterator begin(){return(parameters_.begin());}
    const_iterator begin() const {return(parameters_.begin());}
    const const_iterator cbegin() const {return(parameters_.cbegin());}
    iterator end(){return(parameters_.end());}
    const_iterator end() const {return(parameters_.end());}
    const const_iterator cend() const {return(parameters_.cend());}
    
protected:
private:
    std::map<std::string, std::string> parameters_;
    std::set<std::string> requestedParameters_;
    std::map<std::string, std::string> requestedDefaultParameters_;
    std::string rootKey_;
    bool StripRoot();
};

template<class T>
void ParameterMap::Set(std::string key,T value)
{
    std::stringstream ss;
    
    ss << value;
    Set(key,ss.str());
}


template<class T>
T ParameterMap::Get(std::string key)
{
    std::map<std::string, std::string>::iterator it;
    
    it = parameters_.find(key);
    
    if(it == parameters_.end())
        throw std::runtime_error(key+std::string(" not found in Parameter Map. Embarrassing, but not for me !\n"));
    else
    {
        requestedParameters_.insert(key);
        return(frizzle::stringtools::StringTo<T>(it->second));
    }
}


template<class T>
std::vector<T> ParameterMap::GetArray(std::string key)
{
    std::map<std::string, std::string>::iterator it;
    
    it = parameters_.find(key);
    
    if(it == parameters_.end())
        throw std::runtime_error(key+std::string(" not found in Parameter Map. Embarrassing, but not for me !\n"));
    else
    {
        requestedParameters_.insert(key);
        
        std::vector<T> r;
        std::istringstream ss(it->second);
        
        std::string str;
        
        while(std::getline(ss,str,','))
            r.push_back(frizzle::stringtools::StringTo<T>(str));
        
        return(r);
    }
}


template<class T>
T ParameterMap::Get(std::string key, T defaultValue)
{
    std::map<std::string, std::string>::iterator it;
    
    it = parameters_.find(key);
    
    if(it == parameters_.end())
    {
        std::stringstream ss;
        ss << defaultValue;
        requestedDefaultParameters_.insert({key,ss.str() + " [DEFAULT]"});
        return(defaultValue);
    }
    else
    {
        // TODO: handle case that tag exists but no value given?
        requestedParameters_.insert(key);
        return frizzle::stringtools::StringTo<T>(it->second);
    }
}


template<class T>
std::vector<T> ParameterMap::GetArray(std::string key, std::vector<T> defaultArray)
{
    std::map<std::string, std::string>::iterator it;
    
    it = parameters_.find(key);
    
    if(it == parameters_.end())
    {
        std::stringstream ss;
        
        for(typename std::vector<T>::iterator it2 = defaultArray.begin(); it2 != defaultArray.end(); )
        {
            ss << *it2;
            it2++;
            
            if(it2 != defaultArray.end())
                ss << ",";
        }
        
        requestedDefaultParameters_.insert({key,ss.str() + " [DEFAULT]"});
        return(defaultArray);
    }
    else
    {
        requestedParameters_.insert(key);
        
        std::vector<T> r;
        std::istringstream ss(it->second);
        std::string str;
        
        while(std::getline(ss,str,','))
            r.push_back(frizzle::stringtools::StringTo<T>(str));
        return(r);
    }
}
#endif
