/*
 * File: ParameterMap.cpp
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


#include "ParameterMap.h"

#include <algorithm>


static int GetCurrentLineNumber(std::ifstream& i)
{
    int a = i.tellg();
    int b = 0;
    int n = 0;
    
    std::string dummy;
    i.seekg(0);
    
    while(b < a)
    {
        std::getline(i,dummy);
        b = i.tellg();
        n++;
    }
    
    i.seekg(a);
    return n;
}


void ParameterMap::tokenize(const std::string& str, std::vector<std::string>& tokens, const std::string& delimiters = ".")
{
    std::string::size_type lastPos = str.find_first_not_of(delimiters, 0);
    std::string::size_type pos     = str.find_first_of(delimiters, lastPos);
    
    while(std::string::npos != pos || std::string::npos != lastPos)
    {
        tokens.push_back(str.substr(lastPos, pos - lastPos));
        lastPos = str.find_first_not_of(delimiters, pos);
        pos     = str.find_first_of(delimiters, lastPos);
    }
}

std::string tabs(int n)
{
    std::string str;
    for(int i = 0; i < n; i++)
        str += std::string("\t");
    return(str);
}

void ParameterMap::Set(std::string key, std::string value)
{
    std::map<std::string, std::string>::iterator it;
    it = parameters_.find(key);     // If key already exists -> delete element!
    if(it != parameters_.end())
        parameters_.erase(key);
    parameters_.insert(std::pair<std::string, std::string>(key, value));
}


// added by Steffen Schuler
std::vector<std::string> ParameterMap::GetChildNodes(std::string parentNode)
{
    std::vector<std::string> parentTokens;
    tokenize(parentNode, parentTokens);
    
    std::vector<std::string> childNodes;
    
    for(auto it = parameters_.begin(); it != parameters_.end(); it++)
    {
        if(it->first.find(parentNode + ".") != std::string::npos)
        {
            std::vector<std::string> childTokens;
            tokenize(it->first, childTokens);
            
            std::string childNode = parentNode + "." + childTokens.at(parentTokens.size());
            
            if(std::find(childNodes.begin(), childNodes.end(), childNode) == childNodes.end())
                childNodes.push_back(childNode);
        }
    }
    return(childNodes);
}


void ParameterMap::Print()
{
    std::cout << "Parameters: " << std::endl;
    for(std::map<std::string, std::string>::iterator it = parameters_.begin(); it != parameters_.end(); it++)
        std::cout << it->first << ":" << it->second << "\n";
    std::cout << std::endl;
}

void ParameterMap::PrintRequestedParameters()
{
    std::cout << "Requested parameters:\n";
    for(auto it = requestedParameters_.begin(); it != requestedParameters_.end(); it++)
    {
        std::cout << *it << ": " << parameters_[*it] << "\n";
    }
    std::cout << std::endl;
}


void ParameterMap::PrintIgnoredParameter()
{
    std::cout << "Ignored parameters: \n";
    for(std::map<std::string, std::string>::iterator it = parameters_.begin(); it != parameters_.end(); it++)
        if(requestedParameters_.find(it->first) == requestedParameters_.end())
            std::cout << it->first << ":" << it->second;
}

std::vector<std::string> ParameterMap::GetParameters()
{
    std::vector<std::string> s;
    for(auto it = parameters_.begin(); it != parameters_.end(); it++)
    {
        std::stringstream ss;
        ss << it->first << ": " << it->second;
        s.push_back(ss.str());
    }
    return s;
}

std::vector<std::string> ParameterMap::GetIgnoredParameters()
{
    std::vector<std::string> s;
    for(auto it = parameters_.begin(); it != parameters_.end(); it++)
    {
        std::stringstream ss;
        if(requestedParameters_.find(it->first) == requestedParameters_.end())
        {
            ss << it->first << ": " << it->second;
            s.push_back(ss.str());
        }
    }
    return s;
}

std::vector<std::string> ParameterMap::GetRequestedParameters()
{
    std::vector<std::string> s;
    for(auto it = parameters_.begin(); it != parameters_.end(); it++)
    {
        std::stringstream ss;
        if(requestedParameters_.find(it->first) != requestedParameters_.end())
        {
            ss << it->first << ": " << it->second;
            s.push_back(ss.str());
        }
    }
    for(auto i : requestedDefaultParameters_)
    {
        std::stringstream ss;
        ss << i.first << ": "  << i.second;
        s.push_back(ss.str());
    }
    sort(s.begin(),s.end());
    return s;
}

void ParameterMap::WriteXml(std::string Filename)
{
    std::ofstream                      file(Filename.c_str());
    file << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>";
    std::map<std::string, std::string>::iterator it;
    
    std::vector<std::string>                tags;
    for(it = parameters_.begin(); it != parameters_.end(); it++)
    {
        std::string         str   = this->rootKey_ + it->first;
        std::string         value = it->second;
        std::vector<std::string> tokens;
        tokenize(str, tokens);
        int            i;
        for(i = 0; i < tags.size(); i++)
        {
            if(tags.at(i) != tokens.at(i))
            {
                for(int j = tags.size() - 1; j >= i; j--)
                {
                    file << "</" << tags.at(j) << ">";
                    if(j != i)
                        file << "\n" << tabs(j-1);
                }
                tags.resize(i);
                break;
            }
        }
        for(int j = i; j < tokens.size(); j++)
        {
            tags.push_back(tokens.at(j));
            
            file << "\n" << tabs(j) << "<" << tokens.at(j) << ">";
        }
        file << value;
    }
    for(int j = tags.size() - 1; j >= 0; j--)
    {
        file << "</" << tags.at(j) << ">\n" << tabs(j-1);
    }
    file.close();
}

void ParameterMap::ReadXml(std::string Filename)
{
    std::vector<std::string> tags;
    std::ifstream       file(Filename.c_str());
    if(!file.good())
        throw std::runtime_error("File " + Filename + " does not exist, while i'm an existentialist");
    int            length;
    file.seekg(0, std::ios::end);
    length = file.tellg();
    file.seekg(0, std::ios::beg);
    
    char c;
    bool isOpenKey = false;
    bool isEndTag  = false;
    int lineNumber = 0;
    if(file.good())
    {
        while(!file.eof())
        {
            file >> c;
            if(c == char('<'))           //Read Tag
            {
                std::string str;
                file >> c;
                lineNumber = GetCurrentLineNumber(file);
                
                if(c=='?') {
                    // ignore xml declaration, which does NOT have closing tags
                    char d = 0;
                    file >> c >> d;
                    while( c!='>' && !(c=='?'&&d=='>') )
                    {
                        c=d;
                        file >> d;
                        if(file.eof())
                            throw std::runtime_error(std::string("Parameter file: ") + std::string(Filename) + std::string(" - Error in line ") + std::to_string(lineNumber) + ": " + std::string("XML Declaration has to to be closed with ?> ."));
                    }
                    continue;
                }
                
                if(c == ('!'))
                {
                    char d = 0;
                    char e = 0;
                    file >> c >> d;
                    if(c != '-' || d != '-')
                        throw std::runtime_error(std::string("Parameter file: ") + std::string(Filename) + std::string(" - Error in line ") + std::to_string(lineNumber) + ": " + std::string("Comment has to be opened with <!-- .")); 
                    d = e = 0;
                    while(e != '>' || ( e == '>' && (c!='-' || d!='-')))
                    {
                        c=d;
                        d=e;
                        file >> e;
                        if(file.eof())
                            throw std::runtime_error(std::string("Parameter file: ") + std::string(Filename) + std::string(" - Error in line ") + std::to_string(lineNumber) + ": " + std::string("Comment has to to be closed with --> .")); 
                    }
                    continue;
                }
                if(c == char('/') || c == char('\\'))
                {
                    isEndTag = true;
                    file >> c;
                }
                while(c != char('>'))
                {
                    if(file.eof())
                        throw std::runtime_error(std::string("Parameter file: ") + std::string(Filename) + std::string(" - Error in line ") + std::to_string(lineNumber) + ": "+ std::string("Missing [>] behind ") + str + std::string(". Honey, you'd better fix it!"));
                    str += c;
                    file >> c;
                }
                if(isEndTag)
                {
                    if(tags.size() == 0)
                    {
                        throw std::runtime_error(std::string("Parameter file: ") + std::string(Filename) + std::string(" - Error in line ") + std::to_string(lineNumber) + ": "+ std::string("There is no tag to close with [ </") + str + std::string("> ]. Darling, Fix it! NOW !!"));
                    }
                    if(str == tags.back())
                    {
                        tags.pop_back();
                        if(isOpenKey)
                            isOpenKey = false;
                    }
                    else
                    {
                        
                        throw std::runtime_error(std::string("Parameter file: ") + std::string(Filename) + std::string(" - Error in line ") + std::to_string(lineNumber) + ": "+ std::string("Tag [ <") + tags.back() + std::string("> ] is closed with [ <\\") + str + std::string("> ]. What's wrong with you ?"));
                    }
                    isEndTag = false;
                }
                else
                {
                    if(!isOpenKey)
                        tags.push_back(str);
                    else
                        throw std::runtime_error(std::string("Parameter file: ") + std::string(Filename) + std::string(" - Error in line ") + std::to_string(lineNumber) + ": "+ std::string("Missing [ </") + tags.back() + std::string("> ]. Try your best to fix it!"));
                }
            }
            else             //Read Value
            {
                std::string value;
                std::string key;
                bool   skipWhiteSpace = true;
                while(c != char('<') && !file.eof())
                {
                    value += c;
                    if(skipWhiteSpace)
                        file >> c;
                    else
                        file >> std::noskipws >> c >> std::skipws;
                    skipWhiteSpace = false;
                }
                if(tags.size() == 0)
                {
                    if(!file.eof() || value.size() != 0)
                        throw std::runtime_error(std::string("Parameter file: ") + std::string(Filename) + std::string(" - Error in line ") + std::to_string(lineNumber) + ": "+ std::string("Value ") + value + std::string(" is not embraced. Embrace it with your embracing attitude"));
                    else
                        break;
                }
                else
                {
                    file.unget();
                    key = tags.at(0);
                    for(int i = 1; i < tags.size(); i++)
                        key += std::string(".")+tags.at(i);
                    if(parameters_.find(key) != parameters_.end())
                        throw std::runtime_error(std::string("Parameter file: ") + std::string(Filename) + std::string(" - Error in line ") + std::to_string(lineNumber) + ": "+ std::string("Key ") + key + std::string(" does already exist in the parameters file ") + Filename + std::string(". Keys love to be unique!"));
                    
                    else
                    {
                        Set(key, value);
                        isOpenKey = true;
                    }
                }
            }
        }
        file.close();
        this->StripRoot();
    }
    else
        throw std::runtime_error(std::string("Error in line ") + std::to_string(lineNumber)+ std::string(".File ") + Filename + std::string(" does not exist !"));
}

bool ParameterMap::IsAvailable(std::string key)
{
    std::map<std::string, std::string>::iterator it;
    it = parameters_.find(key);
    if(it == parameters_.end())
        return(false);
    else
        return(true);
}

bool ParameterMap::StripRoot()
{
    // test for common root entry
    bool hasCommonRoot = true;
    std::string item      = "";
    std::string item_prev = "";
    for (auto it=parameters_.begin(); it!=++parameters_.begin(); it++) {
        std::string key = it->first;
        std::stringstream ss(key);
        std::getline(ss, item_prev, '.');
    }
    for (auto it=++parameters_.begin(); it!=parameters_.end(); it++) {
        std::string key = it->first;
        std::stringstream ss(key);
        std::getline(ss, item, '.');
        if (item!=item_prev) {
            hasCommonRoot = false;
            break;
        }
        item_prev = item;
    }
    
    // extract root element string to a separate variable
    if (hasCommonRoot) {
        std::string key = parameters_.begin()->first;
        std::stringstream ss(key);
        std::string item = "";
        std::getline(ss, item, '.');
        this->rootKey_ = item + ".";
        
        std::map<std::string, std::string> tmpMap;
        for (auto it=parameters_.begin(); it!=parameters_.end(); it++) {
            std::string key = it->first;
            key.erase(0,this->rootKey_.length()); // remove root string
            tmpMap[key] = it->second;
        }
        parameters_ = tmpMap;
    }
    
    return hasCommonRoot;
}

template<>
void ParameterMap::Set<bool>(std::string key,bool value)
{
    if(value == true)
        Set(key,"true");
    else
        Set(key,"false");
}
