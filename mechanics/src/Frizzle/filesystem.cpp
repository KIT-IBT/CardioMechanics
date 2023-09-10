#include "filesystem.h"
#include<vector>
#include<iostream>
#include<stdexcept>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <string>

namespace frizzle
{
namespace filesystem
{
bool CreateDirectory(const std::string& path)
{
    struct stat st;
    
    if (stat(path.c_str(), &st) == 0 && S_ISDIR(st.st_mode))
    {
        return true;
    }
    else if (stat(path.c_str(), &st) == 0 && !S_ISDIR(st.st_mode))
    {
        throw std::runtime_error("void mkdir_recursively(const std::string& path) " + path + " is not a directory");
    }
    else
    {
        bool root = false;
        int i=0;
        
        while(path.at(i)==' ')
            i++;
        
        if(path.at(i) == '/' || path.at(i) == '\\')
        {
            root=true;
            i++;
        }
        
        std::string str = path.substr(i,std::string::npos);
        std::vector<std::string>dir;
        size_t pos=0;
        do
        {
            pos = str.find_first_of("/\\");
            if(str != "/" && str != "\\" && str != "")
                dir.push_back(str.substr(0,pos));
            
            if(pos != std::string::npos)
                str = str.substr(pos+1,std::string::npos);
        }
        while(pos != std::string::npos);
        
        std::string newPath = "";
        if(root)
            newPath = "/";
        
        for(auto i: dir)
        {
            newPath += i + "/";
            if(mkdir( newPath.c_str(), 0755) != 0 && errno != EEXIST)
                return false;
        }
        return true;
    }
}

std::string FullPath(const std::string& Filename)
{
    // source: http://www.manpagez.com/man/3/realpath/
    const char* file_name = Filename.c_str();
    const char* full_path = realpath(file_name, NULL);
    return std::string(full_path);
}

std::string ParentPath(const std::string& path)
{
    // code copied from 'CreateDirectory', except push_back and mkdir
    bool root = false;
    int i=0;
    
    while(path.at(i)==' ')
        i++;
    
    if(path.at(i) == '/' || path.at(i) == '\\')
    {
        root=true;
        i++;
    }
    std::string str = path.substr(i,std::string::npos);
    std::vector<std::string>dir;
    size_t pos=0;
    do
    {
        pos = str.find_first_of("/\\");
        if(str != "/" && str != "\\" && str != "" && pos != std::string::npos)
            dir.push_back(str.substr(0,pos));
        
        if(pos != std::string::npos)
            str = str.substr(pos+1,std::string::npos);
    }
    while(pos != std::string::npos);
    
    std::string newPath = "";
    if(root)
        newPath = "/";
    
    for(auto i: dir)
    {
        newPath += i + "/";
    }
    return newPath;
}

}
}
