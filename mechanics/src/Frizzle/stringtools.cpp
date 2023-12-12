/*
 * File: stringtools.cpp
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
#include "stringtools.h"

namespace frizzle
{
namespace stringtools
{

template<>
std::string StringTo<std::string>(std::string str)
{
    return str;
}

template<>
const char* StringTo<const char*>(std::string str)
{
    return str.c_str();
}

template<>
int StringTo<int>(std::string str)
{
    return std::stoi(str);
}


template<>
long StringTo<long>(std::string str)
{
    return std::stol(str);
}

template<>
long long StringTo<long long>(std::string str)
{
    return std::stoll(str);
}

template<>
unsigned long StringTo<unsigned long>(std::string str)
{
    return std::stoul(str);
}

template<>
unsigned long long StringTo<unsigned long long>(std::string str)
{
    return std::stoull(str);
}

template<>
float StringTo<float>(std::string str)
{
    return std::stof(str);
}

template<>
double StringTo<double>(std::string str)
{
    return std::stod(str);
}

template<>
long double StringTo<long double>(std::string str)
{
    return std::stold(str);
}

template<>
bool StringTo<bool>(std::string str)
{
    if(str == "true" || str == "True" || str == "TRUE" || str == "1" || str == "yes" || str == "YES" || str == "Yes")
        return true;
    else if(str == "false" || str == "False" || str == "FALSE" || str == "0" || str == "no" || str == "NO" || str == "No")
        return false;
    else
    {
        std::cout << "Hall";
        throw std::runtime_error(std::string("bool frizzle::stringtools::StringTo<bool>(const string& str): Can't convert std::string " + str + " to boolean "));
    }
}

}
}
