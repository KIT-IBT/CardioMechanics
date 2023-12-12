/*
 * File: stringtools.h
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

#include<string>
#include<iostream>
#include<stdexcept>

namespace frizzle
{
namespace stringtools
{
template<class T>
T StringTo(std::string str)
{
    throw std::runtime_error(std::string("T frizzle::StringTo(const string& str): Cannot cast to requested type"));
}
template<>
std::string StringTo<std::string>(std::string str);
template<>
const char* StringTo<const char*>(std::string str);
template<>
int StringTo<int>(std::string str);
template<>
long StringTo<long>(std::string str);
template<>
long long StringTo<long long>(std::string str);
template<>
unsigned long StringTo<unsigned long>(std::string str);
template<>
unsigned long long StringTo<unsigned long long>(std::string str);
template<>
float StringTo<float>(std::string str);
template<>
double StringTo<double>(std::string str);
template<>
long double StringTo<long double>(std::string str);
template<>
bool StringTo<bool>(std::string str);
}
}
