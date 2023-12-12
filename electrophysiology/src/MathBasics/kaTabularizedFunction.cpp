/*
 * File: kaTabularizedFunction.cpp
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


#include <kaTabularizedFunction.h>

kaTabularizedFunction<double, 255> SinTable(sin);
kaTabularizedFunction<double, 255> CosTable(cos);
