/*
 * File: kaTabularizedFunction.h
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


#ifndef KATABULARIZEDFUNCTION_H
#define KATABULARIZEDFUNCTION_H

#include <kaMachineOS.h>

template<class T, int count>
class kaTabularizedFunction {
 public:
  T val[count];
  kaTabularizedFunction(T(*func)(T)) {
    for (int i = 0; i < count; i++)
      val[i] = func((T)i*(T)M_PI/(T)count);
  }
};


extern kaTabularizedFunction<double, 255> SinTable;
extern kaTabularizedFunction<double, 255> CosTable;


#endif  // ifndef KATABULARIZEDFUNCTION_H
