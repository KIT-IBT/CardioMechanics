/*! \file kaTabularizedFunction.h
   \brief Class for handling of precalculated, tabularized functions

   \author cw,fs,os, IBT - Universität Karlsruhe (TH)
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
