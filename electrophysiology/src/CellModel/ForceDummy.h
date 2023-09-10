/**@file ForceDummy.h
 * @brief <Brief (one-line) description here.>
 *
 * Please see the wiki for details on how to fill out this header:
 * https://intern.ibt.uni-karlsruhe.de/wiki/Document_a_IBT_C%2B%2B_tool
 *
 * @version 1.0.0
 *
 * @date Created <Your Name> (yyyy-mm-dd)
 *
 * @author Your Name\n
 *         Institute of Biomedical Engineering\n
 *         Karlsruhe Institute of Technology (KIT)\n
 *         http://www.ibt.kit.edu\n
 *         Copyright yyyy - All rights reserved.
 *
 * @see ...
 */

#ifndef FORCE_DUMMY_H
#define FORCE_DUMMY_H

#include <ForceModelBasis.h>

template<class T>
class ForceDummy : virtual public vbForceModel<T> {
 public:
  ForceDummy() {}

  ~ForceDummy() {}

  virtual void Init() {}

  virtual T Calc(double tinc, T stretch, T velocity, T &Ca, int euler = 1) {return 0.0;}

  virtual void Print(ostream &tempstr) {}

  virtual int GetSize(void) {return 0;}

  virtual T *GetBase(void) {return NULL;}

  virtual int GetNumParameters() {return 0;}

  virtual void WriteForceStatus(FILE *writeto) {}

  virtual void ReadForceStatus(FILE *readfrom) {}

  virtual void GetParameterNames(vector<string> &getpara) {}
};
#endif  // ifndef FORCE_DUMMY_H
