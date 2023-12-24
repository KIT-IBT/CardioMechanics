/*
 * File: ProgressBar.h
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


#ifndef PROGRESSBAR_H
#define PROGRESSBAR_H

#include "Verbosity.h"
#include <cstring>

using namespace std;

class ProgressBar {
 public:
  ProgressBar(string text, long pmax, unsigned int plevel = 0) : level(plevel) {
    int len     = 50;
    float steps = 4.0F*len;

    max = pmax - 1;
    if (steps > pmax) {
      steps = pmax;
    }
    step   = pmax / steps;
    progi  = 100.0F / steps;
    next   = 0;
    prog   = 0;
    slashc = 0;
    strcpy(slash, "-\\|/#");

    pre = text.substr(0, 20);
    pre.append(20-text.length(), ' ');
    pre.append(" [");
    post.assign(len, ' ');
    post.append("] ");

    Verbosity::Stream(level) << '\r' << pre << slash[slashc] << post << setw(3) << (int)prog << '%' << flush;
  }

  ~ProgressBar() {
    if (slashc != 4) {
      prog   = 100;
      slashc = 4;
      set(max);
    }
  }

  inline void set(long cur) {
    if (cur < next)
      return;

    Verbosity::Stream(level) << '\r' << pre << slash[slashc] << post << setw(3) << (int)prog << '%' << flush;
    if (slashc == 4) {Verbosity::Stream(level) << endl; return;}
    prog += progi;

    if (++slashc > 3) {
      pre.push_back(slash[4]);
      slashc = 0;
      post   = post.substr(1);
    }

    next += step;
    if (next > max) {
      next   = max;
      prog   = 100;
      slashc = 4;
    }
  }

 private:
  string pre;
  string post;
  float next;
  long max;
  float step;
  float prog;
  float progi;
  unsigned int level;
  char slash[6];
  int slashc;
};  // class ProgressBar

#endif  // ifndef PROGRESSBAR_H
