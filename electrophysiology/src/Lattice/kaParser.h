/*
 * File: kaParser.h
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

#ifndef KAPARSER_H
#define KAPARSER_H

#include <kaMachineOS.h>

static const int ARG_SIZE = 256;
static const int BUF_NUM  = 256;

static const char *errorstr[] = {
  "no error",
  "unknown option",
  "incomplete command line",
  "no such option selected",
  "maximum number of options reached"
};

typedef enum {
  ERROR_NONE = 0, UNKNOWN_OPTION = 1, INCOMPLETE_COMMAND_LINE = 2, NO_SUCH_OPTION = 3, NO_FURTHER_OPTION = 4
} CLerr;


//! Class for handling of command line arguments.
/*!
   kaParser copies the argments and allows searching by keywords.

   \bugs Invalid arguments are not reported.
 */

class kaParser {
 public:
  kaParser(void);
  kaParser(int c, char **v);
  ~kaParser(void);
  int   HowMany(string opt);
  char *GetOption(string opt, int arg);
  char *GetOption(string opt, int arg, int num);
  int   IsOption(string opt);

  inline const char *GetError(void) {return errorstr[error];}

  inline CLerr IsError(void) {return error;}

  inline char *GetCall(void) {return av[0];}

  inline int NumArgs(void) {return ac;}

  int  GetArgNum(string opt);
  int  GetArgNum(string opt, int num);
  bool FirstAndLast(string opt, int &f, int &l);

  inline char *GetArg(int num) {
    if (num < ac)
      return av[num]; return NULL;
  }

  void SetCLArgs(int c, char **v);

  inline bool AreCLArgsSet(void) {return av != NULL;}

 protected:
  int ac;
  char **av;

 private:
  void Init(void);
  void Delete(void);

  char *optbuf[BUF_NUM];
  bool option;
  CLerr error;

  int bcount;
};  // class kaParser

#endif  // ifndef KAPARSER_H
