/*
 * File: kaParser.cpp
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


#include <kaParser.h>


kaParser::kaParser(void) {
  Init();
}

kaParser::kaParser(int c, char **v) {
  Init();
  SetCLArgs(c, v);
}

void kaParser::Init(void) {
  error  = ERROR_NONE;
  option = false;
  bcount = 0;

  ac = 0;
  av = NULL;

  int i;
  for (i = 0; i < BUF_NUM; i++)
    optbuf[i] = NULL;
}

void kaParser::SetCLArgs(int c, char **v) {
  Delete();
  ac = c;
  av = new char *[ac];
  int i;
  for (i = 0; i < ac; i++) {
    av[i] = new char[ARG_SIZE];
    strncpy(av[i], v[i], ARG_SIZE);
  }
}

kaParser::~kaParser(void) {
  Delete();
}

void kaParser::Delete(void) {
  if (av) {
    int i;
    for (i = 0; i < ac; i++) {
      delete[] av[i];       // aufpassen: delete[] df 22.06.2004
      av[i] = NULL;
    }
    delete[] av; av = NULL;  // same
  }
  int i;
  for (i = 0; i < BUF_NUM; i++) {
    if (optbuf[i]) delete[] optbuf[i]; optbuf[i] = NULL;
  }
  Init();
}

char *kaParser::GetOption(string opt, int arg) {
  option = false;
  error  = ERROR_NONE;

  if (arg == 0) {
    if (IsOption(opt)) {
      return (char *)"option selected";
    }

    error = NO_SUCH_OPTION;
    return NULL;
  }

  int i;
  for (i = 0; i < ac; i++) {
    if (!strcmp(opt.c_str(), av[i])) {
      option = true;
      if (i+arg < ac) {
        optbuf[bcount] = new char[ARG_SIZE];
        strcpy(optbuf[bcount], av[i+1]);
        int j;
        for (j = 2; j <= arg; j++) {
          strcat(optbuf[bcount], " "); strcat(optbuf[bcount], av[i+j]);
        }
        error =  ERROR_NONE;
        return optbuf[bcount++];
      } else {
        error = INCOMPLETE_COMMAND_LINE;
        return NULL;
      }
    }
  }

  error = NO_SUCH_OPTION;
  return NULL;
}  // kaParser::GetOption

int kaParser::IsOption(string opt) {
  option = false;
  error  = ERROR_NONE;

  int i;
  for (i = 0; i < ac; i++)
    if (!strcmp(opt.c_str(), av[i])) {
      option = true;
      return option;
    }

  error = NO_SUCH_OPTION;
  return option;
}

int kaParser::GetArgNum(string opt) {
  if (!IsOption(opt))
    return -1;

  option = false;
  error  = ERROR_NONE;

  int i;
  for (i = 0; i < ac; i++)
    if (!strcmp(opt.c_str(), av[i])) {
      option = true;
      return i;
    }

  error = NO_SUCH_OPTION;
  return option;
}

int kaParser::GetArgNum(string opt, int num) {
  if (!IsOption(opt))
    return -1;

  option = false;
  error  = ERROR_NONE;

  int i, oc = 0;
  for (i = 0; i < ac; i++)
    if (!strcmp(opt.c_str(), av[i])) {
      option = true;
      if (oc == num)
        return i;

      oc++;
    }

  error = NO_SUCH_OPTION;
  return -1;
}

int kaParser::HowMany(string opt) {
  int cnt = 0;

  option = false;
  error  = ERROR_NONE;

  int i;
  for (i = 0; i < ac; i++)
    if (!strcmp(opt.c_str(), av[i])) {
      option = true;
      cnt++;
    }

  if (cnt == 0) {
    error = NO_SUCH_OPTION;
  }
  return cnt;
}

bool kaParser::FirstAndLast(string opt, int &f, int &l) {
  // f will be the number of the first argument past opt
  // l will be the number of the first argument past opt starting with '-'

  if (!IsOption(opt))
    return false;

  int i;
  for (i = 0; i < ac; i++)
    if (!strcmp(opt.c_str(), av[i])) {
      if (i < ac-1) {f = i+1; i = ac;} else {return false;}
    }

  l = ac;
  for (i = f; i < ac; i++) {
    if (!strncmp(av[i], "-", 1)) {l = i; i = ac;}
  }
  return true;
}

char *kaParser::GetOption(string opt, int arg, int num) {
  option = false;
  error  = ERROR_NONE;

  if (num >= HowMany(opt)) {
    error = NO_FURTHER_OPTION;
    return NULL;
  }

  if (arg == 0) {
    if (IsOption(opt)) {
      return (char *)"option selected";
    }

    error = NO_SUCH_OPTION;
    return NULL;
  }

  int samecounter = 0, i;
  for (i = 0; i < ac; i++) {
    if (!strcmp(opt.c_str(), av[i])) {
      if (samecounter == num) {
        option = true;
        if (i+arg < ac) {
          optbuf[bcount] = new char[ARG_SIZE];
          strcpy(optbuf[bcount], av[i+1]);
          int j;
          for (j = 2; j <= arg; j++) {
            strcat(optbuf[bcount], " "); strcat(optbuf[bcount], av[i+j]);
          }
          error =  ERROR_NONE;
          return optbuf[bcount++];
        } else {
          error = INCOMPLETE_COMMAND_LINE;
          return NULL;
        }
      }
      samecounter++;
    }
  }

  error = NO_SUCH_OPTION;
  return NULL;
}  // kaParser::GetOption
