/*
 * File: kaSharedObjectHandling.h
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

#ifndef KASHAREDOBJECTHANDLING
#define KASHAREDOBJECTHANDLING

#include <iostream>
#include <dlfcn.h>
#include <kaExceptions.h>
#include <mach-o/dyld.h>

static void CloseLibrary(void *handle) {  // dw:oldVersion: static void CloseLibrary(void *handle){
#ifdef KADEBUG
  cerr << "closing library in " << handle << "\n";
#endif  // ifdef KADEBUG
  dlclose(handle);
}

static void *LoadLibrary(char *fileName) {
  void *handle = dlopen(fileName, RTLD_NOW);

  if (handle) {
#if KADEBUG
    cerr<<"library " << fileName << " successfully loaded to " << handle << endl;
#endif  // if KADEBUG
    return handle;
  } else {
    throw kaBaseException(dlerror());
  }
}

static void *LoadSymbol(void *handle, string symbol) {  // dw: oldVersion: static void *LoadSymbol(void *handle, char
                                                        // *symbol){
#if KADEBUG
  cerr<<"loading symbol '" << symbol << "'\n";
#endif  // if KADEBUG
  void *retValue = dlsym(handle, symbol.c_str());

  if (!retValue) {
    CloseLibrary(handle);
    throw kaBaseException("cannot load symbol '%s': %s\n", symbol.c_str(), dlerror());
  }

  return retValue;
}

#endif  // ifndef KASHAREDOBJECTHANDLING
