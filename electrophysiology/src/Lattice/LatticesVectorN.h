/*
 * File: LatticesVectorN.h
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


#ifndef LATTICESVECTORN_H
#define LATTICESVECTORN_H

#include <kaLattice.h>

//! Class for handling of Nd vector lattices.
/*!
   lattices with N vector components are combined in a single object.

 */

template<class T, int N> class LatticesVectorN : virtual public nskaGlobal::DataTypeInfo {
 public:
  kaLattice<T> *l[N];
  IndexType xLattice, yLattice, zLattice, xyLattice, xyzLattice;

  inline LatticesVectorN() {for (int i = 0; i < N; i++) l[i] = NULL; }

  inline ~LatticesVectorN() {CleanUp();}

  inline void CleanUp() {
    for (int i = 0; i < N; i++) {
      delete l[i];
      l[i] = NULL;
    }
  }

  bool isnan() {
    for (int i = 0; i < N; i++) {
      T *pu = &l[i]->lat[0], *pue = pu+xyzLattice;

      while (pu < pue) {  // df 15.03.2006: war (pu<pe)!!!! Verdammt noch mal...
        if (::isnan(*pu))
          return true;

        pu++;
      }
    }
    return false;
  }

  virtual void Load(const char *name, kaLatticeHeader *temp, const char *add = "") {
#if KALATTICEDEBUG
    cerr << "LatticesVectorN<T>::Load(" << name << ", " << add << ")\n";
#endif  // if KALATTICEDEBUG
    assert(temp);
    Load(name, temp->lhWidth(), temp->lhHeight(), temp->lhDepth(), &temp->lhMatrix(), add);
  }

  virtual void Load(const char *name, int lx = -1, int ly = -1, int lz = -1, const kaSharedMatrixN<double> *m = NULL,
                    const char *add          = "") {
    CleanUp();

    try {
      kaLatticeCreate cls;
      if (lx >= 0)
        cls.xLattice = lx;
      if (ly >= 0)
        cls.yLattice = ly;
      if (lz >= 0)
        cls.zLattice = lz;

      const char *ext = GetSuffix<T>();
      char buf[256];
      strcpy(buf, name);
      if (strlen(buf) > 3)
        if (!strcmp(name+strlen(buf)-strlen(ext), ext))
          buf[strlen(buf)-strlen(ext)] = 0;

      for (int i = 0; i < N; i++) {
        if (N == 1)
          sprintf(cls.Name, "%s%s%s", buf, add, ext);
        else
          sprintf(cls.Name, "%s%s.%d%s", buf, add, i, ext);
        l[i] = new kaLattice<T>(cls);
        if (lx > 0)
          if (l[i]->xLattice != lx) throw DimensionMismatch(buf, lx, ly, lz);
        if (ly > 0)
          if (l[i]->yLattice != ly) throw DimensionMismatch(buf, lx, ly, lz);
        if (lz > 0)
          if (l[i]->zLattice != lz) throw DimensionMismatch(buf, lx, ly, lz);

        if (m)
          l[i]->m = *m;
      }
      CheckConsistency();

      xLattice   = l[0]->xLattice;
      yLattice   = l[0]->yLattice;
      zLattice   = l[0]->zLattice;
      xyLattice  = xLattice*yLattice;
      xyzLattice = xyLattice*zLattice;
    } catch (...) {
      CleanUp();
      throw;
    }
  }  // Load

  inline void Save(const char *nameArg = NULL) {
#if KALATTICEDEBUG
    cerr << "LatticesVectorN<T>::Save(" << (nameArg ? nameArg : "NULL") << ")\n";
#endif  // if KALATTICEDEBUG
    char buf[256], name[256];
    const char *ext = GetSuffix<T>();

    if (nameArg) {
      strcpy(name, nameArg);
      if (strlen(name) > 3)
        if (!strcmp(name+strlen(name)-strlen(ext), ext))
          name[strlen(name)-strlen(ext)] = 0;
    }
    for (int i = 0; i < N; i++)
      if (l[i]) {
        if (nameArg)
          if (N == 1)
            sprintf(buf, "%s%s", name, ext);
          else
            sprintf(buf, "%s.%d%s", name, i, ext);
        l[i]->Save(nameArg ? buf : NULL);
      }
  }

  inline void FillWith(T d) {
#if KALATTICEDEBUG
    cerr << "LatticesVectorN<T>::FillWith(" << d << ")\n";
#endif  // if KALATTICEDEBUG
    for (int i = 0; i < N; i++)
      if (l[i])
        l[i]->lat.FillWith(d);
  }

  inline void Remove() {
    for (int i = 0; i < N; i++)
      if (l[i]) {
        char buf[256];
        strcpy(buf, l[i]->getName());
        delete l[i];
        remove(buf);
        l[i] = NULL;
      }
  }

 private:
  void CheckConsistency() {
    for (int i = 0; i < N; i++) {
      if (!l[i])
        throw kaBaseException("LatticesVectorN: Lattice %d missing", i);
      l[0]->DimMatch(l[i]);
    }
  }
};  // class LatticesVectorN

#endif  // ifndef LATTICESVECTORN_H
