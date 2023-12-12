/*
 * File: LatticesVector3.h
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


#ifndef LATTICESVECTOR3_H
#define LATTICESVECTOR3_H

#include <kaLattice.h>

//! Class for handling of 3d vector lattices.
/*!
   lattices with 3 vector components are combined in a single object.

 */

template<class T> class LatticesVector3 : virtual public nskaGlobal::DataTypeInfo {
 public:
  kaLattice<T> *x, *y, *z;
  IndexType xLattice, yLattice, zLattice, xyLattice, xyzLattice;

  inline LatticesVector3() {x = y = z = NULL;}

  inline ~LatticesVector3() {CleanUp();}

  inline void CleanUp() {
    delete x; x = NULL;
    delete y; y = NULL;
    delete z; z = NULL;
  }

  bool lv3isnan() {
    T *pux = &x->lat[0], *puxe = pux+xyzLattice;
    T *puy = &y->lat[0];
    T *puz = &z->lat[0];

    while (pux < puxe) {
      if (::isnan(*pux) || ::isnan(*puy) || ::isnan(*puz))
        break;
      pux++;
      puy++;
      puz++;
    }

    return pux != puxe;
  }

  virtual void Load(const char *name, kaLatticeHeader *temp, const char *add = "") {
#if KALATTICEDEBUG
    cerr << "LatticesVector3<T>::Load(" << name << ", " << add << ")\n";
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
          buf[strlen(buf)-strlen(ext)] = 0; // remove the suffix the dot included (ex: .dlat in Mat.dlat).

      sprintf(cls.Name, "%s%s.x%s", buf, add, ext);
      x = new kaLattice<T>(cls);
      if (lx > 0)
        if (x->xLattice != lx) throw DimensionMismatch(buf, lx, ly, lz);
      if (ly > 0)
        if (x->yLattice != ly) throw DimensionMismatch(buf, lx, ly, lz);
      if (lz > 0)
        if (x->zLattice != lz) throw DimensionMismatch(buf, lx, ly, lz);

      sprintf(cls.Name, "%s%s.y%s", buf, add, ext);
      y = new kaLattice<T>(cls);

      sprintf(cls.Name, "%s%s.z%s", buf, add, ext);
      z = new kaLattice<T>(cls);

      CheckConsistency();

      if (m) {
        x->m = *m;
        y->m = *m;
        z->m = *m;
      }
      xLattice   = x->xLattice;
      yLattice   = x->yLattice;
      zLattice   = x->zLattice;
      xyLattice  = xLattice*yLattice;
      xyzLattice = xyLattice*zLattice;
    } catch (...) {
      CleanUp();
      throw;
    }
  }  // Load

  inline void Save(const char *nameArg = NULL) {
#if KALATTICEDEBUG
    cerr << "LatticesVector3<T>::Save(" << (nameArg ? nameArg : "NULL") << ")\n";
#endif  // if KALATTICEDEBUG
    char buf[256], name[256];
    const char *ext = GetSuffix<T>();

    if (nameArg) {
      strcpy(name, nameArg);
      if (strlen(name) > 3)
        if (!strcmp(name+strlen(name)-strlen(ext), ext))
          name[strlen(name)-strlen(ext)] = 0;
    }
    if (x) {
      if (nameArg)
        sprintf(buf, "%s.x%s", name, ext);
      x->Save(nameArg ? buf : NULL);
    }
    if (y) {
      if (nameArg)
        sprintf(buf, "%s.y%s", name, ext);
      y->Save(nameArg ? buf : NULL);
    }
    if (z) {
      if (nameArg)
        sprintf(buf, "%s.z%s", name, ext);
      z->Save(nameArg ? buf : NULL);
    }
  }  // Save

  inline void FillWith(T d) {
#if KALATTICEDEBUG
    cerr << "LatticesVector3<T>::FillWith(" << d << ")\n";
#endif  // if KALATTICEDEBUG
    if (x)
      x->lat.FillWith(d);
    if (y)
      y->lat.FillWith(d);
    if (z)
      z->lat.FillWith(d);
  }

  inline void Remove() {
    char buf[256];

    if (x) {
      strcpy(buf, x->getName());
      delete x;
      remove(buf);
      x = NULL;
    }

    if (y) {
      strcpy(buf, y->getName());
      delete y;
      remove(buf);
      y = NULL;
    }
    if (z) {
      strcpy(buf, z->getName());
      delete z;
      remove(buf);
      z = NULL;
    }
  }

 private:
  void CheckConsistency() {
    if (!x)
      throw kaBaseException("LatticesVector3: Lattice x missing");
    if (!y)
      throw kaBaseException("LatticesVector3: Lattice y missing");
    if (!z)
      throw kaBaseException("LatticesVector3: Lattice z missing");

    x->DimMatch(z);
    x->DimMatch(y);
  }
};  // class LatticesVector3

#endif  // ifndef LATTICESVECTOR3_H
