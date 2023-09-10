/**@file LatticesSymmetricTensor2.h
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

#ifndef LATTICESSYMMETRICTENSOR2_H
#define LATTICESSYMMETRICTENSOR2_H

#include <kaLattice.h>

//! Class for handling of lattices with elements of type 2. order symmetric tensor
/*!
   lattices with 6 tensor elements are combined in a single object.
   \n\n
   \author cw,fs,idb, IBT - Universität Karlsruhe (TH)
 */

template<class T> class LatticesSymmetricTensor2 : virtual public nskaGlobal::DataTypeInfo {
 public:
  kaLattice<T> *a, *b, *c, *d, *e, *f;
  IndexType xLattice, yLattice, zLattice, xyLattice, xyzLattice;

  inline LatticesSymmetricTensor2() {a = b = c = d = e = f = NULL;}

  inline ~LatticesSymmetricTensor2() {CleanUp();}

  void CleanUp() {
    delete a; a = NULL;
    delete b; b = NULL;
    delete c; c = NULL;
    delete d; d = NULL;
    delete e; e = NULL;
    delete f; f = NULL;
  }

  virtual void Load(const char *name, kaLatticeHeader *temp, const char *add = "") {
#if KALATTICEDEBUG
    cerr << "LatticesSymmetricTensor2<T>::Load(" << name << ", " << add << ")\n";
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

      sprintf(cls.Name, "%s%s.a%s", buf, add, ext);
      a = new kaLattice<T>(cls);

      if (lx >= 0) if ((a->xLattice != cls.xLattice) || !a->xLattice)
          throw kaBaseException("LatticesSymmetricTensor2<T>::Load(%s) failed", name);
      if (ly >= 0) if ((a->yLattice != cls.yLattice) || !a->yLattice)
          throw kaBaseException("LatticesSymmetricTensor2<T>::Load(%s) failed", name);
      if (lz >= 0) if ((a->zLattice != cls.zLattice) || !a->zLattice)
          throw kaBaseException("LatticesSymmetricTensor2<T>::Load(%s) failed", name);

      sprintf(cls.Name, "%s%s.b%s", buf, add, ext);
      b = new kaLattice<T>(cls);

      sprintf(cls.Name, "%s%s.c%s", buf, add, ext);
      c = new kaLattice<T>(cls);

      sprintf(cls.Name, "%s%s.d%s", buf, add, ext);
      d = new kaLattice<T>(cls);

      sprintf(cls.Name, "%s%s.e%s", buf, add, ext);
      e = new kaLattice<T>(cls);

      sprintf(cls.Name, "%s%s.f%s", buf, add, ext);
      f = new kaLattice<T>(cls);

      CheckConsistency();

      if (m) {
        a->m = *m;
        b->m = *m;
        c->m = *m;
        d->m = *m;
        e->m = *m;
        f->m = *m;
      }
      xLattice   = a->xLattice;
      yLattice   = a->yLattice;
      zLattice   = a->zLattice;
      xyLattice  = xLattice*yLattice;
      xyzLattice = xyLattice*zLattice;
    } catch (...) {
      CleanUp();
      throw;
    }
  }  // Load

  inline void Save(const char *nameArg = NULL) {
#if KALATTICEDEBUG
    cerr << "LatticesSymmetricTensor2<T>::Save(" << nameArg << ")\n";
#endif  // if KALATTICEDEBUG
    char buf[256], name[256];
    const char *ext = GetSuffix<T>();

    if (nameArg) {
      strcpy(name, nameArg);
      if (strlen(name) > 3)
        if (!strcmp(name+strlen(name)-strlen(ext), ext))
          name[strlen(name)-strlen(ext)] = 0;
    }
    if (a) {
      if (nameArg)
        sprintf(buf, "%s.a%s", name, ext);
      a->Save(nameArg ? buf : NULL);
    }
    if (b) {
      if (nameArg)
        sprintf(buf, "%s.b%s", name, ext);
      b->Save(nameArg ? buf : NULL);
    }
    if (c) {
      if (nameArg)
        sprintf(buf, "%s.c%s", name, ext);
      c->Save(nameArg ? buf : NULL);
    }
    if (d) {
      if (nameArg)
        sprintf(buf, "%s.d%s", name, ext);
      d->Save(nameArg ? buf : NULL);
    }
    if (e) {
      if (nameArg)
        sprintf(buf, "%s.e%s", name, ext);
      e->Save(nameArg ? buf : NULL);
    }
    if (f) {
      if (nameArg)
        sprintf(buf, "%s.f%s", name, ext);
      f->Save(nameArg ? buf : NULL);
    }
  }  // Save

  inline void FillWith(T d) {
#if KALATTICEDEBUG
    cerr << "LatticesSymmetricTensor2<T>::FillWith(" << d << ")\n";
#endif  // if KALATTICEDEBUG
    if (a)
      a->lat.FillWith(d);
    if (b)
      b->lat.FillWith(d);
    if (c)
      c->lat.FillWith(d);
    if (d)
      d->lat.FillWith(d);
    if (e)
      e->lat.FillWith(d);
    if (f)
      f->lat.FillWith(d);
  }

  inline void Remove() {
    char buf[256];

    if (a) {
      strcpy(buf, a->Name);
      delete a;
      remove(buf);
      a = NULL;
    }
    if (b) {
      strcpy(buf, b->Name);
      delete b;
      remove(buf);
      b = NULL;
    }
    if (c) {
      strcpy(buf, c->Name);
      delete c;
      remove(buf);
      c = NULL;
    }
    if (d) {
      strcpy(buf, d->Name);
      delete d;
      remove(buf);
      d = NULL;
    }
    if (e) {
      strcpy(buf, e->Name);
      delete e;
      remove(buf);
      e = NULL;
    }
    if (f) {
      strcpy(buf, f->Name);
      delete f;
      remove(buf);
      f = NULL;
    }
  }  // Remove

 private:
  void CheckConsistency() {
    if (!a)
      throw kaBaseException("LatticesSymmetricTensor2: Lattice a missing");
    if (!b)
      throw kaBaseException("LatticesSymmetricTensor2: Lattice b missing");
    if (!c)
      throw kaBaseException("LatticesSymmetricTensor2: Lattice c missing");
    if (!d)
      throw kaBaseException("LatticesSymmetricTensor2: Lattice d missing");
    if (!e)
      throw kaBaseException("LatticesSymmetricTensor2: Lattice e missing");
    if (!f)
      throw kaBaseException("LatticesSymmetricTensor2: Lattice f missing");

    a->DimMatch(b);
    a->DimMatch(c);
    a->DimMatch(d);
    a->DimMatch(e);
    a->DimMatch(f);
  }
};  // class LatticesSymmetricTensor2


#endif  // ifndef LATTICESSYMMETRICTENSOR2_H
