/**@file LatticesMaterial.h
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

#ifndef LATTICESMATERIAL_H
#define LATTICESMATERIAL_H

#include <kaLattice.h>

//! Class for handling of orientation lattices
/*!
   lattices with information concerning distribution of 2D- and 3D-orientation are combined in a single object.
   \n\n
   \author cw,fs,idb, IBT - Universität Karlsruhe (TH)
 */

template<class T> class LatticesOrientation : virtual public nskaGlobal::DataTypeInfo {
 public:
  kaLattice<T> *phi;  //!< lattice with phi-angle of orientation
  kaLattice<T> *theta;  //!< lattice with theta-angle of orientation
#ifndef ANISOTROP2D
  kaLattice<T> *rho;  //!< lattice with rho-angle of orientation
#endif  // ifndef ANISOTROP2D

  inline LatticesOrientation() {
    phi   = NULL;
    theta = NULL;
#ifndef ANISOTROP2D
    rho = NULL;
#endif  // ifndef ANISOTROP2D
  }

  inline ~LatticesOrientation() {
    CleanUp();
  }

  void CleanUp() {
    if (phi) {delete phi; phi = NULL;}
    if (theta) {delete theta; theta = NULL;}
#ifndef ANISOTROP2D
    if (rho) {delete rho; rho = NULL;}
#endif  // ifndef ANISOTROP2D
  }

  virtual void Load(const char *name, kaLatticeHeader *temp, const char *add = "", bool noshm = false) {
#if KALATTICEDEBUG
    cerr << "LatticesOrientation<T>::Load(" << name << ", " << add << ")\n";
#endif  // if KALATTICEDEBUG
    assert(temp);
    Load(name, temp->lhWidth(), temp->lhHeight(), temp->lhDepth(), &temp->lhMatrix(), add, noshm);
  }

  virtual void Load(const char *name, IndexType lx   = -1, IndexType ly = -1, IndexType lz = -1,
                    const kaSharedMatrixN<double> *m = NULL, const char *add = "", bool noshm = false) {
#if KALATTICEDEBUG
    cerr << "LatticesOrientation<T>::Load(" << name << ", " << add << ")\n";
#endif  // if KALATTICEDEBUG
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
      if (strlen(name) > 3)
        if (!strcmp(name+strlen(name)-strlen(ext), ext))
          buf[strlen(name)-strlen(ext)] = 0;

      sprintf(cls.Name, "%s%s.phi%s", buf, add, ext);
      phi = new kaLattice<T>(cls, noshm);

      sprintf(cls.Name, "%s%s.theta%s", buf, add, ext);
      theta = new kaLattice<T>(cls, noshm);

#ifndef ANISOTROP2D
      sprintf(cls.Name, "%s%s.rho%s", buf, add, ext);
      rho = new kaLattice<T>(cls, noshm);
#endif  // ifndef ANISOTROP2D

      CheckConsistency();

      if (m) {
        phi->m   = *m;
        theta->m = *m;
#ifndef ANISOTROP2D
        rho->m = *m;
#endif  // ifndef ANISOTROP2D
      }
    } catch (...) {
      CleanUp();
      throw;
    }
  }  // Load

  inline void Remove() {
    char buf[256];

    if (phi) {
      strcpy(buf, phi->getName());
      delete phi;
      remove(buf);
      phi = NULL;
    }

    if (theta) {
      strcpy(buf, theta->getName());
      delete theta;
      remove(buf);
      theta = NULL;
    }
#ifndef ANISOTROP2D
    if (rho) {
      strcpy(buf, rho->getName());
      delete rho;
      remove(buf);
      rho = NULL;
    }
#endif  // ifndef ANISOTROP2D
  }

  inline void Save(const char *nameArg = NULL) {
#if KALATTICEDEBUG
    cerr << "LatticesOrientation::Save(" << nameArg << ")\n";
#endif  // if KALATTICEDEBUG
    char buf[256], name[256];

    const char *ext = GetSuffix<T>();

    if (nameArg) {
      strcpy(name, nameArg);
      if (strlen(name) > 3)
        if (!strcmp(name+strlen(name)-strlen(ext), ext))
          name[strlen(name)-strlen(ext)] = 0;
    }
    if (phi) {
      if (nameArg)
        sprintf(buf, "%s.phi%s", name, ext);
      phi->Save(nameArg ? buf : NULL);
    }
    if (theta) {
      if (nameArg)
        sprintf(buf, "%s.theta%s", name, ext);
      theta->Save(nameArg ? buf : NULL);
    }
#ifndef ANISOTROP2D
    if (rho) {
      if (nameArg)
        sprintf(buf, "%s.rho%s", name, ext);
      rho->Save(nameArg ? buf : NULL);
    }
#endif  // ifndef ANISOTROP2D
  }  // Save

 private:
  void CheckConsistency() {
    if (!phi)
      throw kaBaseException("Lattice phi missing");
    if (!theta)
      throw kaBaseException("Lattice theta missing");
    theta->DimMatch(phi);
#ifndef ANISOTROP2D
    rho->DimMatch(phi);
#endif  // ifndef ANISOTROP2D
  }
};  // class LatticesOrientation


//! Class for handling of material lattices
/*!
   lattices with information concerning distribution of material, 2D- and 3D-orientation are combined in a single
      object.
   \n\n
   \author cw,fs,idb, IBT - Universität Karlsruhe (TH)
 */


template<class T> class LatticesMaterial : virtual public nskaGlobal::DataTypeInfo
#ifdef ANISOTROP
                                           , public LatticesOrientation<uint8_t>
#endif  // ifdef ANISOTROP
{
 public:
  typedef kaLattice<T> LTMatTyp;

  LTMatTyp *p;  //!< lattice with material distribution
  IndexType xLattice, yLattice, zLattice, xyLattice, xyzLattice;

  inline LatticesMaterial() {
    p = NULL;
  }

  inline ~LatticesMaterial() {
    CleanUp();
  }

  void CleanUp() {
    if (p) {delete p; p = NULL;}
  }

  virtual void Load(const char *name, kaLatticeHeader *temp, const char *add = "", bool noshm = false) {
#if KALATTICEDEBUG
    cerr << "LatticesMaterial<T>::Load(" << name << ", " << add << ")\n";
#endif  // if KALATTICEDEBUG
    assert(temp);
    Load(name, temp->lhWidth(), temp->lhHeight(), temp->lhDepth(), &temp->lhMatrix(), add, noshm);
  }

  virtual void Load(const char *name, IndexType lx   = -1, IndexType ly = -1, IndexType lz = -1,
                    const kaSharedMatrixN<double> *m = NULL, const char *add = "", bool noshm = false) {
#if KALATTICEDEBUG
    cerr << "LatticesMaterial<T>::Load(" << name << ", " << add << ")\n";
#endif  // if KALATTICEDEBUG
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
      if (strlen(name) > 3)
        if (!strcmp(name+strlen(name)-strlen(ext), ext))
          buf[strlen(name)-strlen(ext)] = 0;

      sprintf(cls.Name, "%s%s%s", buf, add, ext);
      p = new LTMatTyp(cls, noshm);

      if (!p->xLattice || !p->yLattice || !p->zLattice)
        throw kaBaseException("LatticesMaterial<T>::Load(%s) failed", name);

      if (m)
        p->m = *m;

      xLattice   = p->xLattice;
      yLattice   = p->yLattice;
      zLattice   = p->zLattice;
      xyLattice  = xLattice*yLattice;
      xyzLattice = xyLattice*zLattice;

#ifdef ANISOTROP
      LatticesOrientation<uint8_t>::Load(name, lx, ly, lz, m, add, noshm);
#endif  // ifdef ANISOTROP
      CheckConsistency();
    } catch (...) {
      CleanUp();
#ifdef ANISOTROP
      LatticesOrientation<uint8_t>::CleanUp();
#endif  // ifdef ANISOTROP
      throw;
    }
  }  // Load

  inline void Remove() {
    char buf[256];

    if (p) {
      strcpy(buf, p->getName());
      delete p;
      remove(buf);
      p = NULL;
    }
#ifdef ANISOTROP
    LatticesOrientation<uint8_t>::Remove();
#endif  // ifdef ANISOTROP
  }

 private:
  void CheckConsistency() {
    if (!p)
      throw kaBaseException("Lattice Material missing");

#ifdef ANISOTROP
    phi->DimMatch(p);
#endif  // ifdef ANISOTROP
  }
};  // class LatticesMaterial

#endif  // ifndef LATTICESMATERIAL_H
