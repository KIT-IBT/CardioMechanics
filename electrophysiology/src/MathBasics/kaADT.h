/*
 * File: kaADT.h
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

#ifndef _kaADT_h
#define _kaADT_h

#include <kaExceptions.h>
#include <limits>

//! binary tree node
template<class T> class kaBTNode {
 public:
  T val;
  kaBTNode<T> *left;
  kaBTNode<T> *right;

  kaBTNode() {left = right = NULL;}

  kaBTNode(const T v) : val(v) {left = right = NULL;}

  // is a terminal node?
  // inline bool terminal() { return (!left)&&(!right); }

  //! print the tree starting from this node
  void print(const int level = 0) const {
    const kaBTNode<T> *node = this;

    if (node)
      node->right->print(level+1);

    for (int spc = 0; spc < level; spc++)
      cout << "   ";

    if (node)
      cout<<node->val<<'\n';
    else
      cout<<"@\n";

    if (node)
      node->left->print(level+1);
  }
};  // class kaBTNode

//! hypercube subregion for the ADT
template<unsigned int Dim = 3>
class kaADTInterval {
  double a[2*Dim];

 public:
  kaADTInterval(double *bb) {
    for (int i = 0; i < 2*Dim; i++) a[i] = bb[i];
  }

  //! begin for dimension j
  inline double get_begin(const unsigned int j) const {return a[2*j];}

  //! end for dimension j
  inline double get_end(const unsigned int j) const {return a[2*j+1];}

  //! set begin for dimension j
  inline void set_begin(const unsigned int j, const double c) {a[2*j] = c;}

  //! set end for dimension j
  inline void set_end(const unsigned int j, const double c) {a[2*j+1] = c;}

  inline double center(const unsigned int j) const {return 0.5*(a[2*j]+a[2*j+1]);}

  inline void debug_print() {
    for (int i = 0; i < Dim; i++)
      cerr<<"dim="<<i<<" begin:"<<get_begin(i)
          <<" end:"<<get_end(i)<<'\n';
  }

  /*   //!check if the coordinate is inside for dimension j */
  /*   inline bool is_inside(const double coord, const unsigned int j) */
  /*   { */
  /*     return ( coord >= get_begin(j) ) && ( coord < get_end(j) ); */
  /*   } */
};  // class kaADTInterval


//! alternating digital tree for the multidimensional search
template<class T, unsigned int Dim = 3> class kaADT {
 public:
  kaADT();
  ~kaADT();

  inline bool empty() {return !root;}

  inline unsigned int get_dimension() const {return Dim;}

  //! insert with a given tolerance
  inline bool insert(const T &v, double tol = 0);

  //! find closest object
  // with distance less than a given tolerance
  // if the tolerance is negative, perform a general search for the closest object
  inline const T find_closest(const T &v, double tol = -1.0);

  //! the last find_closest operation succeeded
  inline bool find_closest_ok() const {return close_found_flag;}

  //! set bounding box from array
  inline void set_bbox(double *bb);

  //! print tree contents (for debug purposes) \\
  // object must provide operator <<
  void print() const {
    cout<<"======== kaADT ===========\n";
    cout<<"Bounding box:\n";
    for (int i = 0; i < 2*Dim; i++) cout<<bbox[i]<<" ";
    cout<<"\n =====================================\n\n";
    if (root)
      root->print();
  }

  int size() const {return count;}

 private:
  kaBTNode<T> *root;
  double bbox[2*Dim];

  int count;

  kaBTNode<T> *terminal_found;  // free node found by find_closest
  bool place_to_left;  // place in the left node?(y/n)
  double dist2_min;  // minimal squared distance by find_closest
  T val_closest;  // closest object found

  bool close_found_flag;  // shows if the last find_closest_internal succedded

  // wrapper function to get the object's coordinates
  // defined for kaPoint and for Punkt

  template<class U>

  inline double get_coord(const U &p, const unsigned int j) {
    if (j == 0)
      return (double)p.x;
    else if (j == 1)
      return (double)p.y;
    else
      return (double)p.z;
  }

  //! wrapper function to get the distance betw the objects
  // defined for kaPoint and for  Punkt

  template<class U>

  inline double distance2(const U &p1, const U &p2) {
    return (p1.x-p2.x)*(p1.x-p2.x)+(p1.y-p2.y)*(p1.y-p2.y)+(p1.z-p2.z)*(p1.z-p2.z);
  }

  // internal search function
  inline void find_closest_internal(const T &v, kaBTNode<T> *start_node, kaADTInterval<Dim> &subregion, int level,
                                    bool good_path);

  // destructor implementation
  void cleanup(kaBTNode<T> *);
};  // class kaADT

template<class T, unsigned int Dim>
kaADT<T, Dim>::kaADT() {
  root = NULL;
  for (int i = 0; i < Dim; i++) {
    bbox[2*i]   = 0.0;
    bbox[2*i+1] = 1.0;
  }
  close_found_flag = false;
}

template<class T, unsigned int Dim>
kaADT<T, Dim>::~kaADT() {
  cleanup(root);
}

template<class T, unsigned int Dim>
void kaADT<T, Dim>::cleanup(kaBTNode<T> *start_node) {
  if (!start_node)
    return;

  if (start_node->left) {
    cleanup(start_node->left);
    start_node->left = NULL;
  }
  if (start_node->right) {
    cleanup(start_node->right);
    start_node->right = NULL;
  }

  delete start_node;
}

template<class T, unsigned int Dim>
void kaADT<T, Dim>::set_bbox(double *bb) {
  for (int i = 0; i < 2*Dim; i++)
    bbox[i] = bb[i];
}

template<class T, unsigned int Dim>
bool  kaADT<T, Dim>::insert(const T &v, double tol) {
  if (!root) {
    root = new kaBTNode<T>(v);
    count++;
    return true;
  } else {
    find_closest(v, tol);
    if (!find_closest_ok() ) {
      if (place_to_left)
        terminal_found->left = new kaBTNode<T>(v);
      else
        terminal_found->right = new kaBTNode<T>(v);
      count++;
      return true;
    }
  }
  return false;
}

template<class T, unsigned int Dim>
const T  kaADT<T, Dim>::find_closest(const T &v, const double tol) {
  if (tol < 0)
    dist2_min = numeric_limits<double>::max();
  else
    dist2_min = tol*tol;

  if (!root) {throw kaBaseException("kaADT::find_closest: the tree is empty");} else {
    kaADTInterval<Dim> subregion(bbox);
    close_found_flag = false;

    find_closest_internal(v, root, subregion, 0, true);
  }

  return val_closest;
}

template<class T, unsigned int Dim>
void kaADT<T, Dim>::find_closest_internal(const T &v, kaBTNode<T> *start_node, kaADTInterval<Dim> &subregion,
                                          int level, bool good_path) {
  //  subregion.debug_print();


  int j = level % Dim;

  // cerr<<"level: "<<level<<" DIM: "<<j<<'\n';
  // cerr<<"level: "<<level<<" good_path:"<<(int)good_path<<'\n';

  double d2 = distance2(start_node->val, v);

  if (d2 <= dist2_min) {
    dist2_min        = d2;
    val_closest      = start_node->val;
    close_found_flag = true;
  }

  // XXX I'm 90% sure the followin if blocks are bugs. The else if does not
  // what the indentation implies.

  double center  = subregion.center(j);
  double coord_j = get_coord(v, j);
  if (coord_j - dist2_min < center) // check left son
    if (start_node->left) {
      double sub_end = subregion.get_end(j);
      subregion.set_end(j, center);
      find_closest_internal(v, start_node->left, subregion, level+1, good_path && (coord_j < center) );
      subregion.set_end(j, sub_end);
    } else if (good_path && (coord_j < center) ) {
      place_to_left  = true;
      terminal_found = start_node;
      return;
    }


  if (coord_j + dist2_min >= center) // check right son
    if (start_node->right) {
      double sub_begin = subregion.get_begin(j);
      subregion.set_begin(j, center);
      find_closest_internal(v, start_node->right, subregion, level+1, good_path && (coord_j >= center) );
      subregion.set_begin(j, sub_begin);
    } else if (good_path && (coord_j >= center) ) {
      place_to_left  = false;
      terminal_found = start_node;
      return;
    }
}  // >::find_closest_internal

#endif  // ifndef _kaADT_h
