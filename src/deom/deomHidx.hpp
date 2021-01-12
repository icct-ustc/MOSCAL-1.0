/**
 * git clone https://github.com/hou-dao/deom.git
 * ---
 * Written by Houdao Zhang 
 * mailto: houdao@connect.ust.hk
 */
#ifndef DEOMHIDX_H_
#define DEOMHIDX_H_

#include "armadillo"
#include "json11.hpp"
#include "trie.hpp"
#include <cstdio>
#include <map>
#include <string>


using namespace std;
using namespace arma;
using namespace json11;

class hnod {

public:
  // data
  cx_double gams;
  ivec key;
  int info;

  // constructor
  hnod() : gams(0) {}

  // constructor
  hnod(const hnod &rhs) : gams(rhs.gams), key(rhs.key), info(rhs.info) {}

  // constructor
  hnod(const cx_double _gams, const ivec &_key) : gams(_gams), key(_key) { info = 0; }

  // constructor
  hnod(const cx_double _gams, const ivec &_key, const int &_info) : gams(_gams), key(_key), info(_info) {}

  // overload = operator
  hnod &operator=(const hnod &rhs) {
    if (this != &rhs) {
      gams = rhs.gams;
      key = rhs.key;
      info = rhs.info;
    }
    return *this;
  };

  ~hnod(){};
};

typedef field<hnod> hkey;

class hidx {

public:
  int nind;
  int lmax;
  int nmax;
  int lddo;
  int nddo;
  double ferr;
  hkey keys;
  Trie tree;
  cx_vec expn;

  hidx(const Json &json) {

    const string expnFile = json["expnFile"].string_value();
    if (expn.load(expnFile, arma_ascii)) {
    } else {
      printf("Error: expn is not loaded!\n");
    }

    nind = expn.n_rows;
    lmax = json["lmax"].int_value();
    nmax = json["nmax"].int_value();
    ferr = json["ferr"].number_value();

    nddo = 1;
    lddo = 0;

    int nmax_kl = static_cast<int>(get_nmax(nind, lmax));
    if (nmax_kl < nmax) {
      nmax = nmax_kl;
    }

    keys.set_size(nmax);
    keys(0) = hnod(0, zeros<ivec>(1));

    tree.init(lmax + 1);
    if (!tree.try_insert(keys(0).key, 0)) {
      printf("Error in inserting first key to tree!\n");
    }

    printf("$InitHidx\n");
    printf("nind = %d\n", nind);
    printf("lddo = %d\n", lddo);
    printf("nddo = %d\n", nddo);
    printf("lmax = %d\n", lmax);
    printf("nmax = %d\n", nmax);
    printf("ferr = %g\n", ferr);
    printf("$InitHidx\n\n");
  }

  hidx(const hidx &rhs) : nind(rhs.nind), lmax(rhs.lmax), nmax(rhs.nmax),
                          lddo(rhs.lddo), nddo(rhs.nddo), ferr(rhs.ferr), keys(rhs.keys),
                          tree(rhs.tree), expn(rhs.expn) {}

  hidx(const int _nind, const int _lmax, const int _nmax, const int _lddo, const int _nddo,
       const double _ferr, const hkey &_keys, const Trie &_tree, const cx_vec &_expn) : nind(_nind), lmax(_lmax), nmax(_nmax),
                                                                                        lddo(_lddo), nddo(_nddo), ferr(_ferr), keys(_keys),
                                                                                        tree(_tree), expn(_expn) {}

  ~hidx() {}

  hidx &operator=(const hidx &rhs) {
    if (this != &rhs) {
      nind = rhs.nind;
      lmax = rhs.lmax;
      nmax = rhs.nmax;
      lddo = rhs.lddo;
      nddo = rhs.nddo;
      ferr = rhs.ferr;
      keys = rhs.keys;
      tree = rhs.tree;
      expn = rhs.expn;
    }
    return *this;
  }

  /* The maximal number of ddos given K and L
     *
     */
  unsigned long get_nmax(const int K, const int L) const {
    unsigned long ntot = 1;
    for (int k = 1; k <= K; ++k) {
      ntot *= L + k;
      ntot /= k;
      if (ntot > 300000) {
        printf("Be careful! too many elements");
        return nmax;
      }
    }
    return ntot;
  }

  ivec gen_key(const ivec &str, const int pos, const int chg) const {
    int len0 = str.n_rows;
    ivec key;
    if (len0 >= pos + 1 && chg > 0) {
      key = str;
      key(pos) += chg;
    } else if (len0 < pos + 1 && chg > 0) {
      key = zeros<ivec>(pos + 1);
      for (int i = 0; i < len0; ++i) {
        key(i) = str(i);
      }
      key(pos) += chg;
    } else if (len0 >= pos + 1 && chg < 0 && str(pos) + chg >= 0) {
      ivec tmp(str);
      tmp(pos) += chg;
      int npos = len0;
      while (tmp(npos - 1) == 0 && npos != 1) {
        npos -= 1;
      }
      key = tmp.head(npos);
    } /* else {
            cout << "Error in generating new key" << endl;
        } */
    return key;
  }

  void init_hierarchy() {
    int ni = 0;
    int nf = nddo;
    for (int l = 0; l < lmax; ++l) {
      for (int iddo = ni; iddo < nf; ++iddo) {
        const hnod &nod = keys(iddo);
        ivec key0(nod.key);
        for (int mp = 0; mp < nind; ++mp) {
          ivec key1 = gen_key(key0, mp, 1);
          auto *ptr = tree.find(key1);
          if (ptr && ptr->rank >= 0) {
            // printf("Error! Insert node fails\n");
          } else {
            tree.try_insert(key1, nddo);
            keys(nddo) = hnod(nod.gams + expn(mp), key1);
            nddo += 1;
          }
        }
      }
      ni = nf;
      nf = nddo;
    }
    lddo = lmax;
    nmax = nddo;
  }

  ivec gen_nbar(ivec nvec) {
    int k = 0;
    int len = nvec.n_rows;
    ivec nbar = zeros<ivec>(len + 1);
    while (k < len) {
      if (abs(imag(expn(k))) < 1.0e-20) {
        nbar(k) = nvec(k);
        k += 1;
      } else {
        if (k + 1 < len)
          nbar(k) = nvec(k + 1);
        else
          nbar(k) = 0;
        nbar(k + 1) = nvec(k);
        k += 2;
      }
    }
    return (nbar(len) == 0) ? nbar.head_rows(len) : nbar;
  }

  void status() {
    ivec tcount = zeros<ivec>(lddo + 1);
    for (int i = 0; i < nddo; ++i) {
      int l = tree.find(keys(i).key)->tier;
      tcount(l) += 1;
    }
    printf("DDO index histogram\n");
    for (int l = 0; l <= lddo; ++l) {
      int pcen = (int)((double)tcount(l) / (double)nddo * 200);
      printf("Tier %3d: ", l);
      for (int i = 0; i < pcen; ++i) {
        printf("*");
      }
      printf("(%d)\n", (int)tcount(l));
    }
  }
};

#endif
