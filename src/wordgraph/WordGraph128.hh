#ifndef WORDGRAPH128_HH
#define WORDGRAPH128_HH

#include <algorithm>
#include <limits.h>
#include <map>
#include <set>
#include <stdint.h>
#include <vector>
//#include <unordered_set>

#include <algebra/SGPair.hh>
#include <algebra/SRSet.hh>
#include <algebra/ringmatrix.hh>

typedef unsigned __int128 uint128_t;

namespace WordGraph128 {
typedef uint8_t node_t;
using namespace std;

class WordGraph {
public:
  enum { zero_value = 0L, one_value = -1 };
  static const int MaxVertices = 25;
  uint128_t data;

  typedef enum { ZERO, ONE, NORMAL } element_type_t;

  element_type_t element_type() const {
    return data == zero_value ? ZERO : (data == one_value ? ONE : NORMAL);
  }

  WordGraph(const uint128_t &data = MaxVertices) : data(data) {}
  WordGraph(const uint128_t &data, const uint8_t &length)
      : data((data & ~0x1fL) | length) {}

  WordGraph(const vector<node_t> &indices) : data(indices.size()) {
    for (unsigned i = 0; i < indices.size(); i++)
      set(i + 1, indices[i]);
  }
  WordGraph(const vector<int> &indices) : data(indices.size()) {
    for (unsigned i = 0; i < indices.size(); i++)
      set(i + 1, indices[i]);
  }
  WordGraph(int n){		// Ring elements
    assert(n==0||n==1);
    *this = (n==0? zero() : one());
  }
  
  static WordGraph one() { return WordGraph(one_value); }
  static WordGraph zero() { return WordGraph(zero_value); }

  // i runs from 1,...,25, value from 1: 0: no arc, 1,..,25: actual indices
  inline void set(uint8_t i, uint8_t value) {
    data &= ~(0x1fL << (i * 5));
    data |= (value & 0x1fL) << (i * 5);
  }

  // i runs from 0,...,24, value from 0,...,24 -- can't set "no arc"
  inline void setb(uint8_t i, uint8_t value) {
    data &= ~(0x1fL << ((i + 1) * 5));
    data |= ((value + 1) & 0x1fL) << ((i + 1) * 5);
  }

  inline uint8_t get(uint8_t i) const {
    return ((data & ~0x1fL) >> i * 5) & 0x1fL;
  }
  inline uint8_t length() const { return data & 0x1fL; }
  inline void setlength(int length) {
    data &= ~0x1fL;
    data |= (length & 0x1fL);
  }
  inline void setdimension(int rows, int cols) { setlength(cols); }

  inline bool operator==(const WordGraph &y) const { return data == y.data; }
  inline bool operator>(const WordGraph &y) const { return data > y.data; }
  inline bool operator<=(const WordGraph &y) const { return data <= y.data; }
  inline bool operator<(const WordGraph &y) const { return data < y.data; }

  WordGraph operator*(const WordGraph &g) const;

  int countarcs(int max = INT_MAX) const;
  int countcycles(int max = INT_MAX) const;
  bool synchronizes() const;

  uint128_t number() const;
  static WordGraph fromnumber();

  static WordGraph plain_compose(const WordGraph &g, const WordGraph &h) {
    return g * h;
  }

  static bool synchronize_test(const WordGraph &g) { return !g.synchronizes(); }
  static bool true_test(const WordGraph &g) { return true; }

  friend std::ostream &operator<<(std::ostream &s, const WordGraph &g);
};

class WGSet : public set<uint128_t> {
public:
  typedef WordGraph wordgraph_t;
  typedef uint128_t     value_t;
  typedef set<uint128_t> intset;
  typedef WordGraph binop_t(const WordGraph &x, const WordGraph &y);
  typedef bool predicate_t(const WordGraph &x);

  WGSet(const intset &S = intset()) : intset(S.begin(), S.end()) {}
  WGSet(const WordGraph &x) { insert(x.number()); }
  WGSet(const WGSet &X, const WGSet &Y, binop_t compose, predicate_t include) {
    for (intset::const_iterator x(X.begin()); x != X.end(); x++)
      for (intset::const_iterator y(Y.begin()); y != Y.end(); y++) {
        const WordGraph z(compose(*x, *y));
        //	cout << "product " << WordGraph(*x) << " * " << WordGraph(*y) <<
        //" = " << z << endl;
        if (include(z)) {
          //	  cout << z << " passes include test.\n";
          insert(z.number());
        } else {
          //	  cout << z << " doesn't pass include test.\n";
        }
      }
  }

  WordGraph doublecycle() const {
    for (const_iterator w(begin()); w != end(); w++)
      if (WordGraph(*w).countcycles() >= 2)
        return WordGraph(*w);
    return WordGraph::zero();
  }

  WGSet star() const {
    WGSet xstar(one()), xnth(one());

    if (!empty()) {
      do {
        xstar += xnth;
        xnth *= *this;
        xnth -= xstar;
      } while (!xnth.empty());
    }
    return xstar;
  }

  static WGSet zero() { return WGSet(); }
  static WGSet one() { return WGSet(WordGraph::one()); }

  friend WGSet operator*(const WGSet &S, const WGSet &T) {
    return WGSet(S, T, WordGraph::plain_compose, WordGraph::synchronize_test);
  }

  WGSet &operator*=(const WGSet &T) {
    *this = (*this) * T;
    return *this;
  }

  friend WGSet operator+(const set<uint128_t> &X, const set<uint128_t> &Y) {
    WGSet sum;
    set_union(X.begin(), X.end(), Y.begin(), Y.end(),
              inserter(sum, sum.begin()));
    return sum;
  }

  WGSet &operator+=(const set<uint128_t> &Y) {
    set<uint128_t> &X(*this);
    set_union(X.begin(), X.end(), Y.begin(), Y.end(), inserter(X, X.begin()));
    return *this;
  }

  WGSet &operator+=(const WordGraph &x) {
    insert(x.number());
    return *this;
  }

  friend WGSet operator-(const set<uint128_t> &X, const set<uint128_t> &Y) {
    WGSet sum;
    set_difference(X.begin(), X.end(), Y.begin(), Y.end(),
                   inserter(sum, sum.begin()));
    return sum;
  }

  WGSet &operator-=(const set<uint128_t> &Y) {
    set<uint128_t> &X(*this);
    set_difference(X.begin(), X.end(), Y.begin(), Y.end(),
                   inserter(X, X.begin()));
    return *this;
  }

  friend ostream &operator<<(ostream &s, const WGSet &S) {
    s << "{";
    set<uint128_t>::const_iterator x(S.begin());
    if (x != S.end()) {
      s << WordGraph(*x);
      for (++x; x != S.end(); x++)
        s << "," << WordGraph(*x);
    }
    s << "}";
    return s;
  }
};

typedef ringmatrix<WGSet> WGMatrix;
typedef SRSet<word> wordset;
}




#endif
