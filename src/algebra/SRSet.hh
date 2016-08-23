#ifndef SRSET_HH
# define SRSET_HH

#include <set>
#include <algorithm>
#include <iostream>

using namespace std;

/**********************************************************************
 Given a semigroup G, SRSet<G> is the corresponding powerset semiring,
 defined by the relations  

    A+B = A\cup B
    AB  = {ab | a\in A, b\in B}
     0  = \emptyset
     1  = { 1_G }

***********************************************************************/
template <class G> class SRSet : public set<G> {
public:
  typedef G binop_t(const G& x, const G& y);
  typedef bool predicate_t(const G& x);

  SRSet(const G& x) { this->insert(x); }

  SRSet(const set<G>& S = set<G>()) : set<G>(S.begin(),S.end()) { }

  SRSet(const SRSet &X, const SRSet &Y, binop_t compose, predicate_t include)
  {
    for(typename set<G>::const_iterator  x(X.begin()); x!= X.end(); x++)
      for(typename set<G>::const_iterator y(Y.begin()); y != Y.end(); y++){
	const G z(compose(*x,*y));
	//	cout << "product " << *x << " * " << *y << " = " << z << endl;
	if(include(z)){
	  //	  cout << z << " passes include test.\n";
	  insert(z);
	} else {
	  //	  cout << z << " doesn't pass include test.\n";
	}
      }
  }

  SRSet(int n){
    assert(n==0 || n==1);
    *this = (n==0? zero() : one());
  }

  static SRSet zero(){ return SRSet(); }
  static SRSet one() { return SRSet(G::one()); }

  friend SRSet operator+(const set<G> &X, const set<G> &Y)  {
    SRSet sum;
    set_union(X.begin(),X.end(),Y.begin(),Y.end(),inserter(sum,sum.begin()));
    return sum;
  }

  friend SRSet operator-(const set<G> &X, const set<G> &Y)  {
    SRSet sum;
    set_difference(X.begin(),X.end(),Y.begin(),Y.end(),inserter(sum,sum.begin()));
    return sum;
  }

  friend SRSet operator*(const set<G> &X, const set<G> &Y) {
    SRSet product;
    int i=1;
    for(typename set<G>::const_iterator x(X.begin()); x!=X.end(); x++)
      for(typename set<G>::const_iterator y(Y.begin()); y!=Y.end(); y++)
	product.insert((*x)*(*y));

    return product;
  }

  SRSet& operator+=(const set<G> &Y) {
    set<G>& X(*this);
    set_union(X.begin(),X.end(),Y.begin(),Y.end(),inserter(X,X.begin()));
    return *this;
  }

  SRSet& operator+=(const G& Y) { set<G>::insert(Y); return *this; }
  SRSet& operator-=(const G& Y) { set<G>::remove(Y); return *this; }

  SRSet& operator-=(const set<G> &Y) {
    set<G>& X(*this);
    set_difference(X.begin(),X.end(),Y.begin(),Y.end(),inserter(X,X.begin()));
    return *this;
  }

  friend ostream& operator<<(ostream& s, const SRSet& S)
  {
    s << "{"; 
    typename set<G>::const_iterator x(S.begin());
    if(x != S.end()){
      s << *x;
      for(++x; x != S.end(); x++) s << "," << *x;
    }
    s << "}";
    return s;
  }
};

#endif
