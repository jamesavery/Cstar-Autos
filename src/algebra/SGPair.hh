#ifndef SGPAIR_HH
# define SGPAIR_HH 
#include <list>
#include <string>
#include <auxiliary.hh>

using namespace std;

/**********************************************************************

 Given semigroups G and H, SGPair<G,H> is the direct sum semigroup 
 defined by

  (a,b)*(c,d) = (ac,bd)
  1_{GH} = (1_G,1_H)

 **********************************************************************/
template <typename G, typename H> class SGPair : public pair<G,H> {
public:

  static SGPair zero(){ return SGPair(G::zero(), H::zero()); }
  static SGPair one() { return SGPair(G::one(), H::one()); }

  SGPair(const G& a, const H& b) : pair<G,H>(a,b) {}
  SGPair(int n){
    assert(n==0 || n==1);
    *this = (n==0? zero() : one());
  }

  SGPair operator*(const SGPair& b) const { return SGPair(this->first*b.first, this->second*b.second); }
  SGPair& operator*=(const SGPair& b) { this->first *= b.first; this->second *= b.second; return *this; }
};


/**********************************************************************

   Given an alphabet A, SGString<A> is the free monoid over A.

 **********************************************************************/
template <typename A> class SGString: public basic_string<A> {
 public:
  static SGString one(){ return basic_string<A>(); }
  SGString(const A &x) : basic_string<A>(1,x) {}
  SGString(const basic_string<A>& s) : basic_string<A>(s) {}
  SGString  operator*(const SGString& b) const { return (*this + b);  }
  SGString& operator*=(const SGString& b)      { return basic_string<A>::operator+=(b); }

  friend ostream& operator<<(ostream& S, const SGString& G){S << vector<A>(G.begin(),G.end()); return S;}
};



// TODO: Requires compiling -- should this be put in a .cc-file?
// #include "../auxiliary.hh"
// container_output(SGString);

typedef SGString<char> word;

#endif
