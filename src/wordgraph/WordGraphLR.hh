#ifndef WORDGRAPHLR_HH
# define WORDGRAPHLR_HH


#include <map>
#include <vector>
#include <set>
#include <algorithm>
#include <stdint.h>
#include <limits.h>

#include <algebra/ringmatrix.hh>
#include <algebra/digraph.hh>
#include <algebra/SGPair.hh>
#include <algebra/SRSet.hh>

// This wordgraph implementation allows any number of nodes <= 2^32,
// but requires that the bipartite graph be fan-in free.  It is much,
// much slower than WordGraph64, which should be used when the source-
// and range-set of the bipartite graph has <= 15 nodes.
namespace WordGraphLR {

  typedef int node_t;
  using namespace std;
  
  class WordGraph {
  public:
    // TODO: Figure out what is a better one-value, ~0xff (in band) or -1 (out of band).
    vector<node_t> data;
    template<class Archive> void serialize(Archive & archive){ archive(data); }

    typedef enum {ZERO,ONE,NORMAL} element_type_t;
    element_type_t e_type; 

    element_type_t element_type() const { return e_type; }

    WordGraph(const WordGraph& Y) : data(Y.data.begin(),Y.data.end()), e_type(Y.e_type) {    }
    WordGraph(const int length, const element_type_t e_type) : data(length+1), e_type(e_type) {}
    WordGraph(const vector<node_t>& data_, const element_type_t e_type = NORMAL) : data(data_.size()+1), e_type(e_type) {
      for(int i=0;i<data_.size();i++) data[i+1] = data_[i];
    }
    WordGraph(const node_t value, const int length) : data(length+1,value), e_type(NORMAL) {}
    WordGraph(int n=0){		// Ring elements
      assert(n==0||n==1);
      *this = (n==0? zero() : one());
    }

    static WordGraph one() { return WordGraph(0,ONE);  }
    static WordGraph zero(){ return WordGraph(0,ZERO); }

    // i runs from 1,...,n, value from 1:
    // 0: no arc, 1,..,15: actual indices
    void set(uint8_t i, uint8_t value){  data[i] = value; }
    // i runs from 0,...,n-1, value from 0,...,n-1 -- can't set "no arc"
    void setb(uint8_t i, uint8_t value){ data[i+1] = value+1;  }

  
    uint8_t get(uint8_t i) const   { return data[i]; }
    uint8_t length() const         { return data.size()-1; }
    void setlength(int length)     { data.resize(length+1); }
    void setdimension(int rows, int cols){ setlength(cols); }

    bool operator==(const WordGraph& y) const { return e_type == y.e_type && data == y.data; }
    bool operator>(const WordGraph& y)  const { return e_type > y.e_type  || (e_type == y.e_type && data  > y.data); }
    bool operator<=(const WordGraph& y) const { return e_type <= y.e_type || (e_type == y.e_type && data <= y.data); }
    bool operator<(const WordGraph& y)  const { return e_type < y.e_type  || (e_type == y.e_type && data  < y.data); }

    WordGraph operator*(const WordGraph& g) const;

    int countarcs(int max=INT_MAX) const;
    int countcycles(int max=INT_MAX) const;
    bool synchronizes() const;

    static WordGraph plain_compose(const WordGraph& g, const WordGraph& h) 
    {
      return g*h;
    }

    static bool synchronize_test(const WordGraph& g) { 
      return !g.synchronizes(); 
    }
    static bool true_test(const WordGraph& g){ return true; }

    friend std::ostream& operator<<(std::ostream& s, const WordGraph& g);
  };


class WGSet : public set<WordGraph> {
public:
  typedef WordGraph wordgraph_t;
  typedef WordGraph value_t;
  
  typedef WordGraph binop_t(const WordGraph& x, const WordGraph& y);
  typedef bool predicate_t(const WordGraph& x);

  WGSet(const set<WordGraph>& S = set<WordGraph>()) : set<WordGraph>(S) { }
  WGSet(const WordGraph& x) { insert(x); } 
  WGSet(const WGSet &X, const WGSet &Y, binop_t compose, predicate_t include) 
  {
    for(set<WordGraph>::const_iterator  x(X.begin()); x!= X.end(); x++)
      for(set<WordGraph>::const_iterator y(Y.begin()); y != Y.end(); y++){
	const WordGraph z(compose(*x,*y));
	//	cout << "product " << WordGraph(*x) << " * " << WordGraph(*y) << " = " << z << endl;
	if(include(z)){
	  //	  cout << z << " passes include test.\n";
	  insert(z);
	} else {
	  //	  cout << z << " doesn't pass include test.\n";
	}
      }
  }

  WordGraph doublecycle() const 
  {
    for(const_iterator w(begin()); w!=end();w++)
      if(WordGraph(*w).countcycles() >= 2) return WordGraph(*w);
    return WordGraph::zero();
  }

  WGSet star() const {
    WGSet xstar(one()), xnth(one());
    
    if(!empty()){
      do {
	xstar += xnth;
	xnth  *= *this;
	xnth  -= xstar;
      } while(!xnth.empty());
    }
    return xstar;
  }

  static WGSet zero(){ return WGSet(); }
  static WGSet one() { return WGSet(WordGraph::one()); }

  friend WGSet operator*(const WGSet& S, const WGSet& T)  { 
    return WGSet(S,T,WordGraph::plain_compose,WordGraph::synchronize_test); 
  }  
  
  WGSet& operator*=(const WGSet& T) {
    *this = (*this) * T;
    return *this;
  }

  friend WGSet operator+(const set<WordGraph> &X, const set<WordGraph> &Y)  {
    WGSet sum;
    set_union(X.begin(),X.end(),Y.begin(),Y.end(),inserter(sum,sum.begin()));
    return sum;
  }

  WGSet& operator+=(const set<WordGraph> &Y) {
    set<WordGraph>& X(*this);
    set_union(X.begin(),X.end(),Y.begin(),Y.end(),inserter(X,X.begin()));
    return *this;
  }

  WGSet& operator+=(const WordGraph &x) { insert(x); return *this;}

  friend WGSet operator-(const set<WordGraph> &X, const set<WordGraph> &Y)  {
    WGSet sum;
    set_difference(X.begin(),X.end(),Y.begin(),Y.end(),inserter(sum,sum.begin()));
    return sum;
  }

  WGSet& operator-=(const set<WordGraph> &Y) {
    set<WordGraph>& X(*this);
    set_difference(X.begin(),X.end(),Y.begin(),Y.end(),inserter(X,X.begin()));
    return *this;
  }


  friend ostream& operator<<(ostream& s, const WGSet& S)
  {
    s << "{"; 
    set<WordGraph>::const_iterator x(S.begin());
    if(x != S.end()){
      s << WordGraph(*x);
      for(++x; x != S.end(); x++) s << "," << WordGraph(*x);
    }
    s << "}";
    return s;
  }  
};

typedef ringmatrix<WGSet> WGMatrix;
typedef SRSet<word> wordset;
}


#endif
