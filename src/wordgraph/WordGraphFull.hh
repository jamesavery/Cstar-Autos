#ifndef WORDGRAPHFULL_HH
# define WORDGRAPHFULL_HH

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

namespace WordGraphFull {
  typedef int node_t;
  using namespace std;
  
  class WordGraph {		// TODO: WordGraph: public digraph
  public:
    digraph A;

    typedef enum {ZERO,ONE,NORMAL} element_type_t;
    element_type_t e_type; 

    element_type_t element_type() const { return e_type; }

    WordGraph(const WordGraph& Y) : A(Y.A), e_type(Y.e_type) {    }
    /*
    WordGraph(const int length = 0, const element_type_t e_type = NORMAL) : A(length), e_type(e_type) {}
    WordGraph(const vector<node_t>& data_, const element_type_t e_type = NORMAL) : A(data_.size()), e_type(e_type) {
      for(int i=0;i<data_.size();i++) A[data_[i]].push_back(i);
    }
    WordGraph(const node_t value, const int length) : A(length), e_type(NORMAL) {}
    */		  
    WordGraph(const digraph& A = digraph(0,0), element_type_t e_type = NORMAL) : A(A), e_type(e_type) {}

    WordGraph(const vector<bool>& in, const vector<bool>& out) : A(in.size(),out.size()), e_type(NORMAL)
    {
      for(int i=0;i<in.size();i++)
	if(in[i]) for(int j=0;j<out.size();j++)
		    if(out[j]) A[i].push_back(j);
    }

    WordGraph(int n) : A(0,0) {		// Ring elements
      assert(n==0||n==1);
      *this = (n==0? zero() : one());
    }
    
    static WordGraph one() { return WordGraph(digraph(0,0),ONE);  }
    static WordGraph zero(){ return WordGraph(digraph(0,0),ZERO); }

    WordGraph normalize(int n) const;

    // i runs from 1,...,n, value from 1:
    // 0: no arc, 1,..,15: actual indices
    void set(int j, int i){  
      vector<int>& row(A[i-1]);
      if(find(row.begin(),row.end(),j-1) == row.end())
	row.push_back(j-1); 
    } 

  
    //    const vector<int>& get(uint8_t i) const   { return A[i-1]; }
    uint8_t length() const         { return A.N; }
    void setdimension(int rows, int cols){ 
      A.M  = rows; A.N = cols;
      A.resize(rows); 
    }
    //    void setlength(int length)     { data.resize(length+1); }

    bool operator==(const WordGraph& y) const { return e_type == y.e_type && A == y.A; }
    bool operator>(const WordGraph& y)  const { return e_type > y.e_type  || (e_type == y.e_type && A  > y.A); }
    bool operator<=(const WordGraph& y) const { return e_type <= y.e_type || (e_type == y.e_type && A <= y.A); }
    bool operator<(const WordGraph& y)  const { return e_type < y.e_type  || (e_type == y.e_type && A  < y.A); }

    WordGraph operator*(const WordGraph& g) const;
    WordGraph operator&(const WordGraph& g) const;

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

    vector<bool> innodes() const;
    vector<bool> outnodes() const;

    // TODO: Clean up
    WordGraph transpose() const;
    static WordGraph connected_arcs(const WordGraph& U, const WordGraph& A);
    static bool has_parallel_words(const WordGraph& A, const WordGraph& B);

    static bool is_anchored(const WordGraph& U, const WordGraph& a, const WordGraph V);
    static bool left_anchored(const WordGraph& U, const WordGraph& a);
    static bool right_anchored(const WordGraph& a, const WordGraph& V);

    static bool left_anchored2(const WordGraph& U, const WordGraph& a);
    static bool right_anchored2(const WordGraph& a, const WordGraph& V);

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
