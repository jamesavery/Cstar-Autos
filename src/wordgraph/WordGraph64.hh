#ifndef WORDGRAPH64_HH
# define WORDGRAPH64_HH

#include <map>
#include <vector>
#include <set>
#include <algorithm>
#include <stdint.h>
#include <limits.h>
#include <algebra/SGPair.hh>
#include <algebra/SRSet.hh>
#include <algebra/ringmatrix.hh>


namespace WordGraph64 {
  typedef uint16_t node_t;
  using namespace std;
  
  class WordGraph {
  public:
    // TODO: Figure out what is a better one-value, ~0xff (in band) or -1 (out of band).
    enum { zero_value = 0L, one_value = ~0xffL}; 
    static const int MaxVertices = 15;
    uint64_t   data;

    template<class Archive> void serialize(Archive & archive){ archive(data); }


    typedef enum {ZERO,ONE,NORMAL} element_type_t;

    element_type_t element_type() const { return data == zero_value? ZERO : (data == one_value? ONE : NORMAL); }
 
    WordGraph(const uint64_t& data = 0, const uint8_t& length = MaxVertices) : data((data & ~0xfL) | length) {}

    WordGraph(const vector<node_t>& indices) : data(indices.size()) {
      for(unsigned i=0;i<indices.size();i++) set(i+1,indices[i]);
    }
    WordGraph(const vector<int>& indices): data(indices.size()) {
      for(unsigned i=0;i<indices.size();i++) set(i+1,indices[i]);
    }
    WordGraph(int n){		// Ring elements
      assert(n==0||n==1);
      *this = (n==0? zero() : one());
    }

    static WordGraph from_int(const int64_t& data){
      WordGraph w;
      w.data = data;
      return w;
    }
    
    static WordGraph one() { return from_int(one_value); }
    static WordGraph zero(){ return from_int(zero_value); }

    // i runs from 1,...,15, value from 1: 0: no arc, 1,..,15: actual indices
    inline void set(uint8_t i, uint8_t value){  
      data &= ~(0xfL << (i*4));  
      data |= (value&0xfL) << (i*4);
    }

    // i runs from 0,...,14, value from 0,...,14 -- can't set "no arc"
    inline void setb(uint8_t i, uint8_t value){  
      data &= ~(0xfL << ((i+1)*4));  
      data |= ((value+1)&0xfL) << ((i+1)*4);
    }

  
    inline uint8_t get(uint8_t i) const   {  return ((data&~0xfL) >> i*4) & 0xfL;   }
    inline uint8_t length() const         { return data&0xfL; }
    inline void setlength(int length)     { 
      data &= ~0xfL;
      data |= (length&0xfL); 
    }
    inline void setdimension(int rows, int cols){ setlength(cols); }

    inline bool operator==(const WordGraph& y) const { return data == y.data; }
    inline bool operator>(const WordGraph& y)  const { return data  > y.data; }
    inline bool operator<=(const WordGraph& y) const { return data <= y.data; }
    inline bool operator<(const WordGraph& y)  const { return data < y.data; }

    WordGraph operator*(const WordGraph& g) const;

    int countarcs(int max=INT_MAX) const;
    int countcycles(int max=INT_MAX) const;
    bool synchronizes() const;

    uint64_t number() const;
    static WordGraph fromnumber();


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


class WGSet : public set<uint64_t> {
public:
  typedef WordGraph wordgraph_t;
  typedef uint64_t  value_t;
  typedef set<uint64_t> intset;
  typedef WordGraph binop_t(const WordGraph& x, const WordGraph& y);
  typedef bool predicate_t(const WordGraph& x);

  WGSet(const intset& S = intset()) : intset(S.begin(),S.end()) { }
  WGSet(const WordGraph& x) { insert(x.number()); } 
  WGSet(const WGSet &X, const WGSet &Y, binop_t compose, predicate_t include) 
  {
    for(const auto &x: X)
      for(const auto &y: Y){
	const WordGraph z(compose(WordGraph::from_int(x),WordGraph::from_int(y)));
	//	cout << "product " << WordGraph(*x) << " * " << WordGraph(*y) << " = " << z << endl;
	if(include(z)){
	  //	  cout << z << " passes include test.\n";
	  insert(z.number());
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

  friend WGSet operator+(const set<uint64_t> &X, const set<uint64_t> &Y)  {
    WGSet sum;
    set_union(X.begin(),X.end(),Y.begin(),Y.end(),inserter(sum,sum.begin()));
    return sum;
  }

  WGSet& operator+=(const set<uint64_t> &Y) {
    set<uint64_t>& X(*this);
    set_union(X.begin(),X.end(),Y.begin(),Y.end(),inserter(X,X.begin()));
    return *this;
  }

  WGSet& operator+=(const WordGraph &x) { insert(x.number()); return *this; }

  friend WGSet operator-(const set<uint64_t> &X, const set<uint64_t> &Y)  {
    WGSet sum;
    set_difference(X.begin(),X.end(),Y.begin(),Y.end(),inserter(sum,sum.begin()));
    return sum;
  }

  WGSet& operator-=(const set<uint64_t> &Y) {
    set<uint64_t>& X(*this);
    set_difference(X.begin(),X.end(),Y.begin(),Y.end(),inserter(X,X.begin()));
    return *this;
  }


  friend ostream& operator<<(ostream& s, const WGSet& S)
  {
    s << "{"; 
    set<uint64_t>::const_iterator x(S.begin());
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
