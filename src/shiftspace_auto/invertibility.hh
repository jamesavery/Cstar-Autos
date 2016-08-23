#pragma once

#include <algebra/digraph.hh>
#include <wordgraph/WordGraphFull.hh>

using namespace WordGraphFull;

typedef WordGraphFull::WGSet WGSet;
typedef WordGraphFull::WordGraph WordGraph; // TODO: Multiplex som i synchronizer.

/**********************************************************************

  Implementation of threaded concatenation operator through triplets 
  of wordgraphs.

 **********************************************************************/
struct triplet {
  WordGraph U, A, V;

  triplet() {}
  triplet(const WordGraph& U, const WordGraph& A, const WordGraph& V) : U(U.normalize(A.A.M)), A(A), V(V.normalize(A.A.N)) {}

  bool    right_anchored() const { return WordGraph::left_anchored(U,A);  }
  bool    left_anchored()  const { return WordGraph::right_anchored(A,V); }
  bool    is_anchored()    const { return WordGraph::is_anchored(U,A,V);  }
  bool    has_parallel_words() const;

  triplet prune() const;
  bool operator<(const triplet& t) const;
  
  string to_latex(const WordGraph& L=WordGraph::one(), const WordGraph& R=WordGraph::one()) const;
};


inline triplet operator*(const triplet& T, const WordGraph& W){
  return triplet(T.U,T.A,T.V*W).prune();
}

inline triplet operator*(const WordGraph& W, const triplet& T){
  return triplet(W*T.U,T.A,T.V).prune();
}

/**********************************************************************

              Implementation of invertibility algorithm

 **********************************************************************/
class Invertibility {
public:
  const digraph G;
  map<digraph::dedge_t, string> labels;
  int L;

  vector<WordGraph> letters;

  // TODO: Optional map label -> domain (source -> range). When the
  // label is a graph homomorphism, we have nice things.
  Invertibility(int N, const map<digraph::dedge_t, vector<int>> &labels_);


  bool invertible() const {
    set<triplet> S, Snext;

    for (int a = 0; a < letters.size(); a++) {
      const triplet t =
          triplet(WordGraph::one(), letters[a], WordGraph::one()).prune();
      if (!t.is_anchored())
        S.insert(t);
    }

    for (int i = 0; !S.empty(); i++) {
      bool parallel_or_doublecycle = breadth_first(S, Snext, i & 1);
      if (parallel_or_doublecycle)
        return false;
      S = Snext;
    }
    return true;
  }

  bool breadth_first(const set<triplet> &workset, set<triplet> &next,
                     const bool even = false) const {
    next.clear();
    for (const auto &e: workset) {
      const triplet &x(e);

      for (const auto &a: letters) {
        triplet t;

        if (even) {
          t = x * a;
          if (WordGraph::has_parallel_words(x.V, a))
            return true; // PW of length n+1 found
        } else {
          t = a * x;
          if (WordGraph::has_parallel_words(a, x.U))
            return true; // PW of length m+1 found
        }

        if (t.has_parallel_words())
          return true; // PW of length m+n+1 found

        if ((t.U * t.A * t.V).countcycles() > 1)
          return true; // DC of length m+n+1 found

        if (!t.is_anchored())
          next.insert(t);
      }
    }
    return false;
  }


  // Step-by-step output of algorithm to LaTeX.
  bool invertible_demo() const;
  bool breadth_first_demo(const set<triplet> &workset, set<triplet> &next,
                          const bool even, const string name) const;

  string labeled_graph_latex(const string name) const;
  static string triplet_set_latex(const string name, const set<triplet> &S);
  static string triplet_vector_latex(const string name,
                                     const vector<triplet> &S);
  string triplet_set_product_latex(const string name, const set<triplet> &S,
                                   bool left) const;
  static string WG_latex(const WordGraph &W);
  static string WGprod_latex(const WordGraph &L, const WordGraph &R);
};


template <typename WGSet> class InvertibleLR {
  typedef ringmatrix<WGSet> WGMatrix;
  int M;

  bool is_synchronizing(const WGMatrix &A)  const
  {
    //    if(synchronizes_quickly(A,3)) return true;
    WGMatrix power(A), sum(A.m,A.m);

    for(int k=1;k<M*(M-1)+1;k++){
      unsigned int oldsize = sum.total_size();
      sum += power;

      if(sum.total_size() == oldsize)
	return true;

      power = A*power;
      for(int i=0;i<A.m;i++)
	if(power.getelement(i,i).doublecycle().element_type() != WordGraph::ZERO){
	  // Length-(k+1) word i->i with double-cycle.
	  return false;
	}
    }
    return false;	
  }

};
