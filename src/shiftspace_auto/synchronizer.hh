#pragma once

#include <shiftspace_auto/permutation_dual.hh>
#include <algebra/digraph.hh>
#include <wordgraph/WordGraph.hh>

/* 
   s--beta---p
   |         |
   e         f
   |         |
   r--alpha--q
 */

class Synchronizer {
public:
  bool LR;
  bool fits64bit;

  vector< set<digraph::dedge_t> > WGEdges;
  

  const permutation_metadata& P;
  
  const bool left;
  vector< IDCounter<int> > Psources, Pranges;
  vector<int> local_name;
  int M;
  NFAi AL, AR;

  typedef digraph::dedge_t dedge_t;
  Synchronizer(const permutation_metadata &P, const permutation_dual &phi, bool left);

  bool synchronizes() const;
  template <typename WGSet> bool synchronizes(const ringmatrix<WGSet>& A) const;

  template <typename WG> vector<WG> get_letters() const;
  bool degenerate() const { return !AL.equivalent_to(AR); }

  template <typename WGSet> ringmatrix<WGSet> letter_matrix() const;
};

class DualSynchronizer {
public:
  const permutation_metadata& P;
  const permutation_dual &phi;
  const permutation_metadata  E;
  bool left_resolving, right_resolving;
  bool fits64bit;


  int dL; 

  typedef int node_t;
  vector<int> sources, ranges;

  vector< set<digraph::dedge_t> > WGEdgesL, WGEdgesR;
  NFAi AL, AR;
  ringmatrix<WordGraphFull::WGSet> DL,DR;

  typedef digraph::dedge_t dedge_t;
  
  // TODO: Should LR also be for dual? I think so.
  DualSynchronizer(const permutation_metadata& P, const permutation_dual& phi0);

  bool degenerate() const {
    return !AL.equivalent_to(AR);
  }

  bool synchronizes(bool left) const {
    return synchronizes<WordGraphFull::WGSet>(left?DL:DR);
  }

  template <typename WGSet> bool synchronizes(const ringmatrix<WGSet>& A) const;
};


