#pragma once

#include <vector>
#include <deque>
#include <algorithm>
#include <auxiliary.hh>

using namespace std;

/// Simple implementation of sparse unlabeled directed graph.
class digraph : public vector< vector<int> > { 
public:
  typedef int node_t;
  typedef pair<node_t,node_t> dedge_t;
  int M,N;

  /// Construct directed bipartite MxN graph from sparse-row matrix
  digraph(int M, int N, const vector< vector<node_t> >& edges = vector<vector<node_t> >());

  /// Construct directed bipartite MxN graph from list of directed edges
  digraph(int M, int N, const vector<dedge_t>& edges);						

  /// Reverse all arcs in graph / transpose adjacency matrix.
  digraph transpose() const;

  /// Composition / sparse matrix multiplication
  digraph operator*(const digraph& B) const;
  
  
  /// Is the digraph fan-in free? (I.e., in-degree of each vertex is at most one).
  bool fan_in_free() const;
  /// Is the digraph fan-out free? (I.e., out-degree of each vertex is at most one).
  bool fan_out_free() const;

  /// Compute strongly connected components using Tarjan's linear-time algorithm
  vector< vector<node_t> > TarjanSCC() const;
  void TarjanSCC(node_t v, deque<node_t>& S, vector<int>& SCCindex, vector<int>& lowlink, int& index, vector<vector<node_t> >& SCCs) const;

  /// Insert a directed egde (idempotent).
  void insert(const dedge_t& e); 

  /// Representation as list of directed edges.
  vector<dedge_t> edge_list() const;

  string to_latex() const;
};


