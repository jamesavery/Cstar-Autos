#pragma once
#include <inttypes.h>
#include <vector>

#include <wordgraph/WordGraph.hh>
using namespace WordGraphLR;

struct permutation_graph_edge { node_t alpha, beta, e, f; };
typedef std::vector<permutation_graph_edge>                         permutation_graph;
typedef std::pair< std::vector<WordGraph>, std::vector<WordGraph> > permutation_graph_wg;

typedef enum { UNDEFINED, OK, NO_REMAINING_ARCS, SHORT_DC_LEFT, SHORT_DC_RIGHT, NONSYNC_LEFT, NONSYNC_RIGHT} arcstep_type;
extern const char* arcstep_txt[7];

struct permutation_metadata {
  ringmatrix<int>     E;  // E is a directed multigraph. We are constructing endomorphisms on C*(E)
  int k, m, L;		  // m = |E^0|, L = |E^1|, M = |E^0_\tau| = |E^{k-1}|, 

  ringmatrix<int>        Esum;        // Earcs(s,t)+j is the edge name of the j'th edge e:s->t in E
  vector<node_t> E_source, E_range;    // s(e)   and r(e) for arcs e in E^1 (result in E^0)
  vector<vector<node_t>> E_from, E_to; // E_from[v] = E_{v->*}
                                       // E_to[v]   = E_{*->v}
  vector<int> indegree,  outdegree;    // indegree[u] = |{(u,v) \in E}|, outdegree[v] = |{(u,v) \in E}|
  
  ringmatrix<int>     Etau0, Ek, Etau0sum;    // Etau0(s,t) = |E^{k-1}_{s->t}|, Ek(s,t) = |E^k_{s->t}|
  int M;

  // Given node v in original graph, Ns[v] = |E^{k-1}_{v->*}| and Nr[v] = |E^{k-1}_{*->v}|
  vector<int>  Ns, Nr;  
  vector<node_t> Etau_source, Etau_range;    // s(\mu) and r(\mu) for arcs \mu in E^1_\tau (result in E^{k-1})

  vector<vector<node_t>> Etau_from, Etau_to; // Etau_from[v] = E^{k-1}_{v->*}
                                             // Etau_to[v]   = E^{k-1}_{*->v}


  // Etau_from[s] = { \alpha_j | \alpha_j: s\to v for some v }
  // Etau_to[r]   = { \alpha_k | \alpha_k: v\to r for some v }
  // L1local[\alpha] = j such that \alpha = \alpha^j in Etau_from[source(\alpha)]
  // L2local[\alpha] = k such that \alpha = \alpha^k in Etau_to[range(\alpha)]
  //
  //   beta -e,f-> alpha
  //
  //   e: s->r, f:q->t
  //   beta: s->q
  //   alpha: r->t
  //
  //   beta\in  Etau_from[s], Etau_to[q]
  //   alpha\in Etau_from[r], Etau_to[t]
  
  vector<node_t> L1local, L2local; // global to local node id's // TODO: Rename, check proper definition
  vector<string> node_names;

  permutation_metadata(const ringmatrix<int> &E, int k=2);

  permutation_graph identity() const;

  void     initialize_letters(permutation_graph_wg& letters) const;
  WGMatrix letter_matrix(const vector<WordGraph>& letters) const;

  void                print_paths(ostream &S, int k) const;
  vector<string>      path_names(int k)  const;
  
  ringmatrix<wordset> letter_strings()   const;
  
};


struct partial_permutation_graph { 
  const permutation_metadata &P;
  permutation_graph    graph;
  permutation_graph_wg graph_wg;

  std::vector< vector<int> > remaining;
  std::vector<int> lowest_unused;   // lowest_unused[u*m+v]   is smallest $E_\tau$-node u->v with out-degree 0.
  std::vector<int> lowest_unused_f; // lowest_unused_f[u*m+v] is lowest E-arc  u->v that has not already been used as an f

  partial_permutation_graph(const permutation_metadata &P);

  bool is_synchronizing(arcstep_type &reason) const;
  void push_edge(const permutation_graph_edge& mu);
  void pop_edge(const permutation_graph_edge& mu);
};



bool operator <(const permutation_graph_edge& e1,const permutation_graph_edge& e2);
bool operator==(const permutation_graph_edge& e1,const permutation_graph_edge& e2);
ostream& operator<<(ostream& stream,const permutation_graph_edge& e);
// TODO: Organization




