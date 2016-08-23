#pragma once

#include <vector>
#include <set>

#include <auxiliary.hh>
#include <algebra/digraph.hh>
#include <algebra/ringmatrix.hh>

using namespace std;


/**********************************************************************

  NFAi is an NFA with every state both initial and accepting, and with
  no sink- or source-states.

 **********************************************************************/
class NFAi {
public:
  int Nstates, Nletters;
  /// step[node0,a+1] = {node1,node2,...,nodem} is set of nodes reachable in a single 'a'-step. epsilon is index a=0.
  vector< vector< set<int> > > step;	

  NFAi(int Nstates, int Nletters) : Nstates(Nstates), Nletters(Nletters), step(Nstates,vector< set<int> >(Nletters+1)) {}

  void insert_transition(int state0, int state1, int letter){ 
    step[state0][letter+1].insert(state1);
  }

  void insert_epsilon_transition(int state0, int state1){             
    step[state0][0].insert(state1);
  }

  void epsilon_closure(int state0, set<int>& closure) const;
  set<int> epsilon_closure(int state0) const;

  void Delta(const set<int> &P, int a, set<int> &Q) const;
  set<int> Delta(const set<int>& P, int a) const;

  bool is_deterministic() const;
  bool epsilon_free() const;
  NFAi minimize() const;
  NFAi reverse() const;
  NFAi dual() const;
  NFAi disjoint_union(const NFAi& B) const;

  digraph skeleton() const;
  map<digraph::dedge_t,vector<int> > label_map() const;
  ringmatrix<int> adjacency_matrix() const;

  void trim(); 
  NFAi trim() const { NFAi A(*this); A.trim(); return A; }


  bool equivalent_to(const NFAi& B) const;
  bool operator==(const NFAi& B) const { return equivalent_to(B); }
  vector<string> default_nodenames() const;
  vector<string> default_labels(bool pairs=true) const;
  string to_latex(vector<string> labels = vector<string>(), vector<string> nodenames = vector<string>()) const;

  u32string hash_string() const;
 
  size_t hash() const {
    return std::hash<u32string>()(hash_string());
  }

  friend ostream& operator<<(ostream& S, const NFAi &nfa);
};

