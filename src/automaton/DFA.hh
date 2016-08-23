#pragma once

#include <automaton/NFA.hh>
#include <algebra/digraph.hh>

/*
 * DFAi is an DFA with every state both initial and accepting, and with no sink- or source-states.
 */
class DFAi {
public:
  int Nstates, Nletters;
  vector< vector<int> > step;		// step[node0,letter] = 0|node1+1
  vector<string> alphabet;

  DFAi(int Nstates, int Nletters) : Nstates(Nstates), Nletters(Nletters), step(Nstates, vector<int>(0,Nletters)) {}
  DFAi(const DFAi& dfa, const vector<string>& alphabet) 
    : Nstates(dfa.Nstates), Nletters(dfa.Nletters), step(dfa.step), alphabet(alphabet) {}

  DFAi(const NFAi& nfa, const vector<string>& alphabet_ = vector<string>()) : Nletters(nfa.Nletters) {
    NFAi nfat = nfa.trim();
    if(nfat.epsilon_free())
      *this = determinizeB(nfat);
    else 
      *this = determinizeA(nfat);

    alphabet = alphabet_;
  }

  static DFAi determinizeA(const NFAi& nfa); 
  static DFAi determinizeB(const NFAi& nfa);

  NFAi NFA() const;     // NFA-representation of DFA
  NFAi reverse() const; // rev(A) = {q--a-->p | p--a-->q \in A}
  NFAi dual() const;    // A* = {a--q-->b | p--a-->q--b-->r \in A}
  void trim();          // Remove all sinks, sources, and dead states.


  digraph skeleton() const;

  DFAi minimize() const {
    return DFAi(DFAi(reverse()).reverse(),alphabet);
  }
  // TODO: Add transformation to canonical form (order nodes) to equivalence test as direct comparison.

  static DFAi minimize(const NFAi& nfa) {				
    // Skips one determinization. 
    return DFAi(DFAi(nfa.reverse()).reverse());
  }

  DFAi              reorder(const vector<int>& order) const;
  size_t            hash() const { 
    u32string h(Nstates*Nletters,0); 
    for(int p=0;p<Nstates;p++) 
      for(int a=0;a<Nletters;a++) h[p*Nletters+a] = step[p][a];
    return std::hash<u32string>()(h);
  
  }

  friend ostream& operator<<(ostream& S, const DFAi &dfa);
};

