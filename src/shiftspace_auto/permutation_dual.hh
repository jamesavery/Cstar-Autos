#pragma once

#include <iostream>
#include <auxiliary.hh>
#include <algebra/SGPair.hh>
#include <algebra/SRSet.hh>
#include <automaton/NFA.hh>
#include <automaton/DFA.hh>
#include <shiftspace_auto/permutation_graph.hh>

typedef SGPair< SGString<int>, SGString<int> > PDualG;
typedef SRSet<PDualG> PDualR;

// A[e,f] = (beta,alpha) <=> e -(beta,alpha)-> f \in osv.
// (AB)[e,g] \ni  A[e,f]*B[f,g] = e-(b1,a1)->f-(b2,a2)->g = e-(b1 b2, a1 a2)->g

class permutation_dual: public ringmatrix< PDualR > {
public:

  permutation_dual(int n=1) : ringmatrix<PDualR>(n,n) {}
  permutation_dual(const ringmatrix<PDualR>& Y) : ringmatrix<PDualR>(Y) {}
  permutation_dual(const permutation_metadata& P, const permutation_graph& Etau): ringmatrix<PDualR>(P.L,P.L)
  {
    for(auto edge: Etau){
      PDualR &Aef((*this)(edge.e,edge.f));
      Aef += PDualG(edge.beta,edge.alpha);
    }
  }

  permutation_dual(const DFAi& dfa) : ringmatrix<PDualR>(int(round(sqrt(dfa.Nletters))),int(round(sqrt(dfa.Nletters))))
  {
    assert(n*n == dfa.Nletters);
    
    for(int p=0;p<dfa.Nstates;p++)
      for(int e=0;e<n;e++)
	for(int f=0;f<n;f++){
	  int q = dfa.step[p][e*n+f];

	  if(q!=0)
	    (*this)(e,f) += PDualG(p,q-1);
	}
  }

  permutation_dual(const NFAi& nfa) : ringmatrix<PDualR>(int(round(sqrt(nfa.Nletters))),int(round(sqrt(nfa.Nletters))))
  {
    assert(n*n == nfa.Nletters);

    for(int p=0;p<nfa.Nstates;p++)
      for(int e=0;e<n;e++)
	for(int f=0;f<n;f++){
	  const set<int> &Q(nfa.step[p][e*n+f+1]);

	  for(set<int>::const_iterator q(Q.begin()); q!=Q.end();q++)
	    (*this)(e,f) += PDualG(p,*q);
	}
  }

  permutation_dual inverse() const {
    permutation_dual invG(n);

    for(int i=0;i<n;i++)
      for(int j=0;j<n;j++){
	const PDualR &Gij((*this)(i,j));
	for(auto g: Gij) invG(i,j) += PDualG(g.second,g.first);
      }
    return invG;
  }

  void power(permutation_dual& A, const unsigned int n, const size_t max_size=INT_MAX) const {
    if(A.Nstates() > max_size){
      // TODO: Return value that tells whether operation was successful
      fprintf(stderr,"ERROR: Permutation graph textile system size %d larger than %ld, aborting %d'th power.\n",A.Nstates(), max_size, n);
      return;
    }
    if(n==0){ A = one(m); return; }
    if(n==1){ A = *this; return; }
    if(n&1){
      power(A,n-1);
      A = (*this)*A;
    } else {
      power(A,n>>1);
      A = A*A;
    }
    A = A.minimize();
  }


  permutation_dual power(const unsigned int N, const size_t max_size=INT_MAX) const 
  {
    if(N==0) { 
      //      std::cout << "A = " << *this << std::endl;
      //      std::cout << "A^0 = " << one(m) << std::endl; 
      return one(m);
    }
    if(N==1) return *this;

    permutation_dual A(one(m));
    power(A,N,max_size);
    return A;
  }

  static permutation_dual identity(const ringmatrix<int>& O)
  {
    // Cleaner to construct directly from O, but this is simpler to write.
    permutation_dual sigma(shift(O,1));
    return permutation_dual(sigma*sigma.transpose()).minimize();
  }

  static permutation_dual shift(const ringmatrix<int>& E, int k) {
    if(k==0) return identity(E);
    if(k<0)  return shift(E,-k).transpose();

    permutation_metadata P(E,k+1);
    permutation_graph Id(P.identity());
    return permutation_dual(P,Id);
  }

  bool degenerate() const {
    NFAi left(automaton_component(true)), right(automaton_component(false));
    return !left.equivalent_to(right);
  }

  permutation_dual flatten() const {
    const permutation_dual& A(*this);
    permutation_dual B(n);

    IDCounter< SGString<int> > states;

    // Assign numbers to the composite states
    for(int e=0;e<n;e++)
      for(int f=0;f<n;f++){
	const PDualR &Aef(A(e,f));
	for(PDualR::const_iterator ab(Aef.begin()); ab!=Aef.end(); ab++){
	  states.insert(ab->first);
	  states.insert(ab->second);
	}
      }

    // Rename states
    for(int e=0;e<n;e++)
      for(int f=0;f<n;f++){
	const PDualR &Aef(A(e,f));
	PDualR& Bef(B(e,f));

	for(PDualR::const_iterator ab(Aef.begin()); ab!=Aef.end(); ab++)
	  Bef.insert(PDualG(states(ab->first), states(ab->second)));
      }

    return B;
  }

  int Nstates() const {
    int nstates = 1;
    for(int e=0;e<n;e++)
      for(int f=0;f<n;f++){
	const PDualR &Aef((*this)(e,f));
	for(PDualR::const_iterator ab(Aef.begin()); ab!=Aef.end(); ab++){
	  SGString<int> beta_string(ab->first), alpha_string(ab->second);

	  for(SGString<int>::const_iterator b(beta_string.begin()); b!=beta_string.end(); b++)
	    nstates = max(*b+1,nstates);

	  for(SGString<int>::const_iterator a(alpha_string.begin()); a!=alpha_string.end(); a++)
	    nstates = max(*a+1,nstates);
	}
      }
    return nstates;
  }

  NFAi automaton() const {
    permutation_dual Ad(flatten());
    NFAi A(Ad.Nstates(),n*n);
    
    for(int e=0;e<n;e++)
      for(int f=0;f<n;f++){
	const PDualR &Aef(Ad(e,f));

	for(PDualR::const_iterator ab(Aef.begin()); ab!=Aef.end(); ab++) if(!ab->first.empty() && !ab->second.empty()){
	    int beta_global = *ab->first.begin(),
	       alpha_global = *ab->second.begin();
	    A.insert_transition(beta_global,alpha_global,e*n+f);
	  }
      }
    A.trim();
    return A;
  }

  NFAi automaton_component(const bool left) const {
    permutation_dual Ad(flatten());
    NFAi A(Ad.Nstates(),n);
    
    for(int e=0;e<n;e++)
      for(int f=0;f<n;f++){
	const PDualR &Aef(Ad(e,f));

	for(PDualR::const_iterator ab(Aef.begin()); ab!=Aef.end(); ab++) if(!ab->first.empty() && !ab->second.empty()){
	    int beta_global = *ab->first.begin(),
	       alpha_global = *ab->second.begin();
	    A.insert_transition(beta_global,alpha_global,left?e:f);
	  }
      }
    A.trim(); 
    return A;
  }



  bool is_diagonal() const {
    NFAi mA(automaton());

    for(int s=0;s<mA.Nstates;s++){
      for(int e=0;e<n;e++)
	for(int f=0;f<n;f++)
	  if(e!=f && !mA.step[s][e*n+f+1].empty()) return false;
    }
    return true;
  }

  static NFAi identity(const permutation_metadata& P) 
  {
    NFAi I(P.m,P.L*P.L);
    
    for(int s=0,l=0;s<P.m;s++)
      for(int r=0;r<P.m;r++)
	for(int i=0;i<P.E(s,r);i++,l++)
	  I.insert_transition(s,r,l*(P.L+1));

    return I;
  }

  bool is_identity(const permutation_metadata& P) const {
    NFAi I(identity(P)), A(automaton());
    return A.equivalent_to(I);
  }
 
  permutation_dual minimize() const {
    DFAi mA(DFAi::minimize(automaton()));
    return permutation_dual(mA);
  }

  bool operator==(const permutation_dual& B) const {
    return automaton().equivalent_to(B.automaton());
  }

  size_t hash() const { return automaton().hash(); }
  
  // Right resolving:
  // Does every node emit at most one f? Equivalent: Is right automaton deterministic?
  // Left resolving:
  // Does every node receive at most one e? Equivalent: Is reverse of left automaton deterministic?
  bool resolving(bool left) const {
    if(left) return automaton_component(true).reverse().is_deterministic();
    else return automaton_component(false).is_deterministic();
  }

  friend ostream &operator<<(ostream& S, const permutation_dual& A)
  {
    for(int e=0;e<A.n;e++)
      for(int f=0;f<A.n;f++){
	vector<PDualG> pairs(vectorize(A(e,f)));
	S << "(" << char(e+'a') << "," << char(f+'a') << ") -> " << pairs << endl;
      }
    return S;
  }
};



