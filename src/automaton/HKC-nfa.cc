#include "HKC-nfa.hh"
#include "NFA.hh"
#include "DFA.hh"
#include <string>
#include <stdlib.h>
#include <cassert>
#include <set>
#include <vector>
#include <map>
#include <algorithm>

using namespace std;

//Finds the set of all states that can reached from the 
//set p through the letter a.
set<int> Delta(const NFAi& N, const set<int>& p, int a) {
  set<int> states;
  for(set<int>::iterator it(p.begin()); it != p.end(); ++it) {
    states.insert(N.step[*it][a].begin(), N.step[*it][a].end());
  }
  return states;
}

bool check_equal_HKC(const NFAi& A,const NFAi& B) {
  NFAi N = combine(A, B);
  set<int> IA, IB; //Initial states in A and B
  for(int i = 0;i < A.Nstates; ++i) IA.insert(IA.end(), i);
  for(int i = 0;i < B.Nstates; ++i) IB.insert(IB.end(), A.Nstates+i);

  return HKC_NFA(N, IA, IB);
}

//Returns true iff X \subset Y
bool subset(const set<int>& X, const set<int>& Y) {
  return includes(Y.begin(),Y.end(), X.begin(), X.end());
  // set<int> Z(Y.begin(), Y.end());
  // Z.insert(X.begin(), X.end());
  // return Z.size() == Y.size();
}

bool HKC_NFA(const NFAi& N, const set<int>& X,const set<int>& Y) {
  deque< pair< set<int>, set<int> > > S;
  Relation R;
    
  S.push_back(make_pair(X, Y));
    
  int i = 0;
  while(!S.empty()) {
    ++i;
    //    printf("%d: |S| = %d\n",i,S.size());
    const pair< set<int>, set<int> > &pq(S.front()); 
    const set<int> p(pq.first), q(pq.second);
    S.pop_front();

    if(R.getNormalForm(p) == R.getNormalForm(q)) continue;
    if(p.empty() != q.empty()) return false;

    for(int a = 0; a < N.Nletters; ++a) {
      set<int> p2(Delta(N, p, a)), q2(Delta(N, q, a));
      S.push_back(make_pair(p2, q2));
    }
    R.updateRelation(p, q);
  }
  cout << "Iterations: " << i << endl;
  return true;
}

bool HK_NFA(const NFAi& A,const NFAi& B) {
  const NFAi& N(combine(A, B));
  UnionFind UF;
  queue< pair< set<int>, set<int> > > S;

  set<int> IA, IB; //Initial states in A and B
  for(int i = 0;i < A.Nstates; ++i) IA.insert(IA.end(), i);
  for(int i = 0;i < B.Nstates; ++i) IB.insert(IB.end(), A.Nstates+i);
  UF.make(IA); //line 1
  UF.make(IB); //line 2
  UF.unite(IA,IB); //line 4
  S.push(make_pair(IA, IB)); //line 5

  int i = 0, j = 0;
  while(!S.empty()) { //line 6
    ++i;
    const pair< set<int>, set<int> > &pq = S.front();
    set<int> p = pq.first, q = pq.second;
    S.pop();

    //this is equivalent to eps(p) != eps(q) since all states are final.
    if(p.empty() != q.empty()) return false; //line 7+8
    for(int a = 0; a < N.Nletters; ++a) { //line 9
      int p2Id = UF.find(Delta(N, p, a)), q2Id = UF.find(Delta(N, q, a)); //line 10+11
      if(p2Id != q2Id) { //line 12 
	set<int> p2 = UF.intToSet(p2Id), q2 = UF.intToSet(q2Id);
	UF.unite(p2, q2); //line 13
	S.push(make_pair(p2, q2)); //line 14
      } else { 
	// printf("p2Id = %d, q2Id = %d\n",p2Id,q2Id);
	// cout << "Delta(N,p,a) = " << Delta(N,p,a) << endl;
	// cout << "Delta(N,q,a) = " << Delta(N,q,a) << endl;
	++j; 
      }
    }
  }
  cout << "Iterations: " << i + j << endl;
  return true; //line 15
}

NFAi combine(const NFAi& A,const NFAi& B) {
  assert(A.Nletters == B.Nletters); //assume same alphabet
  NFAi N(A.Nstates+B.Nstates, A.Nletters);
  for(size_t n = 0; n < A.step.size(); ++n) {
    for(size_t a = 0; a < A.step[n].size(); ++a) {
      for(set<int>::iterator it(A.step[n][a].begin()); it != A.step[n][a].end(); ++it) {
	if(a == 0) {
	  N.insert_epsilon(n, *it);
	} else {
	  N.insert_transition(n, *it, a-1);
	}
      }
    }
  }
  for(size_t n = 0; n < B.step.size(); ++n) {
    for(size_t a = 0; a < B.step[n].size(); ++a) {
      for(set<int>::iterator it(B.step[n][a].begin()); it != B.step[n][a].end(); ++it) {
	if(a == 0) {
	  N.insert_epsilon(A.Nstates+n, A.Nstates+(*it));
	} else {
	  N.insert_transition(A.Nstates+n, A.Nstates+(*it), a-1);
	}
      }
    }
  }
  return N;
}
