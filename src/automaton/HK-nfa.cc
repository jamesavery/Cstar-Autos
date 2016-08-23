#include <string>
#include <stdlib.h>
#include <cassert>
#include <set>
#include <vector>
#include <map>
#include <algorithm>
#include "NFA.hh"
#include "DFA.hh"

using namespace std;

class UnionFind {
public:
  int n; //number of sets
  vector<int> p; //parent pointer
  vector<int> l; //max-length of subtree
  map<set<int>, int>  m; //mapping sets to their number.
  vector<set<int> > setList; //List of all sets
	
  UnionFind(): n(0), p(0) {}
  void make(const set<int>& s) {
    setToInt(s);
  }
  void unite(const set<int>& s1, const set<int>& s2) {
    unite(setToInt(s1), setToInt(s2));
  }
  int find(const set<int>& s) {
    find(setToInt(s));
  }
  set<int> intToSet(int a) {
    return setList[a];
  }
  int setToInt(const set<int>& s) {
    map<set<int>, int>::iterator it = m.find(s);
    if(it == m.end()) {
      ++n;
      p.push_back(n-1);
      l.push_back(0);
      m[s] = n-1;
      setList.push_back(s);
      return n-1;
    } else {
      return it->second;
    }
  }
  int find(int a) {
    if(p[a] != a) {
      p[a] = find(p[a]);
    }
    return p[a];
  }
  void unite(int a, int b) {
    if(l[a] <= l[b]) {
      p[a] = b;
      l[b] = max(a+1,b);
    } else {
      p[b] = a;
      l[a] = max(a,b+1);
    }
  }
};


NFAi combine(NFAi A, NFAi B);

//Finds the set of all states that can reached from the 
//set p through the letter a.
set<int> Delta(const NFAi& N, set<int>& p, int a) {
  set<int> states;
  for(set<int>::iterator it(p.begin()); it != p.end(); ++it) {
    states.insert(N.step[*it][a].begin(), N.step[*it][a].end());
  }
  return states;
}

bool HK_NFA(NFAi A, NFAi B) {
  NFAi N = combine(A, B);
  UnionFind UF;
  queue< pair< set<int>, set<int> > > S;

  set<int> IA, IB; //Initial states in A and B
  for(int i = 0;i < A.Nstates; ++i) IA.insert(i);
  for(int i = 0;i < B.Nstates; ++i) IB.insert(A.Nstates+i);
  UF.make(IA); //line 1
  UF.make(IB); //line 2
  UF.unite(IA,IB); //line 4
  S.push(make_pair(IA, IB)); //line 5

  while(!S.empty()) { //line 6
    pair< set<int>, set<int> > pq = S.front();
    S.pop();
    set<int> p = pq.first, q = pq.second;
    //this is equivalent to eps(p) != eps(q) since all states are final.
    if(p.empty() != q.empty()) return false; //line 7+8
    for(int a = 0; a < N.Nletters; ++a) { //line 9
      int p2Id = UF.find(Delta(N, p, a)), q2Id = UF.find(Delta(N, q, a)); //line 10+11
      if(p2Id != q2Id) { //line 12 
	set<int> p2 = UF.intToSet(p2Id), q2 = UF.intToSet(q2Id);
	UF.unite(p2, q2); //line 13
	S.push(make_pair(p2, q2)); //line 14
      }
    }
  }
  return true; //line 15
}

NFAi combine(NFAi A, NFAi B) {
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
