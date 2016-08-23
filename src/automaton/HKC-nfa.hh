#ifndef HKC_NFA_HH
#define HKC_NFA_HH
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


NFAi combine(const NFAi& A,const NFAi& B);
set<int> Delta(const NFAi& N, const set<int>& p, int a);
bool HKC_NFA(const NFAi& N,const set<int>& X,const set<int>& Y);

bool check_equal_HKC(const NFAi& A, const NFAi& B);

bool subset(const set<int>& X, const set<int>& Y);

bool HK_NFA(const NFAi& A, const NFAi& B);


class Relation {
public:
  vector<pair<set<int>, set<int> > > m; //Contains (X, X's normal form) w.r.t. current relation
  Relation() {}
  // void updateRelation( set<int> X,  set<int> Y) {
  //   set<int> Z(X.begin(), X.end());
  //   set_union(Z.begin(),Z.end(),Y.begin(),Y.end(),inserter(Z,Z.end()));
  //   Z = getNormalForm(Z);
  //   for(size_t i = 0;i < m.size(); ++i) {
  //     if(subset(m[i].first, Z)  m[i].second.size() < Z.size()) {
  // 	m[i].second = Z;
  //     }
  //   }
  //   m.push_back(make_pair(X, Z));
  //   m.push_back(make_pair(Y, Z));
  // }
  void updateRelation(const set<int>& X, const set<int>& Y) {
    m.push_back(make_pair(X, Y));
    m.push_back(make_pair(Y, X));
  }

  set<int> getNormalForm(const set<int>& X0) const {
    set<int> X(X0);
    size_t oldSize = 0;
    while(oldSize < X.size()) {
      oldSize = X.size();
      for(size_t i = 0;i < m.size(); ++i) {
	pair<set<int>, set<int> > ab = m[i];
	set<int> a = ab.first, b = ab.second;
	if(subset(a, X)) 
	  set_union(X.begin(),X.end(),b.begin(), b.end(),inserter(X,X.end()));
      }
    }
    return X;
  }
};

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
    return find(setToInt(s));
  }
  set<int> intToSet(int a) {
    return setList[a];
  }
  int setToInt(const set<int> &s) {
    map<set<int>, int>::const_iterator it = m.find(s);
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

#endif
