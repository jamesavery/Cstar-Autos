#include "WordGraphFull.hh"

namespace WordGraphFull {

  std::ostream& operator<<(std::ostream& s, const WordGraph& g){
    switch(g.element_type()){
    case WordGraph::ZERO: s << "zero"; break;
    case WordGraph::ONE:  s << "one"; break;
    case WordGraph::NORMAL: 
      s <<  g.A;
      break;
    }
  return s;
}

  // TODO: max isf. max+1 i alle tre ordgraftyper
int WordGraph::countarcs(int max) const {
  if(element_type() != NORMAL) return -1;

  int cnt=0;
  for(int r=0;r<A.size() && cnt < max;r++) cnt += A[r].size();
  return min(cnt,max);
}


int WordGraph::countcycles(int max) const {
  if(element_type() != NORMAL) return -1;

  int cnt=0;
  for(int r=0;r<A.size();r++)
    for(int i=0;i<A[r].size();i++){
      if(A[r][i]==r) cnt++;
      if(cnt >= max) return cnt;
    }
  
  return cnt;
}

bool WordGraph::synchronizes() const {
  bool first = true;
  for(int r=0;r<A.size();r++)
    if(!A[r].empty()){
      if(first)
	first = false;
      else
	return false;
    }
  return true;
}

  WordGraph WordGraph::connected_arcs(const WordGraph& U, const WordGraph& a)
  {
    if(U.e_type == ZERO) return zero();
    if(U.e_type == ONE)  return a;
    if(a.e_type == ZERO) return zero();
    if(a.e_type == ONE)  return U;

    digraph Ut(U.A.transpose());
    digraph c(a.A);
    
    for(int w1=0;w1<Ut.size();w1++)
      if(Ut[w1].empty()) c[w1] = vector<int>();

    return c;
  }

  vector<bool> WordGraph::innodes() const { // Only works for normal-type word graphs
    vector<bool> us(A.size());

    for(int i=0;i<A.size();i++)
      us[i] = !A[i].empty();

    return us;
  }

  vector<bool> WordGraph::outnodes() const { return transpose().innodes(); }

  // Bit-matrix algoritme: 
  // Ut[w1] != 0 && V[w2] != 0 && a[w1][w2] != 0 for mere end eet par (w1,w2) <~>
  // 1. ut = innodes(Ut), v = innodes(V)
  // 2. UAV = ut*A*v
  // 3. test cnt(UaV) <= 1
  //
  // A*v  : kolonner i A, som er i innodes(V) 
  // ut*A : raekker i A som er i outnodes(U) ~> At*u
  //
  // connected_arcs(A,v) = A*v:
  //
  // v[i] = !(~V[i]) -> 0xffffffff hvis V[i] != 0, eller 0 hvis V[i] = 0
  // (A*v)[i] = A[i] & v[i]
  //
  // Slaa op i Knuth:
  // 1. Effektiv bit-matrix transposition
  // 2. Effektiv bit-count 
  //
  // is_anchored(U,A,V) = connected_arcs(connected_arcs(A,V),U^T) = connected_arcs(connected_arcs(U,A),V) 
  // 
  // connected_arcs(U,V) er et produkt med nulelement 0 og mange neutrale elementer, et af hvilke er ordgraf-1. Skriv U.V.
  //
  // Paastand: U*(U.V) = U*V = (U.V)*V, men for sparse repraesentationer er det hurtigere at regne (U.V)*V ud end U*V hvis man allerede har U.V, da der er faerre kanter.
  // Kommutativ? Nej.    A.B = {(a1,a2)\in A| \exists (a2,b2)\in B} != {(b1,b2)\in B| \exists (b2,a2)\in A} = B.A
  // Associativ? Ja. (A.B).C = {(a1,a2)\in A| \exist (a2,b2)\in B og \exists (b2,c2) \in C} = A.(B.C)
  // A.B.C.D = {(a1,a2) | (a1,a2,b2,c2,d2) \in A:B:C:D}
  //
  // og: is_anchored(U,A,V) = cnt({(a1,a2) | (u1,a1,a2,v2) \in U:A:V}) <= 1
  // og: {(a1,a2) | (u1,a1,a2,v2) \in U:A:V} = {(a1,a2) | (a1,a2,v2) \in A:V og (a2,a1,u1) in (U:A)^T} 
  //   = ((A.V)^T . U^T)^T. 
  // 
  // Hvad, hvis vi definerer det den anden vej? A%B = {(b1,b2) | (a1,b1,b2) \in A:B}
  // Da får vi: {(a1,a2) | (u1,a1,a2,v2) \in U:A:V} = {(a1,a2) \in U%A | (v2,a2,a1) \in (V:A)^T} = (V^T % (U%A)^T)^T
  // altså lige bøvlet. Dog: vi behøver ikke foretage den sidste transposition, da den ikke ændrer antallet af kanter.
  // Derfor:
  // is_anchored(U,A,V) = cnt((A.V)^T . U^T) <= 1
  //
  // Vent, meget pænere: 
  //
  // M = outnodes(U) x innodes(V)
  // is_anchored(U,A,V) = cnt(A & M) <= 1
  //
  // Mr = 1 x innodes(V)
  // Ml = outnodes(U) x 1
  // left_anchored(U,A)  = cnt(Ml & A) <= 1
  // right_anchored(A,V) = cnt(Mr & A) <= 1
  //
  // Algoritme:
  //
  // For hvert par (U,V) saa 


  bool WordGraph::has_parallel_words(const WordGraph& WA, const WordGraph& WB)
  {
    const digraph &A(WA.A);
    digraph Bt(WB.A.transpose());
    
    // Exists (u,v,w) and (u,v',w) in A:B
    for(int u=0;u<A.size();u++)
      for(int w=0;w<Bt.size();w++){
	bool connected = false;
	for(int i=0;i<A[u].size();i++)
	  for(int j=0;j<Bt[w].size();j++)
	    if(A[u][i] == Bt[w][j]){
	      if(connected) return true;
	      else connected = true;
	    }
      }
    return false;
  }

  // TODO: Rewrite to use connection matrices. Add handling of formal one() and zero().
  bool WordGraph::is_anchored(const WordGraph& U, const WordGraph& a, const WordGraph V)
  {
    bool connected = false;
    digraph Ut(U.A.transpose()); 

    for(int w1=0;w1<a.A.size();w1++){
      const vector<int>& nw1(a.A[w1]);
      for(int j=0;j<nw1.size();j++){
	int w2 = nw1[j];
	if(!Ut[w1].empty() && !V.A[w2].empty()) {
	  if(connected) return false; // Exists two (w1,w2)\in a such that (u,w1)\in U and (w2,v)\in V for some u,v
	  else connected = true;
	}
      }
    }
    return true;		// exists at most one (w1,w2)\in a such that (u,w1)\in U and (w2,v)\in V for some u,v.
  }

  bool WordGraph::left_anchored(const WordGraph& U, const WordGraph& a)
  {
    bool connected = false;
    digraph Ut(U.A.transpose());

    for(int w1=0;w1<Ut.size();w1++){
      if(!Ut[w1].empty() && !a.A[w1].empty()){
	if(connected) return false; // exists two w1 such that (u,w1)\in U, (w1,w2) \in a for some u, w2.
	else connected = true; 
      }
    }
    return true; // exists at most one w1 such that (u,w1)\in U and (w1,w2) \in a for some u, w2.
  }

  bool WordGraph::right_anchored(const WordGraph& a, const WordGraph& V)
  {
    bool connected = false;
    const digraph& at(a.A.transpose());

    for(int w2=0;w2<at.size();w2++){
      if(!at[w2].empty() && !V.A[w2].empty()){
	if(connected) return false; // exists two w2 such that (w1,w2)\in a and (w2,v)\in V for some w1,v.
	else connected = true;  
      }
    }
    return true; // exists at most one w1 such that (u,w1)\in U and (w1,w2) \in a for some u, w2.	
  }

  WordGraph WordGraph::transpose() const 
  {
    if(e_type == NORMAL) return A.transpose();
    else return *this;
  }

  WordGraph WordGraph::normalize(int n) const 
  {    
    if(e_type == ONE){
      WordGraph W(digraph(n,n),NORMAL);
      for(int i=0;i<n;i++) W.A.insert(digraph::dedge_t(i,i));
      return W;
    }
    if(e_type == ZERO) return WordGraph(digraph(n,n),NORMAL);
    else return WordGraph(*this);
  }

  bool WordGraph::left_anchored2(const WordGraph& U, const WordGraph& a)
  {
    return connected_arcs(U,a).countarcs(2) < 2;
  }

  bool WordGraph::right_anchored2(const WordGraph& a, const WordGraph& V)
  {
    return connected_arcs(V.transpose(),a.transpose()).countarcs(2);
  }

  WordGraph WordGraph::operator&(const WordGraph& mask) const 
  {
    if(e_type == ONE)      return mask;
    if(mask.e_type == ONE) return *this;
    if(e_type == ZERO || mask.e_type == ZERO) return zero();

    digraph C(A.M,A.N);

    assert(A.M == mask.A.M && A.N == mask.A.N);

    for(int u=0;u<A.M;u++){
      const vector<int> &Au(A[u]), &Mu(mask.A[u]);
      for(int i=0;i<Au.size();i++)
	for(int j=0;j<Mu.size();j++)
	  if(Au[i] == Mu[j]) C[u].push_back(Au[i]);
    }

    return C;
  }

// WordGraph semiring multiplication:
//  gh = { (s,r) | (s,t)\in g and (t,r)\in h } for regular elements
//  g1 = 1g = g
//  g0 = 0g = 0
WordGraph WordGraph::operator*(const WordGraph& B) const {
  // Formal elements must be treated separately.
  if(element_type() == ONE)   return WordGraph(B);
  if(B.element_type() == ONE) return WordGraph(*this);
  if(element_type() == ZERO || B.element_type() == ZERO)  return zero();

  return WordGraph(A*B.A);
}

}
