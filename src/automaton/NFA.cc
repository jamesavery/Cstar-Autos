#include "NFA.hh"
#include "DFA.hh"

#include <queue>

typedef digraph::dedge_t edge;

////////////////////////////////////////////////////////////////////
//	      EC(S) = S \cup {q | p--epsilon-->q \in S}           //
//	      closure = fix(EC, {state0})                         //
////////////////////////////////////////////////////////////////////
void NFAi::epsilon_closure(int state0, set<int>& closure) const {
  int size = closure.size();
  const set<int> &neighbours(step[state0][0]);
  if(neighbours.empty()) return;

  set_union(closure.begin(),closure.end(),neighbours.begin(),neighbours.end(),inserter(closure,closure.end()));

  if(size == closure.size()) return;
    
  for(set<int>::const_iterator q(neighbours.begin()); q!=neighbours.end(); q++)
    epsilon_closure(*q,closure);
}

set<int> NFAi::epsilon_closure(int state0) const { 
  set<int> result;
  result.insert(state0);
  epsilon_closure(state0,result);
  return result;
}


////////////////////////////////////////////////////////////////////
//	    Delta(P,a) = \cup_{p\in P} EC(step[p][a])
////////////////////////////////////////////////////////////////////
void NFAi::Delta(const set<int> &P, int transition, set<int> &Q) const {
  assert(transition>=1);
  for(set<int>::const_iterator p(P.begin()); p!=P.end();p++){
    set<int> qs = step[*p][transition];
    for(set<int>::const_iterator q(qs.begin()); q!=qs.end(); q++){
      Q.insert(*q);
      epsilon_closure(*q,Q);
    }
  }
}

set<int> NFAi::Delta(const set<int> &P, int a) const {
  assert(a>=1);
  set<int> Q;
  Delta(P,a,Q);
  return Q;
}


/////////////////////////////////////////////////////////////////
//	    Replace states with their epsilon closures         //
/////////////////////////////////////////////////////////////////
// NB: Der ser ud til at vaere noget galt med remove_epsilons(). Skal vi bruge den?
// NFAi NFAi::remove_epsilons() const {
//   set< set<int> > states;
//   for(int q=0;q<Nstates;q++)
//     states.insert(epsilon_closure(q));
    
//   NFAi A(states.size(),Nletters);
//   int Qi=0;			// Aendret fra 1. Gammel bug?
//   for(auto Q : states){
//     const vector< set<int> > &steps(A.step[Qi]);
//     for(int a=0;a<steps.size();a++) { // XXX:NB: Aendret fra a<Nletters. Gammel bug?
//       set<int> newsteps = steps[a];

//       for(auto q: Q)
// 	set_union(newsteps.begin(),newsteps.end(),
// 		  step[q][a].begin(),step[q][a].end(),
// 		  inserter(newsteps,newsteps.end()));
//     }
//     Qi++;
//   }
//   return A;
// }

digraph NFAi::skeleton() const
{
  //  set<edge>  edges;	// TODO: Can we use vector instead of set?
  vector<edge> edges;
  edges.reserve(Nstates*Nletters*2);
  for(int p=0;p<Nstates;p++)
    for(int a=0;a<=Nletters;a++){
      const set<int>& Q(step[p][a]);
      for(auto q: Q) 
	edges.push_back(make_pair(p,q));
	//	edges.insert(make_pair(p,*q));
    }
  sort(edges.begin(),edges.end());
  unique(edges.begin(),edges.end());
  return digraph(Nstates,Nstates,edges);    
}

map<edge,vector<int>> NFAi::label_map() const
{
  map<edge,vector<int> > labels;
    
  for(int p=0;p<Nstates;p++)
    for(int a=0;a<=Nletters;a++)
      for(int q: step[p][a])
	labels[{p,q}].push_back(a-1);

  return labels;
}

ringmatrix<int> NFAi::adjacency_matrix() const
{
  ringmatrix<int> A(Nstates,Nstates);
  map<edge,vector<int> > labels = label_map();

  for(auto &ea: labels){
    const edge &e = ea.first;
    int multiplicity = ea.second.size();
    
    A(e.first,e.second) = multiplicity;
  }

  return A;
}
/////////////////////////////////////////////////////////////////
// Removes all sinks, sources, and dead states from automaton. //
/////////////////////////////////////////////////////////////////
void NFAi::trim() {
  digraph G = skeleton();

  vector< vector<digraph::node_t> > SCCs(G.TarjanSCC());

  // if(debug){
  //   cerr << "graph = " << G << ";\n";
  //   cerr << "SCCs  = " << SCCs << ";\n";
  //   cerr << G.to_latex() << endl;
  // }

  vector<bool> removed(Nstates,false);
  for(int scc=0;scc<SCCs.size();scc++)
    if(SCCs[scc].size() == 1){	// Remove unless SCC is a loop.
      int u = SCCs[scc][0];
      const vector<int>& u_edges(G[u]);
      if(find(u_edges.begin(),u_edges.end(),u) == u_edges.end()) 
	removed[SCCs[scc][0]] = true;
    }
  
  // Reorder remaining states
  IDCounter<int> newid;
  for(int q=0;q<Nstates;q++) if(!removed[q]) newid.insert(q);
  
  //    printf("Trimming %d to %d states\n",Nstates,next_id);
  
  vector< vector< set<int> > > new_steps(newid.size(),vector<set<int> >(Nletters+1));
  for(int p=0;p<Nstates;p++)
    if(!removed[p]){
      int new_p = newid(p);
      
      for(int a=0;a<=Nletters;a++){
	const set<int>& Q(step[p][a]);
	for(auto q: Q)
	  if(!removed[q])
	    new_steps[new_p][a].insert(newid(q));
      }
    }
    
  step = new_steps;
  Nstates = newid.size();
}


//////////////////////////////////////////////////////////////////////
//		  NFA-representation of minimal DFA.                //
//   NOTE: Does *NOT* compute minimal NFA, which is PSPACE hard.    //
//////////////////////////////////////////////////////////////////////
NFAi NFAi::minimize() const 
{
  DFAi Dresult(DFAi(reverse()).reverse());
  NFAi result = Dresult.NFA();
  // if(NFAi::debug){
  //   vector<string> letters;
  //   for(int e=0;e<floor(sqrt(Nletters));e++)
  //     for(int f=0;f<floor(sqrt(Nletters));f++)
  // 	letters.push_back({char('a'+e),'/',char('a'+f)});

  //   Dresult.alphabet = letters;
  //   cerr << "Minimized NFA (DFA): "  << Dresult << endl
  // 	 << "Minimized NFA: " << result << "\n\n";
  // }
  return result;
}

/////////////////////////////////////////////////////////////////////
// Construct reverted automaton r(A) = { q--a-->p | p--a-->q \in A}//
/////////////////////////////////////////////////////////////////////
NFAi NFAi::reverse() const {
  //    cout << "Reverse\n";
  NFAi nfa(Nstates,Nletters);

  for(int p=0;p<Nstates;p++)
    for(int a=0;a<=Nletters;a++)
      for(set<int>::const_iterator q(step[p][a].begin()); q!=step[p][a].end(); q++)
	if(a==0)
	  nfa.insert_epsilon_transition(*q,p);
	else
	  nfa.insert_transition(*q,p,a-1);

  nfa.trim();
  return nfa;
}

//////////////////////////////////////////////////////////////////
// Construct automaton A* = { a--q-->b  | p--a-->q--b-->r \in A}//
//////////////////////////////////////////////////////////////////
NFAi NFAi::dual() const {
  //  NFAi A(remove_epsilons());
  assert(epsilon_free());
  NFAi A(*this);
  NFAi nfa(Nletters,Nstates);

  for(int p=0;p<A.Nstates;p++)
    for(int a=1;a<=Nletters;a++)
      for(auto q: A.step[p][a])
	for(int b=1;b<=Nletters;b++)
	  for(auto r: A.step[q][b])
	    nfa.insert_transition(a,b,q);
    
  return nfa;
}



/////////////////////////////////////////////////////////////////
//	     Construct combined automaton A \cupdot B          //
/////////////////////////////////////////////////////////////////
NFAi NFAi::disjoint_union(const NFAi& B) const {
  const NFAi& A(*this);

  assert(A.Nletters == B.Nletters); //assume same alphabet
  NFAi N(A.Nstates+B.Nstates, A.Nletters);
  for(size_t n = 0; n < A.step.size(); ++n) {
    for(size_t a = 0; a < A.step[n].size(); ++a) {
      for(set<int>::const_iterator it(A.step[n][a].begin()); it != A.step[n][a].end(); ++it) {
	if(a == 0) {
	  N.insert_epsilon_transition(n, *it);
	} else {
	  N.insert_transition(n, *it, a-1);
	}
      }
    }
  }
  for(size_t n = 0; n < B.step.size(); ++n) {
    for(size_t a = 0; a < B.step[n].size(); ++a) {
      for(set<int>::const_iterator it(B.step[n][a].begin()); it != B.step[n][a].end(); ++it) {
	if(a == 0) {
	  N.insert_epsilon_transition(A.Nstates+n, A.Nstates+(*it));
	} else {
	  N.insert_transition(A.Nstates+n, A.Nstates+(*it), a-1);
	}
      }
    }
  }
  return N;
}

// 
bool NFAi::is_deterministic() const
{
  if(!epsilon_free()) return false;

  for(int p=0;p<Nstates;p++)
    for(int a=1;a<=Nletters;a++)
      if(step[p][a].size() >= 2)
	return false;

  return true;
}

bool NFAi::epsilon_free() const
{
  for(int p=0;p<Nstates;p++)
    if(!step[p][0].empty())
      return false;

  return true;
}


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

///
bool NFAi::equivalent_to(const NFAi& B) const {
  const NFAi& A(*this);
  const NFAi& N(disjoint_union(B));
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
    set<int> P = pq.first, Q = pq.second;
    S.pop();

    // if(debug) 
    //   cerr << "pq = " << pq << endl;

    //this is equivalent to eps(p) != eps(q) since all states are final.
    if(P.empty() != Q.empty()){
      // if(debug){
      // 	cerr << "P = " << P << endl;
      // 	cerr << "Q = " << Q << endl;
      // }
      return false; //line 7+8
    }
    for(int a = 1; a <= N.Nletters; ++a) { //line 9 
      int p2Id = UF.find(N.Delta(P, a)), q2Id = UF.find(N.Delta(Q, a)); //line 10+11
      if(p2Id != q2Id) { //line 12 
	set<int> P2 = UF.intToSet(p2Id), Q2 = UF.intToSet(q2Id);
	UF.unite(P2, Q2); //line 13
	S.push(make_pair(P2, Q2)); //line 14
      } else 
	++j; 
    }
  }
  return true; //line 15
}

vector<string> NFAi::default_nodenames() const {
    vector<string> nodenames(Nstates); 
    for(int i=0;i<Nstates;i++) nodenames[i] = to_string(i);
    return nodenames;
  }

vector<string> NFAi::default_labels(bool pairs) const {
  vector<string> labels(Nletters+1);

  labels[0] = "\\epsilon";
  if(pairs){
    int sqrtL = round(sqrt(Nletters));
    //    for(int i=0;i<Nletters;i++) labels[i+1] = string({char('a'+(i%sqrtL)),'/',char('a'+(i/sqrtL))});
    for(int i=1;i<=Nletters;i++) labels[i] = string({char('a'+((i-1)%sqrtL)),char('a'+((i-1)/sqrtL))});
  } else
    for(int i=1;i<=Nletters;i++) labels[i] = to_string(i);

  return labels;
}

string NFAi::to_latex(vector<string> labels, vector<string> nodenames) const {
  if(labels.empty())    labels = default_labels();
  if(nodenames.empty()) nodenames = default_nodenames();

  ostringstream S;
  for(int p=0;p<Nstates;p++) S << "\\node[state] ("<<p<<"){};\n";

  vector<string> transitions,loops;

  for(int p=0;p<Nstates;p++)
    for(int a=0;a<=Nletters;a++){
      const set<int>& Q(step[p][a]);
      for(auto q: Q)
	if(p==q) loops.push_back(to_string(p)+"/"+labels[a]);
	else transitions.push_back(to_string(p)+"/"+labels[a]+"/"+to_string(q));
    }

  S << "\n\\foreach \\u/\\t/\\v in " << transitions << "{\n"
    << "\\edef\\body{\n"
    << "\t\\noexpand\\draw (\\u) edge[->] node[lbl]{\\t} (\\v);\n"
    << "}\\body\n"
    << "}\n";

  S << "\n\\foreach \\u/\\t in " << loops << "{\n"
    << "\\edef\\body{\n"
    << "\t\\noexpand\\draw (\\u) edge[->,loop] node[lbl]{\\t} (\\u);\n"
    << "}\\body\n"
    << "}\n";


  return S.str();
}

// u32string NFAi::hash_string() const {
//   vector<vector<int> > result(Nstates, vector<int>(Nletters+1,0));
//   u32string result_string(2+Nstates*(Nletters+1),0);

//   result_string[0] = Nstates;
//   result_string[1] = Nletters;
    
//   for(int p=0;p<Nstates;p++)
//     for(int a=0;a<=Nletters;a++)
//       result[p][a] = step[p][a].size();

//   sort(result.begin(),result.end());

//   for(int p=0,i=0;p<Nstates;p++)
//     for(int a=0;a<=Nletters;a++,i++)
//       result_string[2+i] = result[p][a];

//   return result_string;
// }

u32string NFAi::hash_string() const {
  vector<vector<int> > result(Nstates, vector<int>(2*(Nletters+1),0));
  u32string result_string(2+Nstates*2*(Nletters+1),0);

  NFAi rev(reverse());

  result_string[0] = Nstates;
  result_string[1] = Nletters;
    
  for(int p=0;p<Nstates;p++)
    for(int a=0;a<=Nletters;a++){
      result[p][2*a]   = step[p][a].size();
      result[p][2*a+1] = rev.step[p][a].size();
    }

  sort(result.begin(),result.end());

  for(int p=0,i=0;p<Nstates;p++)
    for(int a=0;a<2*(Nletters+1);a++,i++) 
      result_string[2+i] = result[p][a];

  return result_string;
}
