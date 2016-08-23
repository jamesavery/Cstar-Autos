#include "DFA.hh"

template <typename T> void merge_sets(list<T>& A, const set<T>& B)
{
  list<int>::iterator ia(A.begin());
  set<int>::const_iterator ib(B.begin());

  for(;ia!=A.end() && ib!=B.end();)
    if(*ia > *ib){		// a is the smallest element of A greater than b
      ia = A.insert(ia,*ib);
      ib++;
    } else if(*ib > *ia) {
      ia++;
    } else {			// *ia == *ib
       ib++;
       ia++;
     }

   if(ia == A.end())
     while(ib != B.end()){
       ia = A.insert(ia,*ib++);
       ia++;
     }
 }

void merge_sets(const vector<int>& A, const vector<int>& B, vector<int>& C)
{
  vector<int>::const_iterator ia(A.begin());
  vector<int>::const_iterator ib(B.begin());
  
  while(ia!=A.end() && ib !=B.end()){
    if(*ib < *ia)
      C.push_back(*ib++);
    else if(*ia < *ib)
      C.push_back(*ia++);
    else {
      C.push_back(*ia++);
      ib++;
    }
  }

  if(ia == A.end())
    while(ib != B.end())
      C.push_back(*ib++);
  
  if(ib == B.end())
    while(ia != A.end())
      C.push_back(*ia++);
}

 // Assumes nfa is epsilon-free
 DFAi DFAi::determinizeB(const NFAi& nfa) {
   int Nletters = nfa.Nletters, Nstates = 0;
   vector< vector<int> > step;		
   int id = 0;
   map<vector<int>, int> Statename;
   map<int,vector<int> > State;

   // digraph G = nfa.skeleton();
   // vector< vector<digraph::node_t> > SCCs(G.TarjanSCC());

   // TODO: Need to start either in every (initial) state, or in root 
   // of every connected component.    
   bool done = false;
   vector<int> state0;
   state0.push_back(0);
   Statename[state0] = 0;
   State[0] = state0;
   //   cout << "Determinizing NFA: " << nfa.Nstates << endl;
   deque<int> work;
   work.push_back(0);

   vector<int> Qnext,Q;
   Qnext.reserve(nfa.Nstates); 
   Q.reserve(nfa.Nstates);
   while(!work.empty()){
     //     cout << "Workset: " << work << "\n\n";
     int p = work.front(); work.pop_front();
     const vector<int> &P(State[p]);

     vector<int> row(Nletters,0);
     for(int a=0;a<Nletters;a++){
       Q.clear();

       for(vector<int>::const_iterator p(P.begin()); p!=P.end(); p++){
	 Qnext.clear();
	 const set<int> &Qp(nfa.step[*p][a+1]);
	 merge_sets(Q,vector<int>(Qp.begin(),Qp.end()),Qnext);
	 Q = Qnext;
       }
       
       if(!Q.empty()) {
	 if(Statename.find(Q) == Statename.end()){ // New state
	   Statename[Q] = ++id;
	   State[id] = Q;
	   work.push_back(id);
	 }
	 row[a] = Statename[Q]+1;
       }
     }

    step.push_back(row);
    Nstates = step.size();
  }

  DFAi dfa(Nstates,Nletters);
  dfa.step = step;
  dfa.trim();
  return dfa;  
}


DFAi DFAi::determinizeA(const NFAi& nfa) {
  cerr << "determinizeA(NFA) called!\n";
  int Nletters = nfa.Nletters, Nstates = 0;
  vector< vector<int> > step;		
  int id = 0;
  map<set<int>, int> Statename;
  map<int,set<int> > State;
  // TODO: Need to start either in every (initial) state, or in root 
  // of every connected component. 

  bool done = false;
  set<int> state0;
  state0.insert(0);
  Statename[state0] = 0;
  State[0] = state0;
  //    cout << "Determinizing NFA: " << nfa.Nstates;
  deque<int> work;
  work.push_back(0);
  while(!work.empty()){
    //      cout << "Workset: " << work << "\n\n";
    int p = work.front(); work.pop_front();
    const set<int> &P(State[p]);

    vector<int> row(Nletters,0);
    for(int a=0;a<Nletters;a++){
      set<int> Q(nfa.Delta(P,a));

      if(!Q.empty()) {
	if(Statename.find(Q) == Statename.end()){ // New state
	  Statename[Q] = ++id;
	  State[id] = Q;
	  work.push_back(id);
	}
	row[a] = Statename[Q]+1;
      }
    }

    step.push_back(row);
    Nstates = step.size();
  }
  
  DFAi dfa(Nstates,Nletters);
  dfa.step = step;
  dfa.trim();
  return dfa;
}


NFAi DFAi::NFA() const {
  NFAi nfa(Nstates,Nletters);

  for(int p=0;p<Nstates;p++)
    for(int a=0;a<Nletters;a++) 
      if(step[p][a] != 0)
	nfa.insert_transition(p,step[p][a]-1,a);

  return nfa;
}

NFAi DFAi::dual() const {
  NFAi nfa(Nletters,Nstates);

  for(int p=0;p<Nstates;p++)
    for(int a=0;a<Nletters;a++){
      int q = step[p][a]-1;
      if(q != -1)
	for(int b=0;b<Nletters;b++){
	  int r = step[q][a]-1;
	  if(r != -1)
	    nfa.insert_transition(a,b,q);
	}
    }
  return nfa;
}


void DFAi::trim() {			
  //    cout << "Trimming.\n";
  digraph G(skeleton());

  vector< vector<digraph::node_t> > SCCs(G.TarjanSCC());

   // TODO: Tjek at der ikke dukker DFA'er med mere end een ikke-triviel SCC op.

  vector<bool> removed(Nstates,false);
  for(int scc=0;scc<SCCs.size();scc++)
    if(SCCs[scc].size() == 1){	// Remove unless SCC is a loop. 
      int u = SCCs[scc][0];
      const vector<int>& u_edges(G[u]);
      if(find(u_edges.begin(),u_edges.end(),u) == u_edges.end()) 
	removed[SCCs[scc][0]] = true;
    }
    
  // Reorder remaining states
  map<int,int> newid;
  int next_id=0;
  for(int q=0;q<Nstates;q++) newid[q] = removed[q]? -1 : next_id++;

  vector< vector<int> > new_steps(next_id,vector<int>(Nletters));
  for(int p=0;p<Nstates;p++){
    if(!removed[p]){
      int new_p = newid[p];
	
      for(int a=0;a<Nletters;a++){
	int q = step[p][a];
	if(q != 0 && !removed[q-1])
	  new_steps[new_p][a] = newid[q-1]+1;
      }
    }
  }
  step = new_steps;
  Nstates = next_id;
}


NFAi DFAi::reverse() const {
  //    cout << "Reverse\n";
  NFAi nfa(Nstates,Nletters);

  for(int p=0;p<Nstates;p++)
    for(int a=0;a<Nletters;a++)
      if(step[p][a] != 0)
	nfa.insert_transition(step[p][a]-1,p,a);

  nfa.trim();
  return nfa;
}

digraph DFAi::skeleton() const
{
  set<digraph::dedge_t>  edges;// TODO: Can we use vector instead of set?
  for(int p=0;p<Nstates;p++)
    for(int a=0;a<Nletters;a++){
      int q = step[p][a];
      if(q!=0)
	edges.insert(make_pair(p,q-1));
    }

  return digraph(Nstates,Nstates,vectorize(edges));  
}


DFAi DFAi::reorder(const vector<int>& order) const
{
  DFAi A(*this,alphabet);

  for(int p=0;p<Nstates;p++)
    A.step[p] = step[order[p]];

  return A;
}

