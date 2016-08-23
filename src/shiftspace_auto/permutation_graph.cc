#include <shiftspace_auto/permutation_graph.hh>

//-- Reorganize
WGMatrix letter_matrix(const vector<WordGraph>& letters, const permutation_metadata &P)  {
  WGMatrix A(P.m,P.m);
  node_t e=P.L-letters.size();	// ???
  for(vector<WordGraph>::const_iterator l(letters.begin()); l!=letters.end();++l,e++){
    const node_t &s(P.E_source[e]), &r(P.E_range[e]); 
    A(s,r) += *l;
  }
  return A;
}


// 

bool synchronizes_quickly(const WGMatrix &A, int N)  
{
  WGMatrix power(A);
  for(int k=1;k<N;k++){		// Test if graph synchronizes for wordlengths <= 2^N
    if(power.total_size() == 0)
      return true;
    
    power = power*power;
  }
  return false;
}

bool no_short_doublecycles(const WordGraph& w) 
{
  WordGraph p(w);
  for(int k=0;k<10;k++){
    if(p.countcycles() >= 2){
      return false;
    }
    p = p*w;
  }
  return true;
}


bool no_shallow_doublecycles(const WGMatrix &A, int N)  
{
  WGMatrix power(A);
  for(int k=1;k<N;k++){

    if(power.total_size() == 0)
      return true;

    power = A*power;
    for(int i=0;i<A.m;i++)
      if(power.getelement(i,i).doublecycle().element_type() != WordGraph::ZERO){
	//	  cout << "Length-" << k+1 << " word "<<i<<"->"<<i<<" with double-cycle. Word-graph is " << *s.word << endl;
	return false;
      }
  }
  return true;
}

bool permutationgraph_is_synchronizing(const WGMatrix& A, size_t max_dc_length)  
{
  if(synchronizes_quickly(A,3)) return true;
  WGMatrix power(A), sum(A.m,A.m);

  for(size_t k=1;k<max_dc_length+1;k++){
    unsigned int oldsize = sum.total_size();
    sum += power;

    if(sum.total_size() == oldsize)
      return true;

    power = A*power;
    for(int i=0;i<A.m;i++)
      if(power.getelement(i,i).doublecycle().element_type() != WordGraph::ZERO){
	// Length-(k+1) word i->i with double-cycle.
	return false;
      }
  }
  return false;	
}


const char *arcstep_txt[] = { "UNDEFINED", "OK", "NO_REMAINING_ARCS", "SHORT_DC_LEFT", "SHORT_DC_RIGHT", "NONSYNC_LEFT", "NONSYNC_RIGHT"};
bool partial_permutation_graph::is_synchronizing(arcstep_type &reason) const
{  
  // Is G\oplus \mu a subgraph of an automorphism permutation graph?
  reason=OK;
  WGMatrix L(letter_matrix(graph_wg.first,P)), R(letter_matrix(graph_wg.second,P).transpose());
  
  if(!no_shallow_doublecycles(L,10))
    reason=SHORT_DC_LEFT;
  else if(!no_shallow_doublecycles(R,10))
    reason=SHORT_DC_RIGHT;
  else if(!permutationgraph_is_synchronizing(L,P.M*P.M)) // Maximal length of shortest DC. M(M-1)? M(M-1)/2? 
    reason=NONSYNC_LEFT;
  else if(!permutationgraph_is_synchronizing(R,P.M*P.M))
    reason=NONSYNC_RIGHT;

  // cout << "partial permutation graph:\n\t pg " << graph
  //      << "\n\t wg " << graph_wg << "\nis" << (reason==OK? " ":" not") << " synchronizing. "
  //      << "Reason: " << arcstep_txt[reason] << ".\n";
  return (reason==OK);
}

void partial_permutation_graph::push_edge(const permutation_graph_edge& mu)
{
  WordGraph
    &e_wg(graph_wg.first[mu.e]),
    &f_wg(graph_wg.second[mu.f]);
  
  const node_t
    &alpha1(P.L1local[mu.alpha]),
    &beta1 (P.L1local[mu.beta]),
    &alpha2(P.L2local[mu.alpha]),
    &beta2 (P.L2local[mu.beta]);

  // assert(alpha1 < e_wg.data.size());
  // assert(beta2  < f_wg.data.size());  

  // cout << "Pushing edge " << mu << "; (alpha1,beta1) = (" << alpha1 <<","<<beta1
  //      << "); (beta2,alpha2) = (" << beta2 << "," << alpha2 << ").\n";
  
  e_wg.set(alpha1,beta1);
  f_wg.set(beta2,alpha2);

  // cout << "e_wg = " << e_wg << "; f_wg = " << f_wg << endl;
  // cout << "Resulting wg: " << graph_wg << endl;
  
  graph.push_back(mu);
}

void partial_permutation_graph::pop_edge(const permutation_graph_edge& mu)
{
  assert(graph.back() == mu);

  graph.pop_back();

  WordGraph
    &e_wg(graph_wg.first[mu.e]),
    &f_wg(graph_wg.second[mu.f]);
  
  const node_t
    &alpha1(P.L1local[mu.alpha]),
    &beta2 (P.L2local[mu.beta]);

  // assert(alpha1 < e_wg.data.size());
  // assert(beta2  < f_wg.data.size());
  
  e_wg.set(alpha1,0);
  f_wg.set(beta2,0);
}

partial_permutation_graph::partial_permutation_graph(const permutation_metadata &P)
  : P(P), remaining(P.M*P.m), lowest_unused(P.m*P.m,0), lowest_unused_f(P.m*P.m,0)
{
  for(node_t beta=0;beta<P.M;beta++){
    node_t s = P.Etau_source[beta], q = P.Etau_range[beta];
    for(node_t t=0;t<P.m;t++)
      for(int i=0;i<P.E(q,t);i++)
	remaining[beta*P.m+t].push_back(1); // TODO: Tjek.
  }

  P.initialize_letters(graph_wg);
}



permutation_metadata::permutation_metadata(const ringmatrix<int> &E, int k)
  : E(E),  k(k), m(E.m), L(E.sum()), 
    Esum(m,m), E_source(L), E_range(L), E_from(m), E_to(m), indegree(m), outdegree(m), 
    Etau0(E.power(k - 1)), Ek(E.power(k)), Etau0sum(m,m), M(Etau0.sum()), Ns(m), Nr(m),
    Etau_source(M), Etau_range(M), Etau_from(m), Etau_to(m),
    L1local(M), L2local(M)
{
  for (node_t v = 0; v<m; v++) {
    indegree[v]  = E.colsum(v); // indegree[v]  = |E_{*->v}|
    outdegree[v] = E.rowsum(v); // outdegree[v] = |E_{v->*}|
    Nr[v] = Etau0.colsum(v);    // Nr[v] = |E^{k-1}_{*->v}|
    Ns[v] = Etau0.rowsum(v);    // Ns[v] = |E^{k-1}_{v->*}|      
  }

  for (node_t s = 0, alpha = 0, e=0; s < m; s++)
    for (node_t r = 0; r<m; r++) {
      int Esr = E(s, r);
      Esum(s,r)     = Esr>0? e : -1;


      for(int i=0;i<Esr; i++, e++){
	E_source[e] = s;
	E_range[e]  = r;
	E_from[s].push_back(e);
	E_to[r].push_back(e);
      }

      int Psr = Etau0(s, r);
      Etau0sum(s,r) = Psr>0? alpha : -1;

      for (int i = 0; i < Psr; i++, alpha++) {
	Etau_source[alpha] = s;
	Etau_range[alpha]  = r;
	Etau_from[s].push_back(alpha);
	Etau_to[r].push_back(alpha);
        L1local[alpha] = Etau_from[s].size();
        L2local[alpha] = Etau_to[r].size();	
      }
    }
  node_names = path_names(k-1);
}


bool operator <(const permutation_graph_edge& e1,const permutation_graph_edge& e2)  { return vector<int>({e1.alpha,e1.beta,e1.e,e1.f}) < vector<int>({e2.alpha,e2.beta,e2.e,e2.f}); }
bool operator==(const permutation_graph_edge& e1, const permutation_graph_edge& e2) { return vector<int>({e1.alpha,e1.beta,e1.e,e1.f}) == vector<int>({e2.alpha,e2.beta,e2.e,e2.f}); }
ostream& operator<<(ostream& stream,const permutation_graph_edge& e){ stream << vector<int>{e.alpha,e.beta,e.e,e.f}; return stream; }


ringmatrix<wordset> permutation_metadata::letter_strings() const {
  ringmatrix<wordset> A(m, m);
  for (node_t s = 0; s < m; s++)
    for (node_t t = 0; t < m; t++)
      for (int i = 0; i < E(s,t); i++)
        A(s, t) += wordset(string(1, Esum(s,t)+i + 'a'));

  return A;
}

vector<string> permutation_metadata::path_names(int k) const {
  ringmatrix<wordset> A(letter_strings());
  ringmatrix<wordset> Ak = A.power(k);

  vector<string> names;
  for (node_t u = 0; u < A.m; u++)
    for (node_t v = 0; v < A.n; v++)
      names.insert(names.end(), Ak(u, v).begin(), Ak(u, v).end());
  return names;
}


permutation_graph permutation_metadata::identity() const
{
  permutation_graph Id;

  ringmatrix<wordset> E(letter_strings());
  ringmatrix<wordset> nodes(E.power(k - 1));
  ringmatrix<wordset> paths(E.power(k));

  map<word, int> global_id;

  for (node_t s = 0, alpha_global = 0; s < m; s++)
    for (node_t r = 0; r < m; r++)
      for (auto w(nodes(s, r).begin()); w != nodes(s, r).end(); w++, alpha_global++)
        global_id[*w] = alpha_global;

  // For each path in paths, find prefix beta and suffix alpha elements in Ps
  // then look up end and start letters f and e end add edges alpha--beta.
  // Don't care about speed, so use simple but slow method.
  for (node_t S = 0; S < m; S++)
    for (node_t R = 0; R < m; R++) {
      const wordset &ps(paths(S, R));
      for (const string &path : ps) {
        string alpha_string = path.substr(1, path.size());
        string beta_string = path.substr(0, path.size() - 1);
        node_t alpha_global = global_id[alpha_string],
	  beta_global = global_id[beta_string];
        char e = path[0] - 'a', f = path[path.size() - 1] - 'a';

        Id.push_back(permutation_graph_edge({alpha_global, beta_global, e, f}));
      }
    }

  sort(Id.begin(), Id.end());
  return Id;  
}


void permutation_metadata::initialize_letters(permutation_graph_wg& letters) const
{
  letters.first.resize(L);
  letters.second.resize(L);

  // e:s->r, f:q->t
  // beta: s->q, alpha: r->t
  for(int e=0;e<L;e++)
    letters.first [e] = WordGraph(0,Ns[E_range[e]]); // Numbers of alpha starting in r
  for(int f=0;f<L;f++)
    letters.second[f] = WordGraph(0,Nr[E_source[f]]); // Number of betas ending in q
}
