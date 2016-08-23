#include <shiftspace_auto/permutation_graph.hh>

template<typename output_type> class PermutativeAutos {
public:
  
//
// r -- alpha -- t
// |             |
// e             f
// |             |   
// s --  beta -- q
//
  void Phi0(int m, int n, partial_permutation_graph &G) 
  {
    //    cerr << "Phi0("<<m<<", " << n << ", " << G.graph << "); remaining = " << G.remaining << ";\n";
    const node_t e    = m;
    const node_t  &s  = P.E_source[e], &r = P.E_range[e];
    const uint8_t &ns = P.Ns[s], &nr = P.Ns[r];
    arcstep_type reason;
    
    if(n==P.Etau_from[r].size()) return Phi1(m+1,G);
    else {
      node_t alpha = P.Etau_from[r][n];
      node_t t     = P.Etau_range[alpha];

      // printf("(e,alpha) = (%c,%s) : %d->%d->%d\n",e+'a',P.node_names[alpha].c_str(),s,r,t);
      // cout << "Etau_from["<<s<<"] = " << P.Etau_from[s] << "\n";
      for(int i=0;i<ns;i++){
	node_t beta_i = P.Etau_from[s][i];
	node_t q      = P.Etau_range[beta_i];

	//	printf("beta_%d (of %d) = %s/%d: %d->%d\n",i,ns,P.node_names[beta_i].c_str(),beta_i,s,q);
	vector<int> &remaining(G.remaining[beta_i*P.m+t]);
	for(int j=0;j<P.E(q,t);j++)
	  if(remaining[j]>0) {
	    node_t f_j    = P.Esum(q,t)+j;

	    permutation_graph_edge mu{alpha,beta_i, e, f_j};

	    G.push_edge(mu);
	    if(G.is_synchronizing(reason)){
	      remaining[j]--;
	      Phi0(m, n+1, G);
	      remaining[j]++;
	    }
	    G.pop_edge(mu);
	  } else {
	    reason = NO_REMAINING_ARCS;
	    // printf("beta,f-pair has already been used: (beta,f) = (%s,%c): %d->%d->%d\n",
	    // 	   P.node_names[beta_i].c_str(),P.Esum(q,t)+j+'a',s,q,t);
	  }
      }
    }
  }

  void Phi1(int m, partial_permutation_graph &G) 
  {
    //    cerr << "Phi1("<<m<<", "<<G.graph<<")\n";
    if(m==P.L) output.push_back(G.graph);
    else Phi0(m,0,G);
  }

  output_type operator()(){
    partial_permutation_graph G0(P);
    Phi1(0,G0);
    return output;
  }
  //--------
  output_type &output; // template type
  const permutation_metadata           &P;

  PermutativeAutos(const permutation_metadata& P, output_type &output) : P(P), output(output) {}
};


template<typename output_type> class PermutativeInnerAutos {
public:
  
//
// r_m -- o_mni -- t_n
//  |                |
// e_m               f
//  |                |   
// s_m --- beta ---  q
//
  void Psi_s(const int m, const int n, const int i, partial_permutation_graph &G) 
  {
    
    //    cerr << "Psi_s("<<m<<", " << n << ", " << i << ", " << G.graph << ")\n";
    const int     e_m  = m;
    const node_t  t_n  = n;
    const node_t &s_m  = P.E_source[e_m], &r_m = P.E_range[e_m];
    const node_t o_mni = P.Etau0sum(r_m,t_n)+i;

    // fprintf(stderr,"e|alpha = %c|%s:%d->%d->%d.\n",
    // 	    e_m+'a',P.node_names[o_mni].c_str(),s_m,r_m,t_n);
	    
    if(i==P.Etau0(r_m,t_n)) return Psi_r(m,n+1,G);
    
    arcstep_type reason;    
    for(node_t q=0,offset=0;q<P.m;offset+=P.Etau0(s_m,q),q++){
      // fprintf(stderr,"E(%d,%d) = %d, and Etau0(%d,%d) = %d.\n",
      // 	      q,t_n,P.E(q,t_n),
      // 	      s_m,q,P.Etau0(s_m,q)
      // 	      );

      if(P.E(q,t_n) > 0){
	const int lowest_unused_beta = G.lowest_unused  [s_m*P.m + q]; // beta:s_m->q
	const int lowest_unused_f    = G.lowest_unused_f[q*P.m + t_n]; // f:   q->t_n

       
	for(int k=0;k<P.Etau0(s_m,q) && k <= lowest_unused_beta;k++){
	  node_t beta = P.Etau0sum(s_m,q)+k;
	  assert(q == P.Etau_range[beta]);
	  
	  vector<int> &remaining(G.remaining[beta*P.m+t_n]);

	  //	  fprintf(stderr,"remaining(%d,%d) = ",beta,t_n); cerr << remaining << endl; \
	  
	  for(int j=0;j<P.E(q,t_n) && j<=lowest_unused_f;j++)
	    if(remaining[j]>0) {	      
	      node_t f    = P.Esum(q,t_n)+j;
	      permutation_graph_edge mu{o_mni,beta, e_m, f};

	      G.push_edge(mu);
	      if(G.is_synchronizing(reason)){
		if(k == lowest_unused_beta) G.lowest_unused  [s_m*P.m + q]++;	      
		if(j == lowest_unused_f)    G.lowest_unused_f[q*P.m + t_n]++;
		remaining[j]--;

		Psi_s(m, n, i+1, G);

		remaining[j]++;
		G.lowest_unused  [s_m*P.m + q] = lowest_unused_beta;
		G.lowest_unused_f[q*P.m + t_n] = lowest_unused_f;
	      }
	      G.pop_edge(mu);
	    } else {
	      reason = NO_REMAINING_ARCS;
	      // printf("beta,f-pair has already been used: (beta,f) = (%s,%c): %d->%d->%d\n",
	      // 	   P.node_names[beta_i].c_str(),P.Esum(q,t)+j+'a',s,q,t);
	    }
	}
      }
    }
  }
  
  void Psi_r(const int m, const int n, partial_permutation_graph &G)
  {
    //    cerr << "Psi_r("<<m<<", " << n << ", " << G.graph << ");\n";
    if(n == P.m) Psi1(m+1,G);
    else Psi_s(m,n,0,G);
  }

  void Psi1(const int m, partial_permutation_graph &G) 
  {
    //    cerr << "Psi_1("<<m<<", " << G.graph << ");\n";    
    if(m==P.L) output.push_back(G.graph);
    else Psi_r(m,0,G);
  }

  
  output_type operator()(){
    partial_permutation_graph G0(P);
    Psi1(0,G0);
    return output;
  }

  //--------
  output_type &output; // template type
  const permutation_metadata           &P;

  PermutativeInnerAutos(const permutation_metadata& P, output_type &output) : P(P), output(output) {}
};
