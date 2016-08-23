#include <sstream>
#include <set>

#include "digraph.hh"

using namespace std;

digraph::digraph(int M, int N, const vector< vector<int> >& edges) : M(M), N(N), vector<vector<int> >(edges) 
{
  resize(M);
}

digraph::digraph(int M, int N, const vector<digraph::dedge_t>& edges) : M(M), N(N), vector<vector<int> >(M)
{
  digraph &A(*this);
  
  for(auto e: edges)
    A[e.first].push_back(e.second);

  for(int u=0;u<M;u++){
    vector<node_t> &row(A[u]);
    sort(row.begin(), row.end());
  }
}
						
digraph digraph::transpose() const {
  const digraph& A(*this);
  digraph At(N,M);

  for(int u=0;u<M;u++)
    for(int j=0;j<A[u].size();j++)
      At[A[u][j]].push_back(u);

  return At;
}

digraph digraph::operator*(const digraph& B) const
{
  const digraph &A(*this), Bt(B.transpose());
  digraph C(A.M,B.N);
    
  set<int> sum;
  for(int u=0;u<A.M;u++){
    const vector<node_t> &Arow(A[u]);
    sum.clear();
    for(int i=0;i<Arow.size();i++){
      const vector<node_t> &Bcol(B[Arow[i]]);
      for(int j=0;j<Bcol.size();j++)
	sum.insert(Bcol[j]);
    }
    C[u] = vector<node_t>(sum.begin(),sum.end());
  }
  return C;
}

vector< vector<int> > digraph::TarjanSCC() const {
  int index = 0;
  deque<int> S;
    
  vector<int> SCCindex(N,-1), lowlink(N,-1);
  vector< vector<node_t> > SCCs;
  for(int v=0;v<N;v++)
    if(SCCindex[v] == -1)
      TarjanSCC(v,S,SCCindex,lowlink,index,SCCs);
    
  return SCCs;
}

void digraph::TarjanSCC(node_t v, deque<node_t> &S, vector<int> &SCCindex,
                        vector<int> &lowlink, int &index,
                        vector<vector<node_t>> &SCCs) const {
  SCCindex[v] = index;
  lowlink[v] = index;
  index++;
  S.push_back(v);

  const vector<int> &kids((*this)[v]);

  for (int i = 0; i < kids.size(); i++) {
    node_t w = kids[i];

    if (SCCindex[w] == -1) {
      TarjanSCC(w, S, SCCindex, lowlink, index, SCCs);
      lowlink[v] = min(lowlink[v], lowlink[w]);
    } else if (find(S.begin(), S.end(), w) != S.end()) {
      lowlink[v] = min(lowlink[v], SCCindex[w]);
    }
  }

  if (lowlink[v] == SCCindex[v]) {
    vector<node_t> scc;
    node_t w = -1;
    while (w != v) {
      w = S.back();
      S.pop_back();
      scc.push_back(w);
    }
    SCCs.push_back(scc);
  }
}

bool digraph::fan_in_free() const {
  const digraph &A(*this);
  vector<bool> inarc(N);

  for(int u=0;u<M;u++) 
    for(int j=0;j<A[u].size();j++){
      const int &v(A[u][j]);
      if(!inarc[v]) inarc[v] = true;
      else return false;	// v receives at least two arcs.
    }
  return true;
}

bool digraph::fan_out_free() const {
  for(int i=0;i<M;i++) if((*this)[i].size() >= 2) return false;
  return true;
}
  

void digraph::insert(const digraph::dedge_t& e)
{
  node_t u = e.first, v = e.second;
  vector<node_t> &nu((*this)[u]);
  
  if(find(nu.begin(),nu.end(),v) == nu.end()) nu.push_back(v);
}


vector<digraph::dedge_t> digraph::edge_list() const 
{
  vector<dedge_t> edgelist;
  for(int u=0;u<M;u++) for(int v: (*this)[u]) edgelist.push_back(dedge_t(u,v));
  return edgelist;
}

  
string digraph::to_latex() const {
  std::ostringstream S("");
  const vector<vector<int> > &A(*this); 
  for(int p=0;p<A.size();p++) S << "\\node ("<<p<<"){"<<to_string(p)<<"};\n";

  vector<string> transitions,loops;
  for(int u=0;u<A.size();u++) 
    for(int v: A[u])
      if(u==v) loops.push_back(to_string(u));
      else     transitions.push_back(to_string(u)+"/"+to_string(v));

  S << "\n\\foreach \\u/\\v in " << transitions << "\n"
    << "\t\\draw (\\u) edge[->] (\\v);\n";

  S << "\n\\foreach \\u in " << loops << "\n"
    << "\t\\draw (\\u) edge[->,loop] (\\u);\n";

  return S.str();
}
