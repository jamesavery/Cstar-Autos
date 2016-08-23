#include <auxiliary.hh>
#include <shiftspace_auto/invertibility.hh>

/**********************************************************************

  Implementation of threaded concatenation operator through triplets
  of wordgraphs.

 **********************************************************************/

triplet triplet::prune() const {
  WordGraph MAV(U.innodes(), (A * V).innodes()), UMV(U.outnodes(), V.innodes()),
      UAM((U * A).outnodes(), V.outnodes());

  return triplet(U & MAV, A & UMV, V & UAM);
}

bool triplet::has_parallel_words() const {
  return WordGraph::has_parallel_words(U, A) ||
         WordGraph::has_parallel_words(A, V) ||
         WordGraph::has_parallel_words(U, A * V);
}

bool triplet::operator<(const triplet &t) const {
  return U < t.U || (U == t.U && (A < t.A || (A == t.A && V < t.V)));
}

string triplet::to_latex(const WordGraph &L, const WordGraph &R) const {
  string s;
  int U_m = U.A.M, U_n = U.A.N, A_m = A.A.M, A_n = A.A.N, V_m = V.A.M,
      V_n = V.A.N;

  assert(U_n == A_m && A_n == V_m);

  vector<int> domains = {U_m, U_n, V_m, V_n};

  vector<string> triplet_edges, L_edges, R_edges;
  for (auto e : U.A.edge_list())
    triplet_edges.push_back("0/" + to_string(e.first) + "/1/" +
                            to_string(e.second));
  for (auto e : A.A.edge_list())
    triplet_edges.push_back("1/" + to_string(e.first) + "/2/" +
                            to_string(e.second));
  for (auto e : V.A.edge_list())
    triplet_edges.push_back("2/" + to_string(e.first) + "/3/" +
                            to_string(e.second));

  if (L.e_type != WordGraph::ONE) { // Visualize left-multiplication
    domains.push_back(L.A.M);
    for (int j = 0; j < L.A.M; j++)
      for (int k : L.A[j])
        L_edges.push_back(to_string(j) + "/" + to_string(k));

    s += "\\DrawtripletmulL{" + to_string(domains) + "}{" +
         to_string(triplet_edges) + "}{" + to_string(L_edges) + "}";
  } else if (R.e_type != WordGraph::ONE) { // Visualize right-multiplication
    domains.push_back(R.A.N);

    for (int j = 0; j < R.A.M; j++)
      for (int k : R.A[j])
        R_edges.push_back(to_string(j) + "/" + to_string(k));

    s += "\\DrawtripletmulR{" + to_string(domains) + "}{" +
         to_string(triplet_edges) + "}{" + to_string(R_edges) + "}";
  } else {
    s += "\\Drawtriplet{" + to_string(domains) + "}{" +
         to_string(triplet_edges) + "}";
  }

  return s;
}

/**********************************************************************

              Implementation of invertibility algorithm

 **********************************************************************/

Invertibility::Invertibility(int N,
                             const map<digraph::dedge_t, vector<int>> &labels_)
    : G(N, N, get_keys(labels_)) {
  L = 0;
  for (auto ea : labels_) {
    vector<string> ls_string;
    for (auto l : ea.second) {
      ls_string.push_back(string(1, 'a' + l));
      if (l + 1 > L)
        L = l + 1;
    }
    labels[ea.first] = to_string(ls_string);
  }

  letters = vector<WordGraph>(L, digraph(N, N));

  for (auto ea : labels_) {
    digraph::dedge_t e = ea.first;
    vector<int> ls = ea.second;

    for (auto a : ls)
      letters[a].set(e.second + 1,
                     e.first + 1); // TODO: set(j,i) -> set(i,j) overalt.
  }
}

string Invertibility::labeled_graph_latex(const string name) const {
  vector<string> edgelist;
  for (auto e : labels)
    edgelist.push_back(to_string(e.first.first) + "/" +
                       to_string(e.first.second) + "/" + e.second + "/" +
                       (e.first.first == e.first.second ? "Eloop" : "Eedge"));

  string s = "\\def\\" + name + "edgelist" + to_string(edgelist) + "\n" +
             "\\def\\" + name + "graph{\\drawgraph{" + to_string(G.N) + "}{\\" +
             name + "edgelist}}\n";

  return s;
}

string Invertibility::triplet_set_latex(const string name,
                                        const set<triplet> &S) {
  string r, s;
  int i = 0;
  if (!S.empty()) {
    for (auto t(S.begin()), t_last = --S.end(); t != S.end(); t++, i++) {
      r += "\\def\\" + name + "elem" + string(1, 'a' + i) + "{" +
           t->to_latex() + "}\n";
      s += "\\" + name + "elem" + string(1, 'a' + i) + (t != t_last ? "," : "");
    }
    return r + "\\def\\" + name + "set{\\left\\{" + s + "\\right\\}}\n";
  } else {
    return "\\def\\" + name + "set{\\emptyset}\n";
  }
}

string Invertibility::triplet_vector_latex(const string name,
                                           const vector<triplet> &S) {
  string r, s;
  int i = 0;
  if (!S.empty()) {
    for (auto t(S.begin()), t_last = --S.end(); t != S.end(); t++, i++) {
      r += "\\def\\" + name + "elem" + string(1, 'a' + i) + "{" +
           t->to_latex() + "}\n";
      s += "\\" + name + "elem" + string(1, 'a' + i) + (t != t_last ? "," : "");
    }
    return r + "\\def\\" + name + "set{\\left\\{" + s + "\\right\\}}\n";
  } else {
    return "\\def\\" + name + "set{\\emptyset}\n";
  }
}

string Invertibility::triplet_set_product_latex(const string name,
                                                const set<triplet> &S,
                                                bool left) const {
  string r, s;
  int i = 0;
  if (!S.empty()) {
    for (int a = 0; a < letters.size(); a++) {
      WordGraph W = letters[a];

      for (auto t(S.begin()), t_last = --S.end(); t != S.end(); t++, i++) {
        r += "\\def\\" + name + "elem" + string(1, 'a' + i) + "{" +
             t->to_latex(left ? W : WordGraph::one(),
                         left ? WordGraph::one() : W) +
             "}\n";
        s += "\\" + name + "elem" + string(1, 'a' + i) +
             (t != t_last ? "," : "");
      }
      s += (a + 1 < letters.size() ? "," : "");
    }
    return r + "\\def\\" + name + "set{\\left\\{" + s + "\\right\\}}\n";
  } else
    return "\\def\\" + name + "set{\\emptyset}\n";
}

string Invertibility::WG_latex(const WordGraph &W) {
  vector<string> edges;
  for (auto e : W.A.edge_list())
    edges.push_back(to_string(e.first) + "/" + to_string(e.second));
  return "\\DrawWG{{" + to_string(W.A.M) + "," + to_string(W.A.N) + "}}{" +
         to_string(edges) + "}";
}

string Invertibility::WGprod_latex(const WordGraph &L,
                                          const WordGraph &R) {
  vector<string> Ledges, Redges;
  for (auto e : L.A.edge_list())
    Ledges.push_back(to_string(e.first) + "/" + to_string(e.second));
  for (auto e : R.A.edge_list())
    Redges.push_back(to_string(e.first) + "/" + to_string(e.second));

  return "\\DrawWGprod{{" + to_string(L.A.M) + "," + to_string(L.A.N) + "," +
         to_string(R.A.N) + "}}" + "{" + to_string(Ledges) + "}{" +
         to_string(Redges) + "}\n";
}

bool Invertibility::invertible_demo() const {
  set<triplet> S, Snext;
  vector<triplet> Sfull;

  for (int a = 0; a < letters.size(); a++) {
    const triplet t =
        triplet(WordGraph::one(), letters[a], WordGraph::one()).prune();
    if (!t.is_anchored())
      S.insert(t);
    Sfull.push_back(t);
  }
  cout << labeled_graph_latex("E");
  for (int a = 0; a < letters.size(); a++)
    cout << "\\def\\WG" + string(1, 'a' + a) + "{" + WG_latex(letters[a]) +
                "}\n";
  cout << triplet_set_latex("SA", S);
  cout << triplet_vector_latex("SAfull", Sfull);

  for (int i = 1; !S.empty(); i++) {
    bool parallel_or_doublecycle =
        breadth_first_demo(S, Snext, i & 1, "S" + string(1, 'A' + i));
    if (parallel_or_doublecycle)
      return false;
    S = Snext;
  }
  return true;
}

bool Invertibility::breadth_first_demo(const set<triplet> &workset,
                                       set<triplet> &next, const bool even,
                                       const string name) const {
  next.clear();
  cout << triplet_set_product_latex(name + "prod", workset, !even);
  bool parallel_words = false, double_cycles = false;
  int p = 0;
  vector<triplet> full_result;
  for (int a = 0; a < letters.size(); a++, p++) {
    for (set<triplet>::const_iterator e = workset.begin(); e != workset.end();
         e++) {
      const triplet &x(*e);

      triplet t;

      if (even) {
        t = x * letters[a];
        if (WordGraph::has_parallel_words(x.V, letters[a])) {
          cerr << "Right parallel word - this should never happen.\n";
          abort();
        }
      } else {
        t = letters[a] * x;
        if (WordGraph::has_parallel_words(letters[a], x.U)) {
          cerr << "Left parallel word - this should never happen.\n";
          abort();
        }
      }

      if (t.has_parallel_words()) {
        parallel_words = true; // length m+n+1 PW
        cout << "\\def\\StoppingTriplet{" + t.to_latex() + "}\n"
             << "\\def\\StoppingReason{parallel word triplet}\n";
      }

      if ((t.U * t.A * t.V).countcycles() > 1) {
        double_cycles = true; // DC of length m+n+1 found
        WordGraph W = t.U * t.A * t.V;
        vector<string> edges;
        for (auto e : W.A.edge_list())
          edges.push_back(to_string(e.first) + "/" + to_string(e.second));

        cout << "\\def\\StoppingTriplet{" << t.to_latex() << "}\n"
             << "\\def\\StoppingWord{\\DrawWG{{" << W.A.M << "," << W.A.N
             << "}}{" + to_string(edges) + "}}\n"
             << "\\def\\StoppingReason{double cycle}\n";
      }

      full_result.push_back(t);
      if (!t.is_anchored())
        next.insert(t);
    }
  }
  cout << triplet_set_latex(name, next);
  cout << triplet_vector_latex(name + "full", full_result);

  return parallel_words || double_cycles;
}
