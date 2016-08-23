#include <shiftspace_auto/synchronizer.hh>

#ifdef VERBOSE_DEBUG
# define debugf(x) fprintf x
# define debug(x) x
#else 
# define debugf(x)
# define debug(x) 
#endif

/*
   r--alpha---t
   |          |
   e          f
   |          |
   s---beta---q
 */

Synchronizer::Synchronizer(const permutation_metadata &P, const permutation_dual &phi,
                           bool left)
    : P(P), left(left), LR(phi.resolving(left)), Psources(P.m), Pranges(P.m),
      WGEdges(P.L), AL(phi.automaton_component(true)),
      AR(phi.automaton_component(false)) {
  permutation_dual phiflat(phi.flatten());
  size_t maxclass = 0;
  M = phiflat.Nstates();
  local_name.resize(M);

  // Translate between local and global names;
  for (int e = 0; e < P.L; e++) {
    int s = P.E_source[e], r = P.E_range[e];
    for (int f = 0; f < P.L; f++) {
      int q = P.E_source[f], t = P.E_range[f];
      // beta|f <-> e|alpha
      // beta: s->q, alpha: r->t
      const PDualR &Aef(phi(e, f));

      for (const auto &g : Aef) {
        int beta_global(g.first[0]), alpha_global(g.second[0]);

        Psources[s].insert(beta_global);
        Pranges[q].insert(beta_global);
        Psources[r].insert(alpha_global);
        Pranges[t].insert(alpha_global);

        maxclass = max(maxclass, left ? Psources[s].size() : Pranges[q].size());

        if (left) {
          int beta  = Psources[s](beta_global)  + 1,
              alpha = Psources[r](alpha_global) + 1;
          WGEdges[e].insert(dedge_t(beta, alpha));
        } else {
          int beta  = Pranges[q](beta_global)  + 1,
              alpha = Pranges[t](alpha_global) + 1;
          WGEdges[f].insert(dedge_t(alpha, beta));
        }
      }
    }
  }
  fits64bit = (maxclass <= 15);
}

template <typename WGSet>
bool Synchronizer::synchronizes(const ringmatrix<WGSet> &A) const {
  typedef typename WGSet::wordgraph_t WordGraph;
  typedef ringmatrix<WGSet> WGMatrix;

  WGMatrix power(A), sum(P.m, P.m);

  for (int k = 1; k < 1000; k++) {
    unsigned int oldsize = sum.total_size();
    sum += power;

    if (sum.total_size() == oldsize)
      return true;

    power = A * power;
    for (int i = 0; i < P.m; i++)
      if (power.getelement(i, i).doublecycle().element_type() !=
          WordGraph::ZERO) {
        debug(cout << "Length-" << k + 1 << " word " << i << "->" << i
                   << " with double-cycle.\n");
        return false;
      }
  }
  return false;
}

bool Synchronizer::synchronizes() const {
  if (LR)
    if (fits64bit) {
      //	printf("WordGraph64\n");
      WordGraph64::WGMatrix A(letter_matrix<WordGraph64::WGSet>());
      return synchronizes(A);
    } else {
      //	printf("WordGraphLR\n");
      WordGraphLR::WGMatrix A(letter_matrix<WordGraphLR::WGSet>());
      return synchronizes(A);
    }
  else {
    //      printf("WordGraphFull\n");
    WordGraphFull::WGMatrix A(letter_matrix<WordGraphFull::WGSet>());
    return synchronizes(A);
  }
}

template <typename WG> vector<WG> Synchronizer::get_letters() const {
  vector<WG> letters(P.L);

  //    printf("get_letters(%s)\n",left?"left":"right");
  for (int e = 0; e < P.L; e++) {
    int s = P.E_source[e], r = P.E_range[e];
    int m = left ? Psources[s].size() : Pranges[r].size(),
        n = left ? Psources[r].size() : Pranges[s].size();

    //      printf("%c:%d->%d - %d x %d\n",e+'a',s,r,m,n);
    letters[e].setdimension(m, n);
    for (set<dedge_t>::const_iterator arc(WGEdges[e].begin());
         arc != WGEdges[e].end(); arc++)
      letters[e].set(arc->second, arc->first);
  }
  //      cout << "arcs: " << WGEdges << "\n"  << "letters: " << letters <<
  //      "\n\n";
  return letters;
}

template <typename WGSet>
ringmatrix<WGSet> Synchronizer::letter_matrix() const {
  typedef typename WGSet::wordgraph_t WordGraph;
  typedef ringmatrix<WGSet> WGMatrix;

  vector<WordGraph> letter_graphs(get_letters<WordGraph>());

  WGMatrix A(P.m, P.m);
  int e = P.L - letter_graphs.size();
   for (auto w: letter_graphs) {
    const node_t &s(P.E_source[e]), &r(P.E_range[e]);
    A.setelement(s, r, w);
    e++;
  }
  return left ? A : A.transpose();
}

/**********************************************************************

   Implementation of DualSynchronizer

**********************************************************************/
DualSynchronizer::DualSynchronizer(const permutation_metadata &P,
                                   const permutation_dual &phi0)
    : P(P), E(P.E, 2), left_resolving(phi.resolving(true)),
      right_resolving(phi.resolving(false)), phi(phi0.flatten()),
      dL(phi.Nstates()), sources(dL), ranges(dL), WGEdgesL(dL), WGEdgesR(dL),
      AL(P.L, dL), AR(P.L, dL), DL(P.m, P.m), DR(P.m, P.m) {
  // A^{L*}_{sp} = {beta(e,f):s->p}
  // A^{R*}_{rq} = {alpha(f,e):r->q}

  permutation_metadata E(P.E, 2);

  for (int e = 0; e < P.L; e++) {
    int s = P.E_source[e], r = P.E_range[e];
    for (int f = 0; f < P.L; f++) {
      int q = P.E_source[f], t = P.E_range[f];
      // beta|f <-> e|alpha
      // beta: s->q, alpha: r->t
      const PDualR &Aef(phi(e, f));

      for (const auto &g: Aef) {
        int beta(g.first[0]), alpha(g.second[0]);

        sources[beta] = s;
        ranges[beta]  = q;

        WGEdgesL[beta].insert(dedge_t(e, f));
        WGEdgesR[alpha].insert(dedge_t(f, e));

        AL.insert_transition(e, f, beta);
        AR.insert_transition(e, f, alpha);
      }
    }
  }

  for (int delta = 0; delta < dL; delta++) {
    int s = sources[delta], q = ranges[delta];

    WordGraphFull::WordGraph
      wl(digraph(E.Ns[s], E.Ns[q])),
      wr(digraph(E.Nr[s], E.Nr[q]));

    for (const auto &arc: WGEdgesL[delta]) {
      int el = E.L1local[arc.first], fl = E.L1local[arc.second];
      wl.set(fl, el);
    }
    for (const auto &arc: WGEdgesL[delta]) {
      int fr = E.L2local[arc.first], er = E.L2local[arc.second];
      wr.set(fr, er);
    }
    DL(s, q).insert(wl);
    DR(s, q).insert(wr);
  }
}

template <typename WGSet>
bool DualSynchronizer::synchronizes(const ringmatrix<WGSet> &A) const {
  typedef typename WGSet::wordgraph_t WordGraph;
  typedef ringmatrix<WGSet> WGMatrix;

  WGMatrix power(A), sum(P.m, P.m);

  for (int k = 1;; k++) {
    unsigned int oldsize = sum.total_size();
    sum += power;

    if (sum.total_size() == oldsize)
      return true;

    power = A * power;
    for (int i = 0; i < P.m; i++)
      if (power.getelement(i, i).doublecycle().element_type() !=
          WordGraph::ZERO) {
        return false;
      }
  }
  return false;
}
