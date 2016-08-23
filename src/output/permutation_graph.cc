#include <output/permutation_graph.hh>

std::ostream& operator<<(std::ostream& S, const permutation_metadata &P)
{
  S << "(**** Permutative endomorphism metadata ****)\n"
    << "{k,m,M,L} = " << vector<int>{P.k,P.m,P.M,P.L} << ";\n"
    << "(**** Original graph E ****)\n"
    << "E = " << P.E << ";\n"
    << "Esum = " << P.Esum << ";\n"
    << "E_source = " << P.E_source << ";\n"
    << "E_range  = " << P.E_range << ";\n"
    << "E_from   = " << P.E_from << ";\n"
    << "E_to     = " << P.E_to   << ";\n"
    << "\n"
    << "(**** Permutation graph Etau ****)\n" 
    << "Etau0    = " << P.Etau0 << ";\n"
    << "Etau0sum = " << P.Etau0sum << ";\n"
    << "Ek       = " << P.Ek << ";\n"
    << "node_names = " << P.node_names <<";\n"							   
    << "Ns       = " << P.Ns <<";\n"
    << "Nr       = " << P.Nr <<";\n"
    << "Etau_source = " << P.Etau_source << ";\n"
    << "Etau_range  = " << P.Etau_range << ";\n"
    << "Etau_from   = " << P.Etau_from << ";\n"
    << "Etau_to     = " << P.Etau_to << ";\n"
    << "L1local     = " << P.L1local << ";\n"
    << "L2local     = " << P.L2local << ";\n";

  return S;
}
