#include <shiftspace_auto/permutation_graph.hh>
#include <shiftspace_auto/permutative_autos_simple.hh>
#include <output/output.hh>

using namespace std;
typedef stream_output<permutation_graph> output_type;
// typedef vector<permutation_graph> output_type;

int main()
{
  ringmatrix<int> E_gm(2,2,vector<int>{1,1,1,0});
  // permutation_metadata P(E_gm,4);

  ringmatrix<int> E_bowtie(3,3,vector<int>{1,1,0, 1,0,1, 0,1,1});
  // permutation_metadata P(E_bowtie,4);

  ringmatrix<int> E_O2(1,1,vector<int>{2});
  //  permutation_metadata P(E_O2,4);

   ringmatrix<int> E_O4(1,1,vector<int>{4});
   permutation_metadata P(E_O4,2);       
  
  partial_permutation_graph G0(P);
  output_type output(cout);
  PermutativeInnerAutos<output_type> autosearch(P,output);
  
  cerr << P << endl;

  autosearch();
  
  return 0;
}
