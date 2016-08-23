#pragma once

#include <output/permutation_graph.hh>

template <typename T> class stream_output {
public:
  ostream &stream;
  
  stream_output(ostream &stream) : stream(stream) {}

  void push_back(const T& x){ stream << x << ";\n"; }

  // TODO: Read?
};

// TODO: Binary checksummed i/o type.
// TODO: DFA i/o type.
