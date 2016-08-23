#include "WordGraph128.hh"

namespace WordGraph128 {

std::ostream &operator<<(std::ostream &s, const WordGraph &g) {
  uint8_t l = g.length();
  switch (g.element_type()) {
  case WordGraph::ZERO:
    s << "zero";
    break;
  case WordGraph::ONE:
    s << "one";
    break;
  case WordGraph::NORMAL:
    s << "{";
    for (int i = 1; i <= l; i++)
      s << int(g.get(i)) << (i + 1 <= l ? "," : "}");
    break;
  }
  return s;
}

int WordGraph::countarcs(int max) const {
  if (element_type() != NORMAL)
    return -1;

  int cnt = 0;
  for (int r = 1; r <= length() && cnt < max; r++)
    cnt += get(r) == 0 ? 0 : 1;
  return cnt;
}

int WordGraph::countcycles(int max) const {
  if (element_type() != NORMAL)
    return -1;

  int cnt = 0;
  for (int r = 1; r <= length() && cnt < max; r++) {
    uint8_t s = get(r);
    if (s == r)
      cnt++;
  }
  return cnt;
}

bool WordGraph::synchronizes() const {
  uint8_t first = 0;
  for (int r = 1; r <= length(); r++) {
    uint8_t s = get(r);
    if (s != 0) {
      if (first == 0)
        first = s;
      else if (s != first)
        return false;
    }
  }
  return true;
}

// WordGraph semiring multiplication:
//  gh = { (s,r) | (s,t)\in g and (t,r)\in h } for regular elements
//  g1 = 1g = g
//  g0 = 0g = 0
WordGraph WordGraph::operator*(const WordGraph &B) const {
  // Formal elements must be treated separately.
  if (element_type() == ONE)
    return WordGraph(B);
  if (B.element_type() == ONE)
    return WordGraph(*this);
  if (element_type() == ZERO || B.element_type() == ZERO)
    return zero();

  const int l = B.length();
  WordGraph product(0, l);
  for (int r = 1; r <= l; r++) {
    uint8_t Bs = B.get(r);
    uint8_t ABs = get(Bs);
    product.set(r, ABs);
  }

  return product;
}

uint128_t WordGraph::number() const { return data; }
}
