#ifndef AUXILIARY_HH
# define AUXILIARY_HH

#include <assert.h>
#include <limits.h>
#include <list>
#include <set>
#include <string>
#include <deque>
#include <iostream>
#include <sstream>
#include <map>
#include <vector>

using namespace std;

template <typename T> string to_string(const T& x) {
  ostringstream l;
  l << x;
  return l.str();
}

#define container_output(container) \
  template <typename T> ostream& operator<<(ostream& s, const container<T>& v) \
  { \
  s << "{"; \
  for(typename container<T>::const_iterator x(v.begin());x!=v.end();){ \
    s << *x; \
    if(++x!=v.end()) s << ","; \
  } \
  s << "}"; \
  return s; \
}

#define container_vectorize(container) \
  template <typename T> vector<T> vectorize(const container<T>& v){ return vector<T>(v.begin(),v.end()); }


container_output(vector);
container_output(list);
container_output(deque);
container_output(set);

container_vectorize(set);

#define container_flatten(c1, c2)\
  template <typename T> vector<T> flatten(const c1< c2<T> >& V){\
  int size = 0;\
  for(typename c1<c2<T> >::const_iterator v(V.begin()); v!=V.end();v++)		\
    for(typename c2<T>::const_iterator x(v->begin()); x!=v->end();x++) size++;	\
  vector<T> result(size); int i=0;					\
  for(typename c1<c2<T> >::const_iterator v(V.begin()); v!=V.end();v++)		\
    for(typename c2<T>::const_iterator x(v->begin()); x!=v->end();x++,i++)	\
      result[i] = *x;							\
  return result;							\
  }

container_flatten(vector,vector);

template <typename T> istream& operator>>(istream &S, vector<T>& xs) 
{ 
  xs.clear();
  char c; 
  T x; 
  S >> c; 
  assert(c=='{');		/* Set fail flag instead */ 
  while(c!='}' && !S.fail()){ 
    S >> x; 
    S >> c; 
    xs.push_back(x); 
  } 
  return S; 
} 

template <typename A, typename B> istream& operator>>(istream &S, pair<A,B>& xs) 
{ 
  char c; 
  S >> c; 
  assert(c=='{');
  S >> xs.first;
  S >> c;
  assert(c==',');
  S >> xs.second;
  S >> c;
  assert(c=='}');
  return S; 
} 



template <typename S, typename T> ostream& operator<<(ostream& s, const pair<S,T>& p)
{
  s << "{" << p.first << "," << p.second << "}";
  return s;
}

template<typename K, typename V> vector<K> get_keys(const map<K,V>& m)
{
  vector<K> keys(m.size());
  int i=0;
  for(typename map<K,V>::const_iterator kv(m.begin()); kv!=m.end(); kv++,i++)
    keys[i] = kv->first;
  return keys;
}

template<typename K, typename V> vector<K> get_keys(const vector<pair<K,V> >& m)
{
  vector<K> keys(m.size());
  int i=0;
  for(typename vector<pair<K,V> >::const_iterator kv(m.begin()); kv!=m.end(); kv++,i++)
    keys[i] = kv->first;
  return keys;
}


template<typename K, typename V> vector<V> get_values(const map<K,V>& m)
{
  vector<V> values(m.size());
  int i=0;
  for(typename map<K,V>::const_iterator kv(m.begin()); kv!=m.end(); kv++,i++)
    values[i] = kv->second;
  return values;
}

template <typename T> int sgn(const T& val) { return (T(0) < val) - (val < T(0)); }


template <typename C> int max(const C &v)
{
  int imax = INT_MIN;
  for(typename C::const_iterator i(v.begin()); i!=v.end();i++)
    if(*i > imax) imax = *i;
  return imax;
}

template <typename C> int min(const C &v)
{
  int imin = INT_MAX;
  for(typename C::const_iterator i(v.begin()); i!=v.end();i++)
    if(*i < imin) imin = *i;
  return imin;
}

template <typename T> class IDCounter: public map<T,int> {
public:
  int nextid;

  IDCounter(int start=0) : nextid(start) {}
  
  int insert(const T& x){
    typename map<T,int>::const_iterator it(map<T,int>::find(x));
    if(it != this->end()) return it->second;
    else {
      map<T,int>::insert(make_pair(x,nextid));
      return nextid++;
    }
  }

  int operator()(const T& x) const {
    typename map<T,int>::const_iterator it(map<T,int>::find(x));
    if(it != this->end()) return it->second;
    else return -1;
  }
};

#endif
