#pragma once

#include <type_traits>
#include <inttypes.h>
#include <map>
#include <iostream>
#include <assert.h>
#include <vector>
#include <stdio.h>

// A matrix over some ring R.
// R type must implement constructor from integers, which initialized with 0 resp 1
// yield the zero- resp. the one-element. No other integers need to be defined. 
template <typename R> class ringmatrix {
 public:
  typedef R value_t;
  typedef std::pair<size_t,size_t> index_t;
  size_t m, n;
  std::vector<R> a;

  static const ringmatrix zero(size_t n){ return ringmatrix(n,n); }
  static const ringmatrix one(size_t n){
    ringmatrix o(n,n);
    const R one_element = R(1);
    for(size_t i=0;i<n;i++) o.setelement(i,i,one_element);
    return o;
  }
  
  ringmatrix(size_t m=1, size_t n=1) : m(m), n(n), a(m*n) {}
  ringmatrix(size_t m, size_t n, const std::vector<R>& a) : m(m), n(n), a(a.begin(),a.end()) {
    assert(a.size() == m*n);
  }

  ringmatrix(const ringmatrix& A) : m(A.m), n(A.n), a(m*n) {
    for(int i=0;i<m*n;i++) a[i] = A.a[i];
  }
  
  void setelement(size_t i, size_t j, const R& aij){
    a[i*n+j] = aij;
  }
  
  const R getelement(size_t i, size_t j) const {
    return a[i*n+j];
  }
  
  const std::vector<R> getrow(size_t i) const {
    return std::vector<R>(&a[i*n],&a[(i+1)*n]);
  }

  const std::vector<R> getcol(size_t j) const {
    std::vector<R> col(m);
    for(int i=0;i<m;i++) col[i] = getelement(i,j);
    return col;
  }

  R& operator()(size_t i, size_t j)             { return a[i*n+j]; }
  const R& operator()(size_t i, size_t j) const { return a[i*n+j]; }

  ringmatrix operator *(const ringmatrix& y) const {
    assert(n == y.m);
    ringmatrix z(m,y.n);

    for(int i=0;i<m;i++)
      for(size_t j=0;j<n;j++){
	R sum(0);
	for(size_t k=0;k<y.n;k++)
	  sum += a[i*n+k]*y.a[k*m+j];
	
	z.setelement(i,j,sum);
      }
    return z;
  }
  
  ringmatrix transpose() const {
    ringmatrix T(n,m);
    for(size_t i=0;i<n;i++)
      for(size_t j=0;j<m;j++)
	T(i,j) = getelement(j,i);
    return T;
  }

  ringmatrix& operator +=(const ringmatrix& y)  {
    for(size_t i=0;i<m;i++)
      for(size_t j=0;j<n;j++)
	a[i*n+j] += y.getelement(i,j);

    return *this;
  }

  ringmatrix& operator -=(const ringmatrix& y)  {
    for(size_t i=0;i<m;i++)
      for(size_t j=0;j<n;j++)
	a[i*n+j] -= y.getelement(i,j);

    return *this;
  }

  ringmatrix operator +(const ringmatrix& y) const  {
    ringmatrix r(*this);
    r += y;
    return r;
  }

  ringmatrix operator -(const ringmatrix& y) const  {
    ringmatrix r(*this);
    r -= y;
    return r;
  }


  ringmatrix& operator *=(const ringmatrix& y) {
    return (*this = (*this)*y);
  }

  R rowsum(unsigned int i) const {
    R r(0);
    for(int j=0;j<n;j++) r += getelement(i,j);
    return r;
  }

  R colsum(unsigned int j) const {
    R r(0);
    for(int i=0;i<n;i++) r += getelement(i,j);
    return r;
  }

  R sum() const {
    R r(0);
    for(int i=0;i<m;i++) r += rowsum(i);
    return r;
  }
  
  size_t total_size() const {
    size_t sum = 0;
    for(int i=0;i<m*n;i++) sum += a[i].size();
    return sum;
  }

  void power(ringmatrix& A, const unsigned int n) const {
    if(n==0){ A = one(m); return; }
    if(n==1){ A = *this; return; }
    if(n&1){
      power(A,n-1);
      A = (*this)*A;
    } else {
      power(A,n>>1);
      A = A*A;
    }
  }

  ringmatrix power(const unsigned int N) const 
  {
    if(N==0) return one(m); 
    if(N==1) return *this;

    ringmatrix A(one(m));
    power(A,N);
    return A;
  }

  bool operator ==(const ringmatrix& y) const  {
    bool equal = true;
    for(size_t i=0;i<m;i++)
      for(size_t j=0;j<n;j++)
 	equal &= (getelement(i,j) == y.getelement(i,j));

    return equal;
  }

  bool operator <(const ringmatrix& y) const  {
    for(size_t i=0;i<m;i++)
      for(size_t j=0;j<n;j++){
	const R &X(getelement(i,j)), &Y(y.getelement(i,j));
	if(X < Y) return true;
	if(Y < X) return false;
      }
    
    return false;
  }

  bool operator !=(const ringmatrix& y) const  { return !(*this == y); }

  friend std::ostream& operator<<(std::ostream&s, const ringmatrix<R>& A){
    using namespace std;

    s << "{";
    for(size_t i=0;i<A.m;i++){
      s << A.getrow(i) << (i+1<A.m? ",":"");
    }
    s << "}";
    return s;
  }
};

