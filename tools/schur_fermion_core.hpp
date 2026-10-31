// Exact characteristic-two Schur-mode queries; see schur-fermion-queries.md.
// Appendix B source: anthropics/formal-math@16a5ce4. Apache-2.0.
#pragma once
#include <algorithm>
#include <array>
#include <cstdint>
#include <stdexcept>
#include <vector>

namespace schur_fermion {
struct Partition {
  std::vector<int> parts, arms, legs;
  explicit Partition(std::vector<int> p = {}) : parts(std::move(p)) {
    std::sort(parts.begin(), parts.end(), std::greater<int>());
    for (int i=0; i<(int)parts.size() && parts[i]>i; ++i) {
      arms.push_back(parts[i]-i-1);
      int count=0;
      for (int value: parts) count += value>i;
      legs.push_back(count-i-1);
    }
  }
  int rank() const { return arms.size(); }
};

inline bool invertible(std::array<uint64_t,64>& rows, int n) {
  std::array<uint64_t,64> basis{};
  for (int i=0;i<n;++i) {
    auto row=rows[i];
    while(row) {
      int lead=63-__builtin_clzll(row);
      if(!basis[lead]) { basis[lead]=row; break; }
      row ^= basis[lead];
    }
    if(!row) return false;
  }
  return true;
}

inline bool mode2(const Partition& parameter,const Partition& output,const Partition& input) {
  const int a=output.rank(), b=parameter.rank(), c=input.rank(), n=a+b+c;
  if(n>64) throw std::runtime_error("native packed determinant supports order at most 64");
  if(a>b+c || b>a+c || c>a+b) return false;
  std::array<uint64_t,64> rows{};
  for(int i=0;i<a;++i) {
    for(int j=0;j<b;++j)
      rows[i] |= uint64_t((parameter.legs[j]&~output.legs[i])==0) << (a+j);
    for(int j=0;j<c;++j)
      rows[i] |= uint64_t(output.legs[i]==input.legs[j]) << (a+b+j);
  }
  for(int i=0;i<b;++i) {
    for(int j=0;j<a;++j)
      rows[a+i] |= uint64_t((parameter.arms[i]&~output.arms[j])==0) << j;
    for(int j=0;j<c;++j)
      rows[a+i] |= uint64_t((parameter.arms[i]&input.legs[j])==0) << (a+b+j);
  }
  for(int i=0;i<c;++i) {
    for(int j=0;j<a;++j)
      rows[a+b+i] |= uint64_t(input.arms[i]==output.arms[j]) << j;
    for(int j=0;j<b;++j)
      rows[a+b+i] |= uint64_t((input.arms[i]&parameter.legs[j])==0) << (a+j);
  }
  return invertible(rows,n);
}

inline void partition_rec(int n,int least,std::vector<int>& p,std::vector<Partition>& out) {
  if(!n) { out.emplace_back(p); return; }
  for(int k=least;k<=n;++k) { p.push_back(k);partition_rec(n-k,k,p,out);p.pop_back(); }
}

inline std::vector<Partition> partitions(int n) {
  std::vector<Partition> result;
  std::vector<int> p;
  partition_rec(n,1,p,result);
  return result;
}
}
