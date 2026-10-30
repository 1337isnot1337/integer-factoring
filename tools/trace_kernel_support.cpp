#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

uint32_t mul(uint32_t a, uint32_t b, uint32_t polynomial, uint32_t order) {
  uint32_t answer=0;
  while(b) {
    if(b&1) answer^=a;
    b>>=1; a<<=1;
    if(a&order) a^=polynomial;
  }
  return answer;
}
int main(int argc,char**argv) {
  if(argc<4) return 2;
  int t=std::stoi(argv[1]),max_k=std::stoi(argv[3]);
  if(t<2||t>24||max_k<1||max_k>=t) return 2;
  uint32_t polynomial=std::stoul(argv[2],nullptr,0),d=1u<<t,n=d-1;
  if(polynomial<d||polynomial>=2*d||!(polynomial&1)) return 2;
  std::vector<uint32_t> coefficients(n),next(n),conjugates(t);
  coefficients[0]=1; conjugates[0]=2;
  for(int i=1;i<t;++i) conjugates[i]=mul(conjugates[i-1],conjugates[i-1],polynomial,d);
  std::vector<uint32_t> factors(t,1);
  for(int j=0;j<max_k;++j) {
    std::fill(next.begin(),next.end(),0);
    for(uint32_t a=0;a<n;++a) if(coefficients[a]) {
      for(int i=0;i<t;++i) {
        uint32_t e=a+(1u<<i); if(e>=n)e-=n;
        next[e]^=mul(coefficients[a],factors[i],polynomial,d);
      }
    }
    coefficients.swap(next);
    uint64_t count=0; std::vector<uint64_t> weights(t+1);
    for(uint32_t e=0;e<n;++e)if(coefficients[e]) {++count; ++weights[__builtin_popcount(e)];}
    std::cout<<"{\"t\":"<<t<<",\"k\":"<<j+1<<",\"rank\":"<<count<<",\"weights\":[";
    for(int i=0;i<=t;++i)std::cout<<(i?",":"")<<weights[i];
    std::cout<<"]}"<<std::endl;
    for(int i=0;i<t;++i)factors[i]=mul(factors[i],conjugates[i],polynomial,d);
  }
  if(argc>4) {
    std::ofstream out(argv[4],std::ios::binary);
    for(uint32_t value:coefficients)for(int j=0;j<4;++j)out.put(char((value>>(8*j))&255));
    if(!out) return 3;
  }
}
