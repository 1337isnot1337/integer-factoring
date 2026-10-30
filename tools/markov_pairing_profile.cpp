#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

#ifndef MAX_PROFILE_BITS
#define MAX_PROFILE_BITS 12
#endif
constexpr int MAXWORDS=1<<((MAX_PROFILE_BITS>6?MAX_PROFILE_BITS:6)-6);
using Kernel=std::array<uint64_t,MAXWORDS>;
struct KernelHash {
  size_t operator()(const Kernel& key) const {
    uint64_t h=0x9e3779b97f4a7c15ULL;
    for(auto word:key) h^=word+0x9e3779b97f4a7c15ULL+(h<<6)+(h>>2);
    return h;
  }
};
struct Ranker {
  int t,n,words;
  std::array<uint64_t,6> lower{};
  std::unordered_map<Kernel,uint16_t,KernelHash> cache;
  uint64_t calls=0,gaussian=0,hits=0;
  Ranker(int bits):t(bits),n(1<<bits),words((n+63)/64) {
    for(int j=0;j<6;j++) for(int h=0;h<64;h++)
      if(!(h&(1<<j))) lower[j]|=uint64_t(1)<<h;
    cache.reserve(1<<15);
  }
  uint64_t permute(uint64_t word,int shift) const {
    for(int j=0;j<6;j++) if(shift&(1<<j)) {
      int d=1<<j;
      word=((word&lower[j])<<d)|((word&~lower[j])>>d);
    }
    return word;
  }
  int rank(const Kernel& kernel) {
    calls++;
    uint64_t nonzero=0;int parity=0;
    for(int w=0;w<words;w++) {nonzero|=kernel[w];parity^=__builtin_parityll(kernel[w]);}
    if(!nonzero) return 0;
    if(parity) return n;
    for(int j=0;j<t;j++) {
      int linear=0;
      for(int w=0;w<words;w++)
        if(j<6) linear^=__builtin_parityll(kernel[w]&~lower[j]);
        else if(w&(1<<(j-6))) linear^=__builtin_parityll(kernel[w]);
      if(linear) return n/2;
    }
    auto it=cache.find(kernel);
    if(it!=cache.end()) {hits++;return it->second;}
    std::array<Kernel,64> permutations{};
    int low_count=std::min(n,64);
    for(int low=0;low<low_count;low++)
      for(int w=0;w<words;w++) permutations[low][w]=permute(kernel[w],low);
    std::vector<Kernel> basis(n);
    std::vector<uint8_t> present(n);
    int answer=0;
    for(int x=0;x<n;x++) {
      Kernel row{};
      int low=x&63,high=x>>6;
      for(int w=0;w<words;w++) row[w]=permutations[low][w^high];
      int top=words-1;
      while(top>=0) {
        if(!row[top]) {top--;continue;}
        int pivot=64*top+63-__builtin_clzll(row[top]);
        if(!present[pivot]) {present[pivot]=1;basis[pivot]=row;answer++;break;}
        for(int w=0;w<=top;w++) row[w]^=basis[pivot][w];
      }
    }
    gaussian++;
    if(cache.size()>(1u<<18)) cache.clear();
    cache.emplace(kernel,uint16_t(answer));
    return answer;
  }
};

int main(int argc,char** argv) {
  assert(argc==16);
  bool polynomial_edges=std::string(argv[1])=="polyprofile";
  int t=std::stoi(argv[2]);assert(t>=1&&t<=MAX_PROFILE_BITS);
  std::array<int,3> vertex{};
  std::array<int,9> edge{};
  for(int j=0;j<3;j++) vertex[j]=std::stoi(argv[4+j]);
  for(int j=0;j<9;j++) edge[j]=std::stoi(argv[7+j]);
  int maximum_edge=*std::max_element(edge.begin(),edge.end());
  if(polynomial_edges) maximum_edge=31-__builtin_clz(maximum_edge);
  int degree=t*(*std::max_element(vertex.begin(),vertex.end()))+(t-1)*maximum_edge;
  assert(degree<128);
  int n=1<<t;
  Ranker ranker(t);
  std::vector<std::vector<Kernel>> rows(degree+1,std::vector<Kernel>(n));
  std::vector<std::vector<Kernel>> cols(degree+1,std::vector<Kernel>(n));
  for(int g=0;g<n;g++) for(int h=0;h<n;h++) {
    using Polynomial=unsigned __int128;
    Polynomial dp[3]={(g&1)?0:Polynomial(1)<<vertex[0],
                      Polynomial(1)<<vertex[1],
                      (h&1)?0:Polynomial(1)<<vertex[2]};
    for(int i=1;i<t;i++) {
      Polynomial next[3]={};
      for(int s=0;s<3;s++) {
        if((s==0&&(g&(1<<i)))||(s==2&&(h&(1<<i)))) continue;
        for(int r=0;r<3;r++) {
          if(polynomial_edges) {
            unsigned factor=edge[3*r+s];
            while(factor) {
              int j=__builtin_ctz(factor);factor&=factor-1;
              next[s]^=dp[r]<<(vertex[s]+j);
            }
          } else next[s]^=dp[r]<<(vertex[s]+edge[3*r+s]);
        }
      }
      std::copy(next,next+3,dp);
    }
    Polynomial polynomial=dp[0]^dp[1]^dp[2];
    assert((__builtin_parityll(uint64_t(polynomial))^__builtin_parityll(uint64_t(polynomial>>64)))==(g==h));
    while(polynomial) {
      int j=uint64_t(polynomial)?__builtin_ctzll(uint64_t(polynomial)):64+__builtin_ctzll(uint64_t(polynomial>>64));
      polynomial&=polynomial-1;
      rows[j][g][h>>6]|=uint64_t(1)<<(h&63);
      cols[j][h][g>>6]|=uint64_t(1)<<(g&63);
    }
  }
  std::vector<int> rr(degree+1),cc(degree+1);
  uint64_t maximum_product=0;long double rank_sum=0;
  for(int j=0;j<=degree;j++) {
    for(int g=0;g<n;g++) {
      rr[j]+=ranker.rank(rows[j][g]);
      cc[j]+=ranker.rank(cols[j][g]);
    }
    uint64_t p=uint64_t(rr[j])*cc[j];
    maximum_product=std::max(maximum_product,p);
    rank_sum+=std::sqrt((long double)p);
    std::cerr<<"degree "<<j<<"/"<<degree<<" ranks "<<rr[j]<<","<<cc[j]
             <<" gaussian "<<ranker.gaussian<<std::endl;
  }
  auto output=[&](std::ostream& out) {
    out<<"{\"t\":"<<t<<",\"polynomial_edges\":"<<(polynomial_edges?"true":"false")<<",\"vertex\":[";
    for(int j=0;j<3;j++) out<<(j?",":"")<<vertex[j];
    out<<"],\"edge\":[";
    for(int j=0;j<9;j++) out<<(j?",":"")<<edge[j];
    out<<"],\"row_ranks\":[";
    for(int j=0;j<=degree;j++) out<<(j?",":"")<<rr[j];
    out<<"],\"column_ranks\":[";
    for(int j=0;j<=degree;j++) out<<(j?",":"")<<cc[j];
    out<<"],\"maximum_rank_product\":"<<maximum_product
       <<",\"rank_sum_estimate\":"<<std::setprecision(18)<<rank_sum
       <<",\"maximum_product_rho\":"<<std::log2((long double)maximum_product)/(4*t)
       <<",\"rank_sum_rho\":"<<std::log2(rank_sum)/(2*t)
       <<",\"rank_calls\":"<<ranker.calls<<",\"gaussian_ranks\":"<<ranker.gaussian
       <<",\"cache_hits\":"<<ranker.hits<<"}\n";
  };
  std::ofstream file(argv[3]);assert(file);output(file);output(std::cout);
}
