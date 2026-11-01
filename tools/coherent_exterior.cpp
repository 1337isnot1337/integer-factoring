// Exact characteristic-two coherent boundary profiles for exterior monodromy.
// Uses the exclusion walk and packed elimination of the repository's
// tools/matrix_unit_monodromy.cpp, extended to every weight-zero boundary.
// Source revision: bb3bf763ed054f6fdaf882649c037ad810e16b97 (PR 26).
#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

static const int states[6]={3,5,6,9,10,12};
static int index_of(int state) {
  for(int i=0;i<6;i++) if(states[i]==state) return i;
  return -1;
}

// Byte i is the polynomial coefficient for boundary basis element i.
// Six diagonal functionals followed by the two directed matching coherences.
static uint64_t entry(int t,uint32_t p,uint32_t q) {
  uint64_t result=0;
  for(int init=0;init<6;init++) {
    unsigned state=states[init];
    int a=0,b=0;
    bool alive=true;
    uint32_t pp=p,qq=q;
    for(int i=0;i<t;i++,pp>>=2,qq>>=2) {
      int x=pp&3,y=qq&3;
      if(x==y) a+=(state>>x)&1;
      else {
        if(!((state>>y)&1)||((state>>x)&1)) {alive=false;break;}
        state^=(1u<<x)|(1u<<y);b++;
      }
    }
    if(!alive) continue;
    int final=index_of(state),coordinate=-1;
    if(final==init) coordinate=init;
    else if(init==2&&final==3) coordinate=6;
    else if(init==3&&final==2) coordinate=7;
    if(coordinate<0) continue;
    unsigned polynomial=1u<<a;
    for(int i=0;i<b;i++) polynomial^=polynomial<<1;
    assert(polynomial<256);
    result^=uint64_t(polynomial)<<(8*coordinate);
  }
  return result;
}

static uint64_t binary_rank(std::vector<std::vector<uint64_t>>& rows,int n) {
  int words=(n+63)/64;
  std::vector<int> pivot(n,-1);
  uint64_t answer=0;
  for(int i=0;i<n;i++) {
    auto& row=rows[i];
    int top=words-1;
    while(top>=0) {
      if(!row[top]) {top--;continue;}
      int k=64*top+63-__builtin_clzll(row[top]);
      if(pivot[k]<0) {pivot[k]=i;answer++;break;}
      const auto& pr=rows[pivot[k]];
      for(int w=0;w<=top;w++) row[w]^=pr[w];
    }
  }
  return answer;
}

int main(int argc,char** argv) {
  assert(argc==5);
  int t=std::stoi(argv[1]);assert(t>=1&&t<=7);
  std::vector<int> boundaries;
  if(std::string(argv[2])=="all") {
    for(int b=0;b<256;b++) if(__builtin_parity(unsigned(b&63))) boundaries.push_back(b);
  } else {
    int b=std::stoi(argv[2]);assert(b>=0&&b<256&&__builtin_parity(unsigned(b&63)));
    boundaries.push_back(b);
  }
  std::vector<int> coefficients;
  if(std::string(argv[3])=="span") {
    for(int c=1;c<(1<<(t+1));c++) coefficients.push_back(c);
  } else {
    assert(std::string(argv[3])=="coeff");
    for(int c=0;c<=t;c++) coefficients.push_back(1<<c);
  }
  int n=1<<t;
  std::vector<unsigned> dilated(n);
  std::vector<std::vector<unsigned>> weight(t+1);
  for(int w=0;w<n;w++) {
    weight[__builtin_popcount(unsigned(w))].push_back(w);
    for(int i=0;i<t;i++) if((w>>i)&1) dilated[w]|=1u<<(2*i);
  }
  std::vector<std::array<std::vector<uint64_t>,2>> ranks(boundaries.size());
  for(auto& r:ranks) for(auto& cut:r) cut.assign(coefficients.size(),0);
  uint64_t queries=0;
  for(int cut=0;cut<2;cut++) for(int a=0;a<=t;a++) for(int b=0;b<=t;b++) {
    int na=weight[a].size(),nb=weight[b].size(),size=na*nb;
    std::vector<uint64_t> raw(size_t(size)*size);
    for(int i=0;i<size;i++) for(int j=0;j<size;j++) {
      unsigned x=weight[a][i/nb],v=weight[a][j/nb];
      unsigned middle_row=weight[b][i%nb],middle_col=weight[b][j%nb];
      unsigned y=cut?middle_col:middle_row,u=cut?middle_row:middle_col;
      raw[size_t(i)*size+j]=entry(t,dilated[x]+2*dilated[y],dilated[v]+2*dilated[u]);
      queries++;
    }
    for(size_t bi=0;bi<boundaries.size();bi++) {
      unsigned boundary=boundaries[bi];
      std::vector<unsigned char> polynomial(raw.size());
      for(size_t ij=0;ij<raw.size();ij++) {
        uint64_t e=raw[ij];unsigned char p=0;
        for(int k=0;k<8;k++) if((boundary>>k)&1) p^=(e>>(8*k))&255;
        polynomial[ij]=p;
      }
      for(size_t ci=0;ci<coefficients.size();ci++) {
        unsigned char table[256];
        for(int p=0;p<256;p++) table[p]=__builtin_parity(unsigned(p&coefficients[ci]));
        int words=(size+63)/64;
        std::vector<std::vector<uint64_t>> rows(size,std::vector<uint64_t>(words));
        for(int i=0;i<size;i++) for(int j=0;j<size;j++)
          if(table[polynomial[size_t(i)*size+j]]) rows[i][j>>6]|=uint64_t(1)<<(j&63);
        ranks[bi][cut][ci]+=binary_rank(rows,size);
      }
    }
  }
  std::ofstream out(argv[4]);assert(out);
  out<<"{\"t\":"<<t<<",\"queries\":"<<queries<<",\"coefficient_masks\":[";
  for(size_t i=0;i<coefficients.size();i++)out<<(i?",":"")<<coefficients[i];
  out<<"],\"profiles\":[";
  for(size_t bi=0;bi<boundaries.size();bi++) {
    out<<(bi?",":"")<<"{\"boundary\":"<<boundaries[bi]<<",\"ranks\":[";
    for(int cut=0;cut<2;cut++) {
      out<<(cut?",":"")<<"[";
      for(size_t ci=0;ci<coefficients.size();ci++)out<<(ci?",":"")<<ranks[bi][cut][ci];
      out<<"]";
    }
    out<<"]}";
  }
  out<<"]}\n";
  std::cout<<"t="<<t<<" boundaries="<<boundaries.size()<<" functionals="<<coefficients.size()<<" queries="<<queries<<"\n";
}
