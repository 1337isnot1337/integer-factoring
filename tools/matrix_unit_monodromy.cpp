#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <string>
#include <vector>

// Auxiliary order four.  K(z;p,q)=[p=q] I+(1+z) E_(q,p),
// with p=(x,y), q=(v,u), and boundary ell(M)=M_00.
// Each polynomial entry is evaluated by its exact auxiliary walk.
uint32_t one_entry(int t,uint32_t p,uint32_t q,int initial) {
  int state=initial,diagonal=0,offdiagonal=0;
  for(int i=0;i<t;i++,p>>=2,q>>=2) {
    int a=p&3,b=q&3;
    if(a==b) { if(a==state) diagonal++; }
    else {
      if(b!=state) return 0;
      state=a; offdiagonal++;
    }
  }
  if(state!=initial) return 0;
  uint32_t answer=0;
  for(int sub=offdiagonal;;sub=(sub-1)&offdiagonal) {
    answer|=uint32_t(1)<<(diagonal+sub);
    if(!sub) break;
  }
  return answer;
}

uint32_t particle_entry(int t,uint32_t p,uint32_t q,unsigned initial,bool coherent) {
  unsigned state=initial;
  int diagonal=0,offdiagonal=0;
  for(int i=0;i<t;i++,p>>=2,q>>=2) {
    int a=p&3,b=q&3;
    if(a==b) {if((state>>a)&1)diagonal++;}
    else {
      if(!((state>>b)&1)||((state>>a)&1))return 0;
      state^=(1u<<a)|(1u<<b);offdiagonal++;
    }
  }
  if(state!=initial&&(!coherent||state!=(initial^15)))return 0;
  uint32_t answer=0;
  for(int sub=offdiagonal;;sub=(sub-1)&offdiagonal) {
    answer|=uint32_t(1)<<(diagonal+sub);
    if(!sub)break;
  }
  return answer;
}

uint32_t entry(int t,uint32_t p,uint32_t q,unsigned boundary,bool particles,bool coherent) {
  if(particles)return particle_entry(t,p,q,boundary,coherent);
  uint32_t answer=0;
  for(int a=0;a<4;a++) if((boundary>>a)&1)answer^=one_entry(t,p,q,a);
  return answer;
}

uint64_t rank(std::vector<std::vector<uint64_t>>& rows,int n) {
  int words=(n+63)/64;
  std::vector<int> pivots(n,-1);
  uint64_t answer=0;
  for(int i=0;i<n;i++) {
    auto& row=rows[i];
    int top=words-1;
    while(top>=0) {
      if(!row[top]) {top--;continue;}
      int p=64*top+63-__builtin_clzll(row[top]);
      if(pivots[p]<0) {pivots[p]=i;answer++;break;}
      const auto& pivot=rows[pivots[p]];
      for(int w=0;w<=top;w++) row[w]^=pivot[w];
    }
  }
  return answer;
}

template<class Evaluate>
void add_block(int t,int n,Evaluate evaluate,std::vector<uint64_t>& ranks,
               uint64_t& entry_queries) {
  std::vector<uint32_t> polynomials(size_t(n)*n);
  for(int i=0;i<n;i++) for(int j=0;j<n;j++) {
    polynomials[size_t(i)*n+j]=evaluate(i,j); entry_queries++;
  }
  int words=(n+63)/64;
  for(int d=0;d<=t;d++) {
    std::vector<std::vector<uint64_t>> rows(n,std::vector<uint64_t>(words));
    bool nonzero=false;
    for(int i=0;i<n;i++) for(int j=0;j<n;j++)
      if((polynomials[size_t(i)*n+j]>>d)&1) {
        rows[i][j>>6]|=uint64_t(1)<<(j&63);nonzero=true;
      }
    if(nonzero) ranks[d]+=rank(rows,n);
  }
}

int main(int argc,char** argv) {
  assert(argc>=3&&argc<=5);
  int t=std::stoi(argv[1]); assert(1<=t&&t<=9);
  unsigned boundary=argc==4?std::stoul(argv[3]):1;
  if(argc==5)boundary=std::stoul(argv[3]);
  bool coherent=argc==5&&std::string(argv[4])=="coherent";
  bool particles=argc==5&&(std::string(argv[4])=="particles"||coherent);
  assert(boundary<16 && (particles||__builtin_parity(boundary)));
  if(coherent)assert(boundary==6||boundary==9);
  int n=1<<t;
  std::vector<uint64_t> rr(t+1),cc(t+1);
  std::map<uint32_t,std::vector<uint32_t>> rowblocks;
  for(uint32_t word=0;word<uint32_t(n)*n;word++) {
    uint32_t key=0,copy=word;
    for(int i=0;i<t;i++,copy>>=2) key+=1u<<(4*(copy&3));
    rowblocks[key].push_back(word);
  }
  uint64_t rowqueries=0,colqueries=0;
  if(!coherent)for(const auto& block:rowblocks) {
    const auto& labels=block.second;
    add_block(t,labels.size(),[&](int i,int j){return entry(t,labels[i],labels[j],boundary,particles,coherent);},rr,rowqueries);
  }
  std::vector<std::vector<uint32_t>> weight(t+1);
  std::vector<uint32_t> dilated(n);
  for(int w=0;w<n;w++) {
    weight[__builtin_popcount(unsigned(w))].push_back(w);
    for(int i=0;i<t;i++) if((w>>i)&1) dilated[w]|=1u<<(2*i);
  }
  for(int a=0;a<=t;a++) for(int b=0;b<=t;b++) {
    int na=weight[a].size(),nb=weight[b].size(),size=na*nb;
    add_block(t,size,[&](int i,int j) {
      uint32_t x=weight[a][i/nb],u=weight[b][i%nb];
      uint32_t v=weight[a][j/nb],y=weight[b][j%nb];
      return entry(t,dilated[x]+2*dilated[y],dilated[v]+2*dilated[u],boundary,particles,coherent);
    },cc,colqueries);
    if(coherent)add_block(t,size,[&](int i,int j) {
      uint32_t x=weight[a][i/nb],y=weight[b][i%nb];
      uint32_t v=weight[a][j/nb],u=weight[b][j%nb];
      return entry(t,dilated[x]+2*dilated[y],dilated[v]+2*dilated[u],boundary,particles,coherent);
    },rr,rowqueries);
  }
  long double sum=0;uint64_t maxproduct=0,integerlower=0;
  for(int d=0;d<=t;d++) {
    uint64_t p=rr[d]*cc[d]; maxproduct=std::max(maxproduct,p);
    sum+=std::sqrt((long double)p);
    uint64_t r=std::sqrt((long double)p);
    while((r+1)*(r+1)<=p)r++;
    while(r*r>p)r--;
    integerlower+=r;
  }
  auto output=[&](std::ostream& out) {
    out<<"{\"t\":"<<t<<",\"external_labels\":"<<n<<",\"auxiliary_order\":"<<(particles?std::vector<int>{1,4,6,4,1}[__builtin_popcount(boundary)]:4)
       <<",\"particles\":"<<(particles?"true":"false")<<",\"coherent\":"<<(coherent?"true":"false")<<",\"boundary_mask\":"<<boundary<<",\"row_ranks\":[";
    for(int d=0;d<=t;d++)out<<(d?",":"")<<rr[d];
    out<<"],\"column_ranks\":[";
    for(int d=0;d<=t;d++)out<<(d?",":"")<<cc[d];
    out<<"],\"maximum_rank_product\":"<<maxproduct<<",\"rank_sum_estimate\":"<<std::setprecision(18)<<sum
       <<",\"integer_rank_sum_lower_bound\":"<<integerlower
       <<",\"maximum_product_rho\":"<<std::log2((long double)maxproduct)/(4*t)
       <<",\"rank_sum_rho\":"<<std::log2(sum)/(2*t)
       <<",\"row_blocks\":"<<(coherent?(t+1)*(t+1):rowblocks.size())<<",\"column_blocks\":"<<(t+1)*(t+1)
       <<",\"row_entry_queries\":"<<rowqueries<<",\"column_entry_queries\":"<<colqueries<<"}\n";
  };
  std::ofstream out(argv[2]); assert(out);output(out);output(std::cout);
}
