#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

unsigned degree(uint32_t a) {return a?31-__builtin_clz(a):0;}
uint32_t remainder(uint32_t a,uint32_t b) {
  while(a && degree(a)>=degree(b))a^=b<<(degree(a)-degree(b));
  return a;
}
bool irreducible(uint32_t p,unsigned t) {
  for(unsigned k=1;k<=t/2;++k)
    for(uint32_t q=(1u<<k)|1u;q<(1u<<(k+1));q+=2)
      if(!remainder(p,q))return false;
  return true;
}
struct Field {
  unsigned t,d,p,g;
  std::vector<unsigned> logarithm,exponential;
  unsigned slow_mul(unsigned a,unsigned b) const {
    unsigned c=0;
    while(b){if(b&1)c^=a;b>>=1;a<<=1;if(a&d)a^=p;}
    return c;
  }
  unsigned power(unsigned a,unsigned n) const {
    unsigned c=1;
    while(n){if(n&1)c=slow_mul(c,a);a=slow_mul(a,a);n>>=1;}
    return c;
  }
  unsigned mul(unsigned a,unsigned b) const {
    return (!a||!b)?0:exponential[logarithm[a]+logarithm[b]];
  }
  unsigned trace(unsigned a) const {
    unsigned c=0;for(unsigned i=0;i<t;++i){c^=a;a=mul(a,a);}assert(c<2);return c;
  }
  explicit Field(unsigned width):t(width),d(1u<<width) {
    for(p=d|1;!irreducible(p,t);p+=2){}
    unsigned order=d-1,n=order;std::vector<unsigned> primes;
    for(unsigned r=2;r*r<=n;++r)if(n%r==0){primes.push_back(r);while(n%r==0)n/=r;}
    if(n>1)primes.push_back(n);
    for(g=2;g<d;++g){bool ok=true;for(auto r:primes)if(power(g,order/r)==1)ok=false;if(ok)break;}
    assert(g<d);
    logarithm.assign(d,~0u);exponential.resize(2*order);
    unsigned value=1;
    for(unsigned i=0;i<order;++i){assert(logarithm[value]==~0u);logarithm[value]=i;exponential[i]=value;value=slow_mul(value,g);}
    assert(value==1);
    for(unsigned i=0;i<order;++i)exponential[i+order]=exponential[i];
  }
};

uint64_t reverse_bits(uint64_t x) {
  x=((x>>1)&0x5555555555555555ULL)|((x&0x5555555555555555ULL)<<1);
  x=((x>>2)&0x3333333333333333ULL)|((x&0x3333333333333333ULL)<<2);
  x=((x>>4)&0x0f0f0f0f0f0f0f0fULL)|((x&0x0f0f0f0f0f0f0f0fULL)<<4);
  x=((x>>8)&0x00ff00ff00ff00ffULL)|((x&0x00ff00ff00ff00ffULL)<<8);
  x=((x>>16)&0x0000ffff0000ffffULL)|((x&0x0000ffff0000ffffULL)<<16);
  return (x>>32)|(x<<32);
}

int main(int argc,char**argv) {
  for(unsigned i=0;i<64;++i)assert(reverse_bits(1ULL<<i)==(1ULL<<(63-i)));
  if(argc!=5){std::cerr<<"usage: trace_weight_restriction t weight maximum_k output\n";return 2;}
  unsigned t=std::stoul(argv[1]),weight=std::stoul(argv[2]),max_k=std::stoul(argv[3]);
  if(t<2||t>18||weight>t||max_k<1||max_k>6||max_k>=t)return 2;
  Field field(t);unsigned d=field.d;
  std::vector<unsigned> labels,linear(max_k),sigma(d),squares(d);
  for(unsigned x=0;x<d;++x)if(unsigned(__builtin_popcount(x))==weight)labels.push_back(x);
  unsigned a=1;
  for(unsigned j=0;j<max_k;++j){
    for(unsigned i=0;i<t;++i)linear[j]|=field.trace(field.mul(a,1u<<i))<<i;
    a=field.mul(a,2);
  }
  for(unsigned x=0;x<d;++x){
    squares[x]=field.mul(x,x);
    for(unsigned j=0;j<max_k;++j)sigma[x]|=(__builtin_parity(x&linear[j])<<j);
  }
  std::vector<uint64_t> images(d);
  std::vector<unsigned> block_sizes(d);
  images[0]=1;block_sizes[0]=labels.size();
  for(unsigned i=0;i<labels.size();++i){
    unsigned x=labels[i];
    for(unsigned j=i+1;j<labels.size();++j){
      unsigned y=labels[j],s=x^y,xy=field.mul(x,y);
      images[s]|=1ULL<<sigma[squares[x]^xy];
      images[s]|=1ULL<<sigma[squares[y]^xy];
      block_sizes[s]+=2;
    }
  }
  std::vector<uint64_t> ranks(max_k+1),caps(max_k+1),image_totals(max_k+1),blocks(max_k+1);
  for(unsigned s=0;s<d;++s){
    uint64_t image=images[s];
    for(unsigned k=max_k;k>0;--k){
      unsigned width=1u<<k;
      unsigned paired=__builtin_popcountll(image&(reverse_bits(image)>>(64-width)));
      ranks[k]+=paired;caps[k]+=std::min(block_sizes[s],width);
      image_totals[k]+=__builtin_popcountll(image);blocks[k]+=paired!=0;
      unsigned half=width/2;uint64_t mask=(1ULL<<half)-1;
      image=(image&mask)|(image>>half);
    }
  }
  std::string path=argv[4],temporary=path+".tmp";
  std::ofstream out(temporary);out.exceptions(std::ios::badbit|std::ios::failbit);
  out<<"{\"t\":"<<t<<",\"weight\":"<<weight<<",\"label_count\":"<<labels.size()
     <<",\"field_polynomial\":"<<field.p<<",\"primitive_element\":"<<field.g
     <<",\"unordered_label_pairs\":"<<uint64_t(labels.size())*(labels.size()-1)/2
     <<",\"profiles\":[";
  for(unsigned k=1;k<=max_k;++k){
    if(k>1)out<<",";
    out<<"{\"k\":"<<k<<",\"row_rank\":"<<ranks[k]<<",\"column_rank\":"<<ranks[k]
       <<",\"block_size_cap\":"<<caps[k]<<",\"syndrome_image_total\":"<<image_totals[k]
       <<",\"nonzero_blocks\":"<<blocks[k]<<"}";
    std::cout<<"t "<<t<<" weight "<<weight<<" labels "<<labels.size()<<" k "<<k
             <<" exact rank "<<ranks[k]<<" cap "<<caps[k];
    if(ranks[k] && labels.size()>1)std::cout<<" exponent "<<std::log(double(ranks[k]))/std::log(double(labels.size()));
    std::cout<<"\n";
  }
  out<<"]}\n";out.flush();out.close();
  if(std::rename(temporary.c_str(),path.c_str()))return 3;
}
