// Exact width-five binary-code pencil screen; ordinary finite certificate.
// No factoring improvement or asymptotic rank obstruction is asserted.
#include <array>
#include <algorithm>
#include <cassert>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <random>
#include <string>
#include <vector>
#include <omp.h>

using U=uint32_t;
using Kernels=std::array<std::array<U,32>,6>;
std::array<int,10> qm,cm;
uint8_t qr[1024],cr[1024],mul3[256][256],split[10][4][4][256];
uint16_t qe[4][256],ce[4][256];
U linear_mask=0,quartic_mask=0;

inline U shift(U x,int p){
 if(p&1)x=((x&0x55555555U)<<1)|((x>>1)&0x55555555U);
 if(p&2)x=((x&0x33333333U)<<2)|((x>>2)&0x33333333U);
 if(p&4)x=((x&0x0f0f0f0fU)<<4)|((x>>4)&0x0f0f0f0fU);
 if(p&8)x=((x&0x00ff00ffU)<<8)|((x>>8)&0x00ff00ffU);
 if(p&16)x=(x<<16)|(x>>16);
 return x;
}
inline U anf(U x){
 x^=(x>>1)&0x55555555U; x^=(x>>2)&0x33333333U;
 x^=(x>>4)&0x0f0f0f0fU; x^=(x>>8)&0x00ff00ffU;
 x^=(x>>16)&0x0000ffffU; return x;
}
int alt_rank(int pattern,const std::array<int,10>& masks,bool complement){
 U rows[5]={};
 for(int j=0;j<10;++j)if(pattern>>j&1){
  int m=complement?31^masks[j]:masks[j];
  int a=__builtin_ctz(m),b=__builtin_ctz(m^(1<<a));
  rows[a]^=1U<<b;rows[b]^=1U<<a;
 }
 int rank=0;
 for(int c=0;c<5;++c){
  int p=rank;while(p<5&&!(rows[p]>>c&1))++p;
  if(p==5)continue;std::swap(rows[p],rows[rank]);
  for(int i=rank+1;i<5;++i)if(rows[i]>>c&1)rows[i]^=rows[rank];
  ++rank;
 }
 return rank;
}
void initialize(){
 int q=0,c=0;
 for(int m=0;m<32;++m){
  int k=__builtin_popcount((U)m);
  if(k==1)linear_mask|=1U<<m;
  if(k==2)qm[q++]=m;
  if(k==3)cm[c++]=m;
  if(k==4)quartic_mask|=1U<<m;
 }
 for(int p=0;p<1024;++p){qr[p]=alt_rank(p,qm,false);cr[p]=alt_rank(p,cm,true);}
 for(int byte=0;byte<4;++byte)for(int value=0;value<256;++value){
  for(int j=0;j<10;++j){
   if(qm[j]/8==byte&&(value>>(qm[j]%8)&1))qe[byte][value]|=1<<j;
   if(cm[j]/8==byte&&(value>>(cm[j]%8)&1))ce[byte][value]|=1<<j;
  }
 }
 for(int a=0;a<256;++a)for(int b=0;b<256;++b){
  int result=0;
  for(int i=0;i<8;++i)if(a>>i&1)for(int j=0;j<8;++j)
   if((b>>j&1)&&!(i&j))result^=1<<(i|j);
  mul3[a][b]=result;
 }
 for(int pair=0;pair<10;++pair){
  int bit_i=qm[pair]&-qm[pair],bit_j=qm[pair]^bit_i;
  std::vector<int> remaining;
  for(int b=1;b<32;b<<=1)if(b!=bit_i&&b!=bit_j)remaining.push_back(b);
  for(int flag=0;flag<4;++flag)for(int small=0;small<8;++small){
   int original=(flag&1?bit_i:0)|(flag&2?bit_j:0);
   for(int j=0;j<3;++j)if(small>>j&1)original|=remaining[j];
   int byte=original/8,bit=original%8;
   for(int value=0;value<256;++value)if(value>>bit&1)
    split[pair][flag][byte][value]|=1<<small;
  }
 }
}
inline int fast_rank(U kernel){
 U a=anf(kernel);
 if(a&1)return 32;
 if(a&linear_mask)return 16;
 int q=qe[0][a&255]|qe[1][a>>8&255]|qe[2][a>>16&255]|qe[3][a>>24];
 if(q){
  if(qr[q]==4)return 12;
  int pair=__builtin_ctz((U)q),s[4];
  for(int f=0;f<4;++f)s[f]=split[pair][f][0][a&255]^split[pair][f][1][a>>8&255]
    ^split[pair][f][2][a>>16&255]^split[pair][f][3][a>>24];
  int residual=s[0]^mul3[mul3[s[1]][s[2]]][s[3]];
  assert(residual==0||residual==128);
  return 8+(residual?2:0);
 }
 int c=ce[0][a&255]|ce[1][a>>8&255]|ce[2][a>>16&255]|ce[3][a>>24];
 if(c)return 2+cr[c];
 if(a&quartic_mask)return 2;
 return a?1:0;
}
int direct_rank(U kernel){
 U basis[32]={};int rank=0;
 for(int x=0;x<32;++x){
  U row=shift(kernel,x);
  while(row){int p=31-__builtin_clz(row);if(basis[p])row^=basis[p];
   else{basis[p]=row;++rank;break;}}
 }
 return rank;
}
void transpose(std::array<U,32>& a){
 U mask=0x0000ffffU;
 for(int j=16;j;j>>=1,mask^=mask<<j){
  for(int k=0;k<32;k=(k+j+1)&~j){
   U t=((a[k]>>j)^a[k+j])&mask;
   a[k]^=t<<j;a[k+j]^=t;
  }
 }
}
std::pair<int,int> piece_ranks(const std::array<U,32>& rows){
 int r=0,c=0;
 for(U a:rows)r+=fast_rank(a);
 auto cols=rows;transpose(cols);
 for(U a:cols)c+=fast_rank(a);
 return {r,c};
}
Kernels base(){
 Kernels k{};
 for(int s=0;s<32;++s)for(int w=0;w<32;++w)if(!(s&~w)){
  int a=__builtin_popcount((U)s),b=__builtin_popcount((U)(w^s));
  for(int i=0;i<=b;++i)if(!(i&~b))k[a+i][s]|=1U<<w;
 }
 return k;
}
void toggle(Kernels& k,int p){
 for(int i=0;i<6;++i)for(int s=0;s<32;++s){
  U next=i<5?k[i+1][s]:0;
  k[i][s]=shift(k[i][s]^next,p)^next;
 }
}
void validate(){
 std::mt19937 gen(261030);
 for(int i=0;i<100000;++i){U x=gen();assert(fast_rank(x)==direct_rank(x));}
 for(int i=0;i<1000;++i){
  std::array<U,32>a{},b{};for(U& x:a)x=gen();
  for(int r=0;r<32;++r)for(int c=0;c<32;++c)b[c]|=((a[r]>>c)&1U)<<r;
  transpose(a);assert(a==b);
 }
 auto k=base();int br[6]={1,12,66,232,496,1024},bc[6]={1024,496,232,66,12,1};
 for(int i=0;i<6;++i)assert(piece_ranks(k[i])==std::make_pair(br[i],bc[i]));
 for(int p=1;p<32;++p){auto other=k;toggle(other,p);toggle(other,p);assert(other==k);}
 std::cerr<<"validated 100000 direct convolution ranks, 1000 transposes, baseline, 31 involutions\n";
}
int main(int argc,char**argv){
 initialize();validate();
 int first=argc>1?std::stoi(argv[1]):0,last=argc>2?std::stoi(argv[2]):1;
 int threads=argc>3?std::stoi(argv[3]):1;
 std::string prefix=argc>4?argv[4]:"quadratic-work/code-width5";
 assert(0<=first&&first<last&&last<=64&&threads>=1);
 std::vector<int> parities;for(int p=1;p<32;++p)if(__builtin_popcount((U)p)>1)parities.push_back(p);
 assert(parities.size()==26);
 auto start=std::chrono::steady_clock::now();
 #pragma omp parallel for schedule(dynamic,1) num_threads(threads)
 for(int chunk=first;chunk<last;++chunk){
  uint64_t begin=uint64_t(chunk)<<20,end=uint64_t(chunk+1)<<20;
  auto k=base();U gray=U(begin^(begin>>1));
  for(int bit=0;bit<26;++bit)if(gray>>bit&1)toggle(k,parities[bit]);
  uint64_t counts[6]={},violations=0,first_bad=0,checksum=0;
  int order[6]={2,3,1,4,0,5};
  for(uint64_t index=begin;index<end;++index){
   if(index!=begin)toggle(k,parities[__builtin_ctzll(index)]);
   bool passed=false;int maximum=0;
   for(int pos=0;pos<6;++pos){
    int i=order[pos];auto [r,c]=piece_ranks(k[i]);int q=r*c;
    maximum=std::max(maximum,q);
    if(q>=15312){++counts[i];checksum+=uint64_t(q)*(index+1);passed=true;break;}
   }
   if(!passed){if(!violations)first_bad=index;++violations;}
  }
  double elapsed=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
  std::string path=prefix+"-chunk-"+std::to_string(chunk)+".json";
  std::ofstream out(path);assert(out);
  out<<"{\"chunk\":"<<chunk<<",\"begin\":"<<begin<<",\"end\":"<<end
     <<",\"covered\":"<<(end-begin)<<",\"threshold\":15312,\"violations\":"<<violations
     <<",\"first_bad_index\":"<<first_bad<<",\"checksum_mod_2_64\":"<<checksum<<",\"witness_degree_counts\":[";
  for(int i=0;i<6;++i)out<<(i?",":"")<<counts[i];out<<"],\"elapsed_seconds\":"<<elapsed<<"}\n";out.close();
  #pragma omp critical
  std::cout<<"{\"chunk\":"<<chunk<<",\"covered\":"<<(end-begin)<<",\"violations\":"<<violations
           <<",\"first_bad_index\":"<<first_bad<<",\"elapsed_seconds\":"<<elapsed<<"}"<<std::endl;
 }
}
