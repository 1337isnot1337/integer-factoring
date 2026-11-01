#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <functional>
#include <iomanip>
#include <iostream>
#include <unordered_map>
#include <vector>
using Poly=unsigned __int128;
using Label=std::vector<int>;
int lowbit(Poly a){uint64_t b=a;return b?__builtin_ctzll(b):64+__builtin_ctzll(uint64_t(a>>64));}
Poly product(Poly a,Poly b){Poly c=0;while(b){c^=a<<lowbit(b);b&=b-1;}return c;}
Poly binomial_poly(int n){Poly a=0;int s=n;while(true){a|=Poly(1)<<s;if(!s)return a;s=(s-1)&n;}}
int rank(std::vector<uint64_t> a,int n,int words){std::vector<int> piv(n,-1);int ans=0;
 for(int i=0;i<n;++i){auto r=&a[size_t(i)*words];int w=words-1;while(w>=0){while(w>=0&&!r[w])--w;if(w<0)break;
 int c=64*w+63-__builtin_clzll(r[w]);if(piv[c]<0){piv[c]=i;++ans;break;}auto p=&a[size_t(piv[c])*words];for(int j=0;j<=w;++j)r[j]^=p[j];}}return ans;}
int main(int argc,char**argv){int h=argc>1?std::stoi(argv[1]):16;int dilation=argc>2?std::stoi(argv[2]):1;
 if(h<1||h>127||dilation<1||dilation>127||h*(1+dilation)>127)return 2;auto started=std::chrono::steady_clock::now();
 std::vector<int>w;for(int j=1;j<=h;j*=2)w.push_back(j);int r=w.size();
 std::vector<Label> labels;Label cur(r);
 std::function<void(int,int)> gen=[&](int j,int n){if(j==0){cur[0]=n;labels.push_back(cur);return;}for(int a=0;a*w[j]<=n;++a){cur[j]=a;gen(j-1,n-a*w[j]);}};gen(r-1,h);
 int d=labels.size();std::vector<int> shifts(r);int shift=0;for(int j=0;j<r;++j){shifts[j]=shift;int bits=1;while((1<<bits)<=h/w[j])++bits;shift+=bits;}
 assert(shift<=63);auto key=[&](const Label&a){uint64_t v=0;for(int j=0;j<r;++j)v|=uint64_t(a[j])<<shifts[j];return v;};
 std::unordered_map<uint64_t,int> ids;for(int i=0;i<d;++i)ids[key(labels[i])]=i;
 std::vector<std::vector<uint32_t>> blocks(1<<r);
 for(int a=0;a<d;++a)for(int b=0;b<d;++b){int parity=0;for(int j=0;j<r;++j)parity|=((labels[a][j]+labels[b][j])&1)<<j;blocks[parity].push_back(a*d+b);}
 size_t source_size=size_t(d)*d*d;std::vector<Poly> cache(source_size);std::vector<uint8_t> known(source_size);
 auto source=[&](int a,int y,int beta){size_t index=(size_t(a)*d+y)*d+beta;if(known[index])return cache[index];Poly answer=1;
   for(int j=0;j<r;++j){Poly local=0;for(int k=0;k<=std::min(labels[a][j],labels[y][j]);++k){int v=labels[a][j]-k;
    if((k&~labels[y][j])==0&&v<=labels[beta][j]&&(v&~labels[beta][j])==0)local^=binomial_poly(w[j]*k);}
    answer=product(answer,local);if(!answer)break;}
   known[index]=1;return cache[index]=answer;};
 auto value=[&](int a,int y,int u,int beta){Label gamma(r),v(r);for(int j=0;j<r;++j){int sigma=labels[beta][j]+labels[u][j],delta=labels[a][j]-labels[y][j];
   if(sigma<std::abs(delta)||((sigma+delta)&1))return Poly(0);gamma[j]=(sigma+delta)/2;v[j]=(sigma-delta)/2;}
   auto i=ids.find(key(gamma)),j=ids.find(key(v));if(i==ids.end()||j==ids.end())return Poly(0);
   Poly A=source(a,y,i->second),B=source(i->second,u,beta),dilated=0;
   while(B){int k=lowbit(B);dilated|=Poly(1)<<(dilation*k);B&=B-1;}
   return product(A,dilated);};
 int degree=(dilation+1)*h;std::vector<std::vector<uint64_t>> ranks;uint64_t checked=0,entries=0;int maxblock=0;
 for(bool row:{true,false}){std::vector<uint64_t> total(degree+1);int done=0;
  for(const auto&pairs:blocks){int n=pairs.size();if(!n)continue;maxblock=std::max(maxblock,n);if(row)entries+=uint64_t(n)*n;int words=(n+63)/64;size_t plane=size_t(n)*words;
   std::vector<std::vector<uint64_t>> matrices(degree+1,std::vector<uint64_t>(plane));
   for(int i=0;i<n;++i)for(int j=0;j<n;++j){int a=pairs[i]/d,second=pairs[i]%d,beta=pairs[j]/d,last=pairs[j]%d;
    int y=row?second:last,u=row?last:second;Poly p=value(a,y,u,beta);int parity=__builtin_parityll(uint64_t(p))^__builtin_parityll(uint64_t(p>>64));assert(parity==(i==j));++checked;
    while(p){int k=lowbit(p);matrices[k][size_t(i)*words+j/64]|=uint64_t(1)<<(j%64);p&=p-1;}}
   for(int k=0;k<=degree;++k)total[k]+=rank(std::move(matrices[k]),n,words);
   std::cerr<<(row?"row":"col")<<" block "<<++done<<" n="<<n<<"\n";
  }ranks.push_back(total);}
 long double S=0;for(int k=0;k<=degree;++k)S+=sqrtl((long double)ranks[0][k]*ranks[1][k]);
 std::cout<<std::setprecision(18)<<"{\"energy\":"<<h<<",\"labels\":"<<d<<",\"dilation\":"<<dilation<<",\"entries_per_cut\":"<<entries<<",\"max_block\":"<<maxblock<<",\"crossing_checks\":"<<checked<<",\"row_ranks\":[";
 for(int k=0;k<=degree;++k)std::cout<<(k?",":"")<<ranks[0][k];std::cout<<"],\"col_ranks\":[";for(int k=0;k<=degree;++k)std::cout<<(k?",":"")<<ranks[1][k];
 std::cout<<"],\"rank_sum\":"<<S<<",\"threshold\":"<<powl(d,sqrtl(2))<<",\"rho\":"<<logl(S)/(2*logl(d))<<",\"seconds\":"<<std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count()<<"}\n";
}
