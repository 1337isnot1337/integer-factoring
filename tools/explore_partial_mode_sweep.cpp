#include <algorithm>
#include <array>
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <functional>
#include <iomanip>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

using Label=std::vector<int>;
using Blocks=std::unordered_map<uint64_t,std::vector<uint32_t>>;
using Poly=std::array<uint64_t,4>;
bool odd_choose(int n,int k){return k>=0&&k<=n&&(k&~n)==0;}
bool nonzero(const Poly&p){return p[0]||p[1]||p[2]||p[3];}
Poly multiply64(const Poly&a,uint64_t b){Poly out{};while(b){int shift=__builtin_ctzll(b);for(int j=0;j<4;++j){out[j]^=a[j]<<shift;if(shift&&j<3)out[j+1]^=a[j]>>(64-shift);}b&=b-1;}return out;}
int rank(std::vector<uint64_t> a,int n,int words){std::vector<int> piv(n,-1);int ans=0;
 for(int i=0;i<n;++i){auto r=&a[size_t(i)*words];int w=words-1;while(w>=0){while(w>=0&&!r[w])--w;if(w<0)break;
 int c=64*w+63-__builtin_clzll(r[w]);if(piv[c]<0){piv[c]=i;++ans;break;}auto p=&a[size_t(piv[c])*words];for(int j=0;j<=w;++j)r[j]^=p[j];}}return ans;}

int main(int argc,char**argv){int h=argc>1?std::stoi(argv[1]):27;std::string order_arg=argc>2?argv[2]:"ascending";bool constant_only=argc>3&&std::string(argv[3])=="--constant-only";
 if(h<3||h>63)return 2;auto started=std::chrono::steady_clock::now();
 std::vector<int>w;for(int a=1;a<=h;a*=3)w.push_back(a);int r=w.size();
 std::vector<int>order;for(int j=0;j<r-1;++j)order.push_back(j);
 if(order_arg=="descending")std::reverse(order.begin(),order.end());
 else if(order_arg!="ascending"){if(order_arg.size()!=order.size())return 2;for(int j=0;j<r-1;++j)order[j]=order_arg[j]-'0';auto sorted=order;std::sort(sorted.begin(),sorted.end());for(int j=0;j<r-1;++j)if(sorted[j]!=j)return 2;}
 int degree=constant_only?0:(r-1)*h;assert(degree<256);
 std::vector<Label>labels;Label cur(r);
 std::function<void(int,int)>gen=[&](int j,int remaining){if(j==0){cur[0]=remaining;labels.push_back(cur);return;}for(int a=0;a*w[j]<=remaining;++a){cur[j]=a;gen(j-1,remaining-a*w[j]);}};gen(r-1,h);int d=labels.size();
 std::vector<int>shifts(r);int shift=0;for(int j=0;j<r;++j){shifts[j]=shift;int bits=1;while((1<<bits)<=2*h/w[j])++bits;shift+=bits;}assert(shift<=64);
 Blocks row_blocks,col_blocks;
 for(int a=0;a<d;++a)for(int b=0;b<d;++b){uint64_t sum=0,difference=0;for(int j=0;j<r;++j){sum|=uint64_t(labels[a][j]+labels[b][j])<<shifts[j];difference|=uint64_t(labels[a][j]-labels[b][j]+h/w[j])<<shifts[j];}row_blocks[difference].push_back(a*d+b);col_blocks[sum].push_back(a*d+b);}
 std::vector<std::vector<uint64_t>>cache(r);std::vector<std::vector<uint8_t>>known(r);std::vector<int>width(r);for(int j=0;j<r;++j){width[j]=h/w[j]+1;size_t n=size_t(width[j])*width[j]*width[j];cache[j].resize(n);known[j].resize(n);}
 std::vector<uint64_t>binomial(h+1);for(int n=0;n<=h;++n){int s=n;while(true){binomial[n]|=uint64_t(1)<<s;if(!s)break;s=(s-1)&n;}}
 auto mode=[&](int j,int a,int y,int beta){size_t index=(size_t(a)*width[j]+y)*width[j]+beta;if(known[j][index])return cache[j][index];uint64_t answer=0;
  if(constant_only)answer=odd_choose(beta+y,a);
  else for(int k=0;k<=std::min(a,y);++k)if(odd_choose(y,k)&&odd_choose(beta,a-k))answer^=binomial[w[j]*k];
  known[j][index]=1;return cache[j][index]=answer;};
 uint64_t endpoint_checks=0,entries=0,valid_paths=0;
 auto entry=[&](int ai,int yi,int bi,int ui){const auto&alpha=labels[ai];const auto&y=labels[yi];Label first=labels[bi],second=labels[ui];std::vector<int>q(r-1);int prior=0;
  for(int j=0;j<r-1;++j){int numerator=alpha[j]-first[j]+prior;assert(numerator%3==0);q[j]=numerator/3;prior=q[j];}
  Poly answer{1,0,0,0};bool exists=true;
  for(int j:order){int old0=first[j],old1=first[j+1];first[j]+=3*q[j];first[j+1]-=q[j];second[j]+=3*q[j];second[j+1]-=q[j];
   if(first[j]<0||first[j+1]<0||second[j]<0||second[j+1]<0){exists=false;answer={};break;}
   uint64_t left=mode(j,first[j],second[j],old0),right=mode(j+1,first[j+1],second[j+1],old1);
   Poly gate=multiply64(Poly{left,0,0,0},right);assert(!gate[1]&&!gate[2]&&!gate[3]);answer=multiply64(answer,gate[0]);
  }
  if(exists){assert(first==alpha&&second==y);++valid_paths;}
  if(!constant_only){int parity=0;for(uint64_t v:answer)parity^=__builtin_parityll(v);assert(parity==(ai==bi&&yi==ui));++endpoint_checks;}
  return answer;};
 std::vector<std::vector<uint64_t>>ranks;
 for(bool rowcut:{true,false}){const auto&blocks=rowcut?row_blocks:col_blocks;std::vector<uint64_t>total(degree+1);
  for(const auto&item:blocks){const auto&pairs=item.second;int n=pairs.size(),words=(n+63)/64;size_t plane=size_t(n)*words;
   std::vector<std::vector<uint64_t>>matrices(degree+1,std::vector<uint64_t>(plane));
   for(int i=0;i<n;++i)for(int j=0;j<n;++j){int a=pairs[i]/d,y=rowcut?pairs[i]%d:pairs[j]/d,beta=rowcut?pairs[j]/d:pairs[j]%d,u=rowcut?pairs[j]%d:pairs[i]%d;
    Poly p=entry(a,y,beta,u);for(int word=0;word<4;++word){uint64_t v=p[word];while(v){int k=64*word+__builtin_ctzll(v);assert(k<=degree);matrices[k][size_t(i)*words+j/64]|=uint64_t(1)<<(j%64);v&=v-1;}}++entries;}
   for(int k=0;k<=degree;++k)total[k]+=rank(std::move(matrices[k]),n,words);
  }ranks.push_back(total);
 }
 long double S=0;int pieces=0;for(int k=0;k<=degree;++k){assert((ranks[0][k]==0)==(ranks[1][k]==0));if(ranks[0][k])++pieces;S+=sqrtl((long double)ranks[0][k]*ranks[1][k]);}
 std::cout<<std::setprecision(18)<<"{\"field\":2,\"energy\":"<<h<<",\"labels\":"<<d<<",\"application_order\":[";for(size_t j=0;j<order.size();++j)std::cout<<(j?",":"")<<order[j];
 std::cout<<"],\"constant_only\":"<<(constant_only?"true":"false")<<",\"degree_bound\":"<<degree<<",\"nonzero_pieces\":"<<pieces<<",\"entry_checks\":"<<entries<<",\"endpoint_checks\":"<<endpoint_checks<<",\"valid_paths\":"<<valid_paths<<",\"row_ranks\":[";
 for(int k=0;k<=degree;++k)std::cout<<(k?",":"")<<ranks[0][k];std::cout<<"],\"column_ranks\":[";for(int k=0;k<=degree;++k)std::cout<<(k?",":"")<<ranks[1][k];
 std::cout<<"],\"rank_sum\":"<<S<<",\"finite_target\":"<<powl(d,sqrtl(2))<<",\"seconds\":"<<std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count()<<"}\n";
 return 0;
}
