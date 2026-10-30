#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

using Label=std::vector<int>;
struct Hash {size_t operator()(const Label& a)const{size_t h=0;for(int x:a)h=h*1000003+uint32_t(x+1000);return h;}};
int matrix_rank(std::vector<uint8_t> A,int n,int p){int result=0;
 for(int c=0;c<n&&result<n;++c){int pivot=result;while(pivot<n&&!A[size_t(pivot)*n+c])++pivot;if(pivot==n)continue;
  for(int j=c;j<n;++j)std::swap(A[size_t(result)*n+j],A[size_t(pivot)*n+j]);
  int inv=1;while(inv*A[size_t(result)*n+c]%p!=1)++inv;
  if(inv!=1)for(int j=c;j<n;++j)A[size_t(result)*n+j]=A[size_t(result)*n+j]*inv%p;
  for(int r=result+1;r<n;++r){int factor=A[size_t(r)*n+c];if(!factor)continue;
   for(int j=c;j<n;++j){int value=A[size_t(r)*n+j]-factor*A[size_t(result)*n+j];value%=p;if(value<0)value+=p;A[size_t(r)*n+j]=value;}}
  ++result;
 }return result;
}

int main(int argc,char**argv){
 if(argc!=4)return 2;int p=std::stoi(argv[1]),h=std::stoi(argv[2]);if(h%p)return 2;int m=h/p;
 auto started=std::chrono::steady_clock::now();std::vector<int> weights{1};for(int w=1;w<=m;w*=p)weights.push_back(w);
 int modes=weights.size();std::vector<Label> labels;Label current(modes);
 std::function<void(int,int)> generate=[&](int j,int energy){if(j==0){current[0]=energy;labels.push_back(current);return;}
  for(int a=0;a*weights[j]<=energy;++a){current[j]=a;generate(j-1,energy-a*weights[j]);}};generate(modes-1,m);
 int d=labels.size();std::ifstream masks_in(argv[3]);int number;masks_in>>number;
 std::vector<std::string> names(number);std::vector<std::vector<int>> masks(number,std::vector<int>(m+1));
 for(int i=0;i<number;++i){masks_in>>names[i];for(int&v:masks[i])masks_in>>v;}if(!masks_in)return 3;
 std::vector<std::vector<int>> binom(m+1,std::vector<int>(m+1));
 for(int i=0;i<=m;++i){binom[i][0]=binom[i][i]=1;for(int j=1;j<i;++j)binom[i][j]=(binom[i-1][j-1]+binom[i-1][j])%p;}
 using Blocks=std::unordered_map<Label,std::vector<uint32_t>,Hash>;Blocks rowblocks,colblocks;
 for(int a=0;a<d;++a)for(int b=0;b<d;++b){Label diff(modes),sum(modes);for(int j=0;j<modes;++j){diff[j]=labels[a][j]-labels[b][j];sum[j]=labels[a][j]+labels[b][j];}
  rowblocks[diff].push_back(a*d+b);colblocks[sum].push_back(a*d+b);}
 std::unordered_map<std::string,std::vector<uint8_t>> cache;
 auto evaluations=[&](int ai,int yi,int bi)->const std::vector<uint8_t>&{
  const auto&a=labels[ai];const auto&y=labels[yi];const auto&b=labels[bi];std::vector<uint8_t> polynomial(m+1),next(m+1);polynomial[0]=1;
  for(int j=0;j<modes;++j){std::fill(next.begin(),next.end(),0);
   for(int k=0;k<=std::min(a[j],y[j]);++k){int q=a[j]-k;if(q>b[j])continue;int c=binom[y[j]][k]*binom[b[j]][q]%p;if(!c)continue;int degree=k*weights[j];
    for(int i=0;i+degree<=m;++i)if(polynomial[i])next[i+degree]=(next[i+degree]+c*polynomial[i])%p;}
   polynomial.swap(next);
  }
  assert(polynomial[0]==(ai==bi));std::string key(reinterpret_cast<char*>(polynomial.data()),polynomial.size());
  auto old=cache.find(key);if(old!=cache.end())return old->second;
  std::vector<uint8_t> value(number);for(int i=0;i<number;++i){int sum=0;for(int j=0;j<=m;++j)sum+=masks[i][j]*polynomial[j];value[i]=sum%p;}
  return cache.emplace(std::move(key),std::move(value)).first->second;
 };
 std::vector<std::vector<uint64_t>> ranks;uint64_t entries=0;
 for(bool row:{true,false}){const auto&blocks=row?rowblocks:colblocks;std::vector<uint64_t> total(number);int processed=0;
  for(const auto& item:blocks){const auto& pairs=item.second;int n=pairs.size();std::vector<std::vector<uint8_t>> matrices(number,std::vector<uint8_t>(size_t(n)*n));
   for(int i=0;i<n;++i)for(int j=0;j<n;++j){int a=pairs[i]/d,y=row?pairs[i]%d:pairs[j]/d,b=row?pairs[j]/d:pairs[j]%d;
    const auto& value=evaluations(a,y,b);for(int k=0;k<number;++k)matrices[k][size_t(i)*n+j]=value[k];++entries;}
   for(int k=0;k<number;++k)total[k]+=matrix_rank(std::move(matrices[k]),n,p);
   if((++processed%1000)==0)std::cerr<<(row?"row":"col")<<" blocks "<<processed<<"/"<<blocks.size()<<"\n";
  }ranks.push_back(std::move(total));
 }
 long double raw=0,projector=0,refined=0;
 std::cout<<std::setprecision(18)<<"{\"prime\":"<<p<<",\"energy\":"<<h<<",\"reduced_energy\":"<<m<<",\"labels\":"<<d
  <<",\"row_blocks\":"<<rowblocks.size()<<",\"column_blocks\":"<<colblocks.size()<<",\"entries\":"<<entries<<",\"distinct_entry_polynomials\":"<<cache.size()<<",\"pieces\":[";
 for(int i=0;i<number;++i){long double cost=sqrtl((long double)ranks[0][i]*ranks[1][i]);if(names[i].rfind("raw_",0)==0)raw+=cost;else if(names[i].rfind("projector_",0)==0)projector+=cost;else if(names[i].rfind("refined_",0)==0)refined+=cost;
  std::cout<<(i?",":"")<<"{\"name\":\""<<names[i]<<"\",\"row\":"<<ranks[0][i]<<",\"col\":"<<ranks[1][i]<<"}";}
 std::cout<<"],\"rank_sums\":{\"raw\":"<<raw<<",\"projectors\":"<<projector<<",\"refined\":"<<refined<<"},\"threshold\":"<<powl(d,sqrtl(2))
  <<",\"seconds\":"<<std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count()<<"}\n";
 return 0;
}
