#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <vector>
using namespace std;
using U=uint64_t;
static uint32_t n,states,rows,projects,precision;
static U mask;
static vector<U> projection_data;
static uint64_t total_spaces=0,total_chars=0,total_survivors=0;
static void walsh(vector<U>& a) {
 for(size_t gap=1;gap<a.size();gap*=2)
  for(size_t start=0;start<a.size();start+=2*gap)
   for(size_t j=0;j<gap;j++) {U x=a[start+j],y=a[start+j+gap];a[start+j]=(x+y)&mask;a[start+j+gap]=(x-y)&mask;}
}
static void inspect(const vector<int>& basis,const vector<int>& support) {
 total_spaces++;
 const size_t width=support.size();total_chars+=width;
 vector<U> vals(width);
 vector<size_t> candidates;
 for(size_t i=0;i<width;i++)vals[i]=projection_data[support[i]];
 walsh(vals);
 for(size_t i=0;i<width;i++)if(!vals[i])candidates.push_back(i);
 for(uint32_t p=1;p<projects&&!candidates.empty();p++) {
  for(size_t i=0;i<width;i++)vals[i]=projection_data[p*states+support[i]];
  walsh(vals);
  candidates.erase(remove_if(candidates.begin(),candidates.end(),[&](size_t a){return vals[a]!=0;}),candidates.end());
 }
 for(auto alpha:candidates) {
  total_survivors++;
  cout<<"{\"basis\":[";
  for(size_t i=0;i<basis.size();i++){if(i)cout<<',';cout<<basis[i];}
  cout<<"],\"character\":"<<alpha<<"}"<<'\n';
 }
}
static void enumerate_pivots(int rank,int at,int next,vector<int>& pivots) {
 if(at<rank) {for(int p=next;p<=int(n)-(rank-at);p++){pivots[at]=p;enumerate_pivots(rank,at+1,p+1,pivots);}return;}
 vector<pair<int,int>> free;
 for(int i=0;i<rank;i++)for(int j=pivots[i]+1;j<int(n);j++)
  if(find(pivots.begin(),pivots.end(),j)==pivots.end())free.push_back({i,j});
 for(U flags=0;flags<(U(1)<<free.size());flags++) {
  vector<int> basis(rank),support(1,0);
  for(int i=0;i<rank;i++)basis[i]=1<<pivots[i];
  for(size_t i=0;i<free.size();i++)if((flags>>i)&1)basis[free[i].first]|=1<<free[i].second;
  for(auto v:basis){int size=support.size();for(int i=0;i<size;i++)support.push_back(support[i]^v);}
  inspect(basis,support);
 }
}
int main(int argc,char**argv) {
 if(argc!=4)return 2;
 n=stoi(argv[1]);precision=stoi(argv[2]);projects=3;states=1<<n;
 if(n>8||precision>63||precision<1)return 3;
 mask=(U(1)<<precision)-1;
 projection_data.resize(projects*states);
 ifstream input(argv[3],ios::binary);input.read((char*)projection_data.data(),projection_data.size()*sizeof(U));
 if(input.gcount()!=int64_t(projection_data.size()*sizeof(U)))return 4;
 for(int rank=0;rank<=int(n);rank++){vector<int>pivots(rank);enumerate_pivots(rank,0,0,pivots);}
 cerr<<"{\"spaces\":"<<total_spaces<<",\"characters\":"<<total_chars<<",\"survivors\":"<<total_survivors<<"}"<<'\n';
}
