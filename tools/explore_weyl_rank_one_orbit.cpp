#include <algorithm>
#include <cassert>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using Vector=std::vector<uint64_t>;
struct Edge {int row,col; uint64_t polynomial;};
uint32_t read32(std::ifstream& f){uint32_t v;f.read(reinterpret_cast<char*>(&v),4);return v;}
uint64_t read64(std::ifstream& f){uint64_t v;f.read(reinterpret_cast<char*>(&v),8);return v;}

struct Data {
 int h=0,d=0,n=0,words=0;
 std::vector<int> global,local;
 std::vector<std::vector<Edge>> edges;
 Vector diagonal;
};

void read_dump(Data& data,const std::string& path,bool first){
 std::ifstream f(path,std::ios::binary);if(!f)throw std::runtime_error("cannot open dump");
 f.exceptions(std::ios::badbit|std::ios::failbit);
 char magic[8];f.read(magic,8);assert(std::string(magic,7)=="WBIN001");
 int h=read32(f),d=read32(f),blocks=read32(f),polywords=read32(f),nw=read32(f),rowflag=read32(f);
 assert(polywords==1&&h<64);std::vector<int> weights(nw);for(auto&x:weights)x=read32(f);
 std::vector<std::vector<int>> labels(d,std::vector<int>(nw));for(auto& a:labels)for(auto&x:a)x=read32(f);
 if(first){data.h=h;data.d=d;data.local.assign(d*d,-1);
  for(int a=0;a<d;++a)for(int y=0;y<d;++y){bool even=true;for(int j=0;j<nw;++j)if((labels[a][j]+labels[y][j])&1)even=false;
   if(even){data.local[a*d+y]=data.global.size();data.global.push_back(a*d+y);}}
  data.n=data.global.size();data.words=(data.n+63)/64;data.diagonal.assign(data.words,0);
  for(int a=0;a<d;++a){int i=data.local[a*d+a];data.diagonal[i/64]|=uint64_t(1)<<(i%64);}
 }else assert(h==data.h&&d==data.d);
 std::vector<Edge> edges;
 for(int block=0;block<blocks;++block){int n=read32(f);std::vector<int> pairs(n);for(auto&x:pairs)x=read32(f);
  bool retained=data.local[pairs[0]]>=0;
  for(int i=0;i<n;++i)for(int j=0;j<n;++j){uint64_t p=read64(f);if(!retained||!p)continue;
   int r=data.local[pairs[i]],code=pairs[j];if(!rowflag)code=(code%d)*d+code/d;
   int c=data.local[code];assert(r>=0&&c>=0);edges.push_back({r,c,p});}
 }
 data.edges.push_back(std::move(edges));
}

struct Orbit {std::vector<Vector> basis;std::vector<int> depth;int maxdepth=0;};
Orbit orbit(const Data& data,bool transpose){
 Orbit out;std::vector<int> pivots(data.n,-1);
 auto insert=[&](Vector value,int depth){for(int w=data.words-1;w>=0;--w)while(value[w]){
  int p=64*w+63-__builtin_clzll(value[w]);if(pivots[p]<0){pivots[p]=out.basis.size();out.basis.push_back(std::move(value));out.depth.push_back(depth);out.maxdepth=std::max(out.maxdepth,depth);return;}
  const auto& b=out.basis[pivots[p]];for(int k=0;k<=w;++k)value[k]^=b[k];
 }};
 insert(data.diagonal,0);
 for(size_t pos=0;pos<out.basis.size();++pos){Vector input=out.basis[pos];int nextdepth=out.depth[pos]+1;
  for(const auto& edges:data.edges){std::vector<uint64_t> polys(data.n,0);
   for(const auto&e:edges){int r=transpose?e.col:e.row,c=transpose?e.row:e.col;if((input[c/64]>>(c%64))&1)polys[r]^=e.polynomial;}
   std::vector<Vector> images(data.h+1,Vector(data.words));
   for(int r=0;r<data.n;++r){uint64_t p=polys[r];while(p){int k=__builtin_ctzll(p);images[k][r/64]|=uint64_t(1)<<(r%64);p&=p-1;}}
   for(auto& v:images)insert(std::move(v),nextdepth);
  }
  if((pos+1)%500==0)std::cerr<<(transpose?"right":"left")<<" processed "<<pos+1<<" basis "<<out.basis.size()<<"\n";
 }
 return out;
}

int pairing_rank(const Orbit& left,const Orbit& right){
 int n=left.basis.size(),words=(n+63)/64;std::vector<Vector> pivots(n);int answer=0;
 for(const auto& b:right.basis){Vector row(words,0);for(int i=0;i<n;++i){int parity=0;for(size_t w=0;w<b.size();++w)parity^=__builtin_parityll(b[w]&left.basis[i][w]);if(parity)row[i/64]|=uint64_t(1)<<(i%64);}
  for(int w=words-1;w>=0;--w)while(row[w]){int p=64*w+63-__builtin_clzll(row[w]);if(pivots[p].empty()){pivots[p]=std::move(row);++answer;w=-1;break;}for(int j=0;j<=w;++j)row[j]^=pivots[p][j];}
 }
 return answer;
}

int main(int argc,char** argv){if(argc!=2)return 2;auto started=std::chrono::steady_clock::now();Data data;
 std::string prefix=argv[1];read_dump(data,prefix+".row.bin",true);read_dump(data,prefix+".col.bin",false);
 auto left=orbit(data,false),right=orbit(data,true);int gram=pairing_rank(left,right);
 auto symmetry_checks=[&](const Orbit& o){int nonsymmetric=0,nonconstant_diagonal=0;
  for(const auto&v:o.basis){bool badsym=false,baddiag=false;int diagonal=-1;
   for(int i=0;i<data.n;++i){int code=data.global[i],a=code/data.d,y=code%data.d;
    int j=data.local[y*data.d+a];bool x=(v[i/64]>>(i%64))&1,z=(v[j/64]>>(j%64))&1;
    if(x!=z)badsym=true;if(a==y){if(diagonal<0)diagonal=x;else if(x!=diagonal)baddiag=true;}}
   nonsymmetric+=badsym;nonconstant_diagonal+=baddiag;
  }return std::pair<int,int>(nonsymmetric,nonconstant_diagonal);};
 auto lc=symmetry_checks(left),rc=symmetry_checks(right);
 std::cout<<"{\"energy\":"<<data.h<<",\"labels\":"<<data.d<<",\"parity_zero_dimension\":"<<data.n
  <<",\"nonzero_polynomial_entries\":["<<data.edges[0].size()<<","<<data.edges[1].size()<<"]"
  <<",\"left_orbit\":"<<left.basis.size()<<",\"right_orbit\":"<<right.basis.size()<<",\"pairing_rank\":"<<gram
  <<",\"maximum_word_lengths\":["<<left.maxdepth<<","<<right.maxdepth<<"]"
  <<",\"nonsymmetric_orbit_vectors\":["<<lc.first<<","<<rc.first<<"]"
  <<",\"nonconstant_diagonal_orbit_vectors\":["<<lc.second<<","<<rc.second<<"]"
  <<",\"seconds\":"<<std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count()<<"}\n";
 return 0;
}
