// Exact projective F4 screen of supplied binary tensor coefficient spaces.
// F4 = F2[t]/(t^2+t+1); every flattening uses all supplied direct-sum blocks.
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

using Matrix=std::vector<uint64_t>;
using Blocks=std::vector<Matrix>;
using Basis=std::vector<Blocks>;

int binary_rank_capped(Matrix rows,int n,int stride,int cap) {
  std::vector<int>pivot(n,-1);int rank=0;
  for(int i=0;i<n;++i){size_t offset=size_t(i)*stride;int block=stride-1;
    while(block>=0){while(block>=0&&!rows[offset+block])--block;if(block<0)break;
      int column=64*block+63-__builtin_clzll(rows[offset+block]);
      if(pivot[column]<0){pivot[column]=i;if(++rank>cap)return rank;break;}
      size_t previous=size_t(pivot[column])*stride;
      for(int j=0;j<=block;++j)rows[offset+j]^=rows[previous+j];
    }
  }
  return rank;
}

int field_rank(Matrix lo,Matrix hi,int n,int stride,int cap) {
  std::vector<int>pivot(n,-1);int rank=0;
  for(int i=0;i<n;++i) {
    size_t offset=size_t(i)*stride;int block=stride-1;
    while(block>=0) {
      while(block>=0&&!(lo[offset+block]|hi[offset+block]))--block;
      if(block<0)break;
      uint64_t both=lo[offset+block]|hi[offset+block];
      int bit=63-__builtin_clzll(both),column=64*block+bit;
      int scalar=int((lo[offset+block]>>bit)&1)+2*int((hi[offset+block]>>bit)&1);
      if(pivot[column]<0) {
        // Normalize the leading entry by scalar^-1; 2^-1=3 and 3^-1=2.
        if(scalar!=1)for(int j=0;j<=block;++j) {
          uint64_t a=lo[offset+j],b=hi[offset+j];
          if(scalar==2){lo[offset+j]=a^b;hi[offset+j]=a;}
          else{lo[offset+j]=b;hi[offset+j]=a^b;}
        }
        pivot[column]=i;if(++rank>cap)return rank;break;
      }
      size_t previous=size_t(pivot[column])*stride;
      for(int j=0;j<=block;++j) {
        uint64_t a=lo[previous+j],b=hi[previous+j];
        if(scalar==1){lo[offset+j]^=a;hi[offset+j]^=b;}
        else if(scalar==2){lo[offset+j]^=b;hi[offset+j]^=a^b;}
        else{lo[offset+j]^=a^b;hi[offset+j]^=a;}
      }
    }
  }
  return rank;
}

int combination_rank(const Basis&basis,const std::vector<int>&sizes,int block,
                     uint64_t code,int cap,bool binary_only) {
  int n=sizes[block],stride=(n+63)/64;size_t words=size_t(n)*stride;
  Matrix lo(words,0),hi(words,0);
  for(size_t i=0;i<basis.size();++i) {
    int scalar=(code>>(2*i))&3;
    if(scalar&1)for(size_t j=0;j<words;++j)lo[j]^=basis[i][block][j];
    if(scalar&2)for(size_t j=0;j<words;++j)hi[j]^=basis[i][block][j];
  }
  if(binary_only)return binary_rank_capped(std::move(lo),n,stride,cap);
  return field_rank(std::move(lo),std::move(hi),n,stride,cap);
}

int main(int argc,char**argv) {
  try {
    uint64_t first=1,count=0;bool binary_only=false;
    for(int i=1;i<argc;++i){std::string arg=argv[i];
      if(arg=="--start"&&i+1<argc)first=std::stoull(argv[++i]);
      else if(arg=="--count"&&i+1<argc)count=std::stoull(argv[++i]);
      else if(arg=="--field2")binary_only=true;
      else throw std::runtime_error("unknown argument: "+arg);
    }
    int labels,dimension,blocks,gate;
    if(!(std::cin>>labels>>dimension>>blocks>>gate))throw std::runtime_error("missing header");
    if(labels<2||dimension<1||dimension>(binary_only?20:12)||blocks<1||gate<2)throw std::runtime_error("invalid dimensions");
    std::vector<int>sizes(blocks);uint64_t total_rows=0;
    for(int&n:sizes){std::cin>>n;if(n<=0)throw std::runtime_error("invalid block size");total_rows+=n;}
    if(total_rows!=uint64_t(labels)*labels)throw std::runtime_error("blocks do not cover a complete cut");
    Basis cuts[2];
    for(Basis&basis:cuts){basis.resize(dimension);
      for(Blocks&matrices:basis){matrices.resize(blocks);
        for(int b=0;b<blocks;++b){int n=sizes[b],stride=(n+63)/64;matrices[b].resize(size_t(n)*stride);
          for(uint64_t&word:matrices[b]){std::string text;if(!(std::cin>>text))throw std::runtime_error("truncated matrix input");word=std::stoull(text,nullptr,16);}
        }
      }
    }
    uint64_t limit=uint64_t(1)<<((binary_only?1:2)*dimension);
    if(!count)count=limit-first;
    if(!first||first>=limit||count>limit-first)throw std::runtime_error("invalid coefficient interval");
    uint64_t product_cap=uint64_t(gate)*gate-1,tested=0,eligible=0;
    auto began=std::chrono::steady_clock::now();
    std::cout<<"{\"event\":\"start\",\"labels\":"<<labels<<",\"dimension\":"<<dimension
             <<",\"field_size\":"<<(binary_only?2:4)<<",\"blocks\":"<<blocks<<",\"integer_cost_gate\":"<<gate<<",\"product_cap\":"<<product_cap
             <<",\"code_start\":"<<first<<",\"code_count\":"<<count<<"}"<<std::endl;
    for(uint64_t index=first;index<first+count;++index) {
      uint64_t code=index;
      if(binary_only){code=0;for(int i=0;i<dimension;++i)code|=((index>>i)&1)<<(2*i);}
      int last=(63-__builtin_clzll(code))/2;
      if(((code>>(2*last))&3)!=1)continue;
      ++tested;uint64_t ranks[2]={0,0};bool rejected=false;
      for(int b=0;b<blocks&&!rejected;++b)for(int cut=0;cut<2;++cut) {
        uint64_t other=ranks[1-cut],maximum=other?product_cap/other:total_rows;
        if(ranks[cut]>maximum){rejected=true;break;}
        int cap=int(std::min<uint64_t>(sizes[b],maximum-ranks[cut]));
        ranks[cut]+=combination_rank(cuts[cut],sizes,b,code,cap,binary_only);
        if(ranks[0]*ranks[1]>product_cap){rejected=true;break;}
      }
      if(!rejected){++eligible;
        std::cout<<"{\"event\":\"eligible\",\"code\":"<<code<<",\"row_rank\":"<<ranks[0]
                 <<",\"column_rank\":"<<ranks[1]<<",\"product\":"<<ranks[0]*ranks[1]<<"}"<<std::endl;
      }
      if(tested%10000==0)std::cout<<"{\"event\":\"progress\",\"projective_tested\":"<<tested<<",\"eligible\":"<<eligible<<",\"code\":"<<code<<"}"<<std::endl;
    }
    auto ms=std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now()-began).count();
    std::cout<<"{\"event\":\"complete\",\"projective_tested\":"<<tested<<",\"eligible\":"<<eligible<<",\"milliseconds\":"<<ms<<"}"<<std::endl;
    return 0;
  }catch(const std::exception&e){std::cerr<<e.what()<<std::endl;return 2;}
}
