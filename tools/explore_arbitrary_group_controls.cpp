#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

// Complete coefficient-space searches over arbitrary D8 control permutations.
// The group coordinates are (a,b,c), with central cocycle a*B.
using Matrix = std::array<uint64_t,64>;

int mul(int x,int y) { return x^y^(((x&1)*((y>>1)&1))<<2); }
int inv(int x) { return x^(((x&1)*((x>>1)&1))<<2); }
int rank_rows(const uint64_t* rows,int n) {
  uint64_t pivots[64]{}; int rank=0;
  for(int i=0;i<n;++i) {
    uint64_t row=rows[i];
    while(row) {
      int lead=63-__builtin_clzll(row);
      if(pivots[lead]) row^=pivots[lead];
      else {pivots[lead]=row; ++rank; break;}
    }
  }
  return rank;
}

double shortest(const std::vector<double>& costs,int target,std::vector<int>* witness=nullptr) {
  int n=costs.size();
  std::vector<double> dist(n,std::numeric_limits<double>::infinity());
  std::vector<int> parent(n,-1),edge(n,0); std::vector<bool> done(n,false);
  dist[0]=0;
  for(int step=0;step<n;++step) {
    int u=-1;
    for(int i=0;i<n;++i) if(!done[i]&&(u<0||dist[i]<dist[u]))u=i;
    if(u==target)break;
    done[u]=true;
    for(int v=1;v<n;++v) {
      int next=u^v; double value=dist[u]+costs[v];
      if(value<dist[next]) {dist[next]=value; parent[next]=u;edge[next]=v;}
    }
  }
  if(witness) {
    witness->clear();
    for(int at=target;at;at=parent[at]) {
      if(parent[at]<0)throw std::runtime_error("unreachable endpoint");
      witness->push_back(edge[at]);
    }
  }
  return dist[target];
}

void print_array(const std::vector<int>& values) {
  std::cout<<'[';
  for(size_t i=0;i<values.size();++i)std::cout<<(i?",":"")<<values[i];
  std::cout<<']';
}

int main(int argc,char**argv) {
  try {
    std::array<int,3> weights{1,1,1};
    if(argc!=1&&argc!=4)throw std::runtime_error("usage: arbitrary-group-controls [w0 w1 w2]");
    if(argc==4)for(int j=0;j<3;++j)weights[j]=std::stoi(argv[j+1]);
    int degree=std::accumulate(weights.begin(),weights.end(),0);
    if(*std::min_element(weights.begin(),weights.end())<1||degree>15)
      throw std::runtime_error("positive weights with sum at most15 required");
    auto start=std::chrono::steady_clock::now();
    std::array<std::array<uint64_t,8>,8> kernels{};
    for(int bits=0;bits<8;++bits) {
      std::array<uint64_t,8> values{}; values[0]=1;
      for(int j=0;j<3;++j)if((bits>>j)&1) {
        std::array<uint64_t,8> next{};
        for(int a=0;a<8;++a) {
          next[a]^=values[a]<<weights[j];
          next[mul(a,1<<j)]^=values[a]^(values[a]<<weights[j]);
        }
        values=next;
      }
      kernels[bits]=values;
      for(int a=0;a<8;++a)if((__builtin_popcountll(values[a])&1)!=(a==0))
        throw std::runtime_error("kernel endpoint failed");
    }
    std::vector<uint64_t> coefficient(degree+1,0);
    for(int bits=0;bits<8;++bits)for(int a=0;a<8;++a)for(int j=0;j<=degree;++j)
      if((kernels[bits][a]>>j)&1)coefficient[j]^=uint64_t(1)<<(8*bits+a);
    uint64_t pivots[64]{}; int coordinates[64]{};
    std::vector<int> indices; std::vector<uint64_t> basis;
    for(int j=0;j<=degree;++j) {
      uint64_t value=coefficient[j]; int label=0;
      while(value) {
        int p=63-__builtin_clzll(value);
        if(pivots[p]) {value^=pivots[p];label^=coordinates[p];}
        else {
          label^=1<<indices.size();pivots[p]=value;coordinates[p]=label;
          indices.push_back(j);basis.push_back(coefficient[j]);break;
        }
      }
    }
    int dimension=indices.size(),count=1<<dimension;
    if(dimension>8)throw std::runtime_error("unexpected coefficient dimension");
    uint64_t endpoint=0;for(auto value:coefficient)endpoint^=value;
    int target=0;
    while(endpoint) {
      int p=63-__builtin_clzll(endpoint);
      if(!pivots[p])throw std::runtime_error("endpoint outside basis");
      endpoint^=pivots[p];target^=coordinates[p];
    }
    std::vector<uint64_t> combined(count,0);
    std::vector<int> rowrank(count,0);
    for(int mask=1;mask<count;++mask) {
      int bit=__builtin_ctz(mask);combined[mask]=combined[mask^(1<<bit)]^basis[bit];
      for(int bits=0;bits<8;++bits) {
        uint64_t rows[8]{};int support=(combined[mask]>>(8*bits))&255;
        for(int y=0;y<8;++y)for(int a=0;a<8;++a)if((support>>a)&1)
          rows[y]^=uint64_t(1)<<mul(y,a);
        rowrank[mask]+=rank_rows(rows,8);
      }
    }
    if(rowrank[target]!=64)throw std::runtime_error("row endpoint rank failed");
    std::cout.precision(17);
    std::cout<<"{\"event\":\"start\",\"labels\":8,\"weights\":["<<weights[0]<<','<<weights[1]<<','<<weights[2]
             <<"],\"dimension\":"<<dimension<<",\"basis_indices\":";
    print_array(indices);std::cout<<",\"target\":"<<target<<",\"row_ranks\":";print_array(rowrank);std::cout<<"}\n";
    std::array<int,8> permutation{0,1,2,3,4,5,6,7},best_permutation{};
    double best=std::numeric_limits<double>::infinity(); int tested=0,min_integer=1000000;
    std::vector<int> best_columns,best_witness;
    uint64_t digest=1469598103934665603ULL;
    std::vector<int> histogram(1000,0);
    do {
      std::vector<Matrix> matrices(dimension);
      for(int S=0;S<8;++S)for(int y=0;y<8;++y) {
        int x=mul(S,y);
        for(int a=0;a<8;++a) {
          uint64_t polynomial=kernels[permutation[S]][a];
          if(!polynomial)continue;
          int u=mul(y,a),v=mul(x,a);
          for(int b=0;b<dimension;++b)if((polynomial>>indices[b])&1)
            matrices[b][8*x+u]^=uint64_t(1)<<(8*y+v);
        }
      }
      Matrix work{};int previous=0;std::vector<int> columns(count,0);
      std::vector<double> costs(count,0),lower(count,0);
      for(int k=1;k<count;++k) {
        int mask=k^(k>>1),changed=__builtin_ctz(mask^previous);
        for(int i=0;i<64;++i)work[i]^=matrices[changed][i];
        columns[mask]=rank_rows(work.data(),64);
        int product=rowrank[mask]*columns[mask];
        costs[mask]=std::sqrt(double(product));
        int floor=int(costs[mask]);while((floor+1)*(floor+1)<=product)++floor;while(floor*floor>product)--floor;
        lower[mask]=floor;
        digest^=uint64_t(mask*65536+rowrank[mask]*256+columns[mask]);digest*=1099511628211ULL;
        previous=mask;
      }
      if(columns[target]!=64)throw std::runtime_error("column endpoint rank failed");
      double value=shortest(costs,target);int integer=int(shortest(lower,target));
      min_integer=std::min(min_integer,integer);
      if(integer>=int(histogram.size()))throw std::runtime_error("histogram overflow");
      ++histogram[integer];
      if(value<best-1e-10) {
        best=value;best_permutation=permutation;best_columns=columns;
        shortest(costs,target,&best_witness);
        std::cout<<"{\"event\":\"improvement\",\"index\":"<<tested<<",\"cost\":"<<best
                 <<",\"integer_lower\":"<<integer<<",\"control_word_permutation\":";
        print_array(std::vector<int>(permutation.begin(),permutation.end()));std::cout<<",\"pieces\":";
        print_array(best_witness);std::cout<<"}\n"<<std::flush;
      }
      ++tested;
    } while(std::next_permutation(permutation.begin(),permutation.end()));
    auto ms=std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now()-start).count();
    std::cout<<"{\"event\":\"complete\",\"permutations_tested\":"<<tested<<",\"best_cost\":"<<best
             <<",\"global_integer_lower\":"<<min_integer<<",\"best_control_word_permutation\":";
    print_array(std::vector<int>(best_permutation.begin(),best_permutation.end()));
    std::cout<<",\"best_column_ranks\":";print_array(best_columns);
    std::cout<<",\"best_pieces\":";print_array(best_witness);
    std::cout<<",\"integer_lower_histogram\":[";
    bool first=true;for(int i=0;i<int(histogram.size());++i)if(histogram[i]) {
      std::cout<<(first?"":",")<<'['<<i<<','<<histogram[i]<<']';first=false;
    }
    std::cout<<"],\"digest_fnv64\":\""<<digest<<"\",\"milliseconds\":"<<ms<<"}\n";
  } catch(const std::exception&e) {std::cerr<<e.what()<<'\n';return 1;}
}
