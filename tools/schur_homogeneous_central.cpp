// Stream exact homogeneous Schur central ranks from channel-space bases.
// This is a finite computation, not a stronger factoring theorem.
#include "schur_fermion_core.hpp"
#include <chrono>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <sstream>
#include <string>

using namespace schur_fermion;
using Bits=std::vector<uint64_t>;

class Basis {
 public:
  int width, words, rank=0;
  std::vector<Bits> rows;
  explicit Basis(int w):width(w),words((w+63)/64),rows(w) {}
  bool insert(Bits value) {
    for(int word=words-1;word>=0;--word) {
      while(value[word]) {
        const int p=word*64+63-__builtin_clzll(value[word]);
        if(rows[p].empty()) {rows[p]=std::move(value);++rank;return true;}
        for(int j=0;j<=word;++j) value[j]^=rows[p][j];
      }
    }
    return false;
  }
  std::vector<const Bits*> vectors() const {
    std::vector<const Bits*> answer;
    for(const auto& row:rows) if(!row.empty()) answer.push_back(&row);
    return answer;
  }
};

int pairing_rank_dense(const Basis& a,const Basis& b) {
  auto av=a.vectors(), bv=b.vectors();
  Basis gram(b.rank);
  for(const Bits* row:av) {
    Bits out(gram.words);
    for(int j=0;j<(int)bv.size();++j) {
      unsigned parity=0;
      for(int k=0;k<a.words;++k) parity^=__builtin_parityll((*row)[k]&(*bv[j])[k]);
      out[j/64]|=uint64_t(parity)<<(j%64);
    }
    gram.insert(std::move(out));
  }
  return gram.rank;
}

std::vector<Bits> perpendicular(const Basis& basis) {
  std::vector<Bits> answer;
  for(int free=0;free<basis.width;++free) {
    if(!basis.rows[free].empty()) continue;
    Bits v(basis.words);
    v[free/64]|=uint64_t(1)<<(free%64);
    for(int pivot=0;pivot<basis.width;++pivot) {
      const auto& row=basis.rows[pivot];
      if(row.empty())continue;
      unsigned parity=0;
      for(int k=0;k<=pivot/64;++k) parity^=__builtin_parityll(row[k]&v[k]);
      if(parity)v[pivot/64]|=uint64_t(1)<<(pivot%64);
    }
    answer.push_back(std::move(v));
  }
  return answer;
}

Bits remainder(Bits value,const Basis& basis) {
  for(int word=basis.words-1;word>=0;--word) {
    uint64_t remaining=value[word];
    while(remaining) {
      const int p=word*64+63-__builtin_clzll(remaining);
      if(!basis.rows[p].empty()) {
        for(int j=0;j<=word;++j)value[j]^=basis.rows[p][j];
        remaining=value[word];
      } else remaining&=~(uint64_t(1)<<(p%64));
    }
  }
  return value;
}

int pairing_rank_small_codimension(const Basis& a,const Basis& b) {
  Basis quotient(a.width);
  for(auto v:perpendicular(b))quotient.insert(remainder(std::move(v),a));
  return a.rank-(b.width-b.rank)+quotient.rank;
}

int pairing_rank(const Basis& a,const Basis& b) {
  if(b.width-b.rank<=16)return pairing_rank_small_codimension(a,b);
  if(a.width-a.rank<=16)return pairing_rank_small_codimension(b,a);
  return pairing_rank_dense(a,b);
}

void print_supports(const std::vector<Bits>& vectors,int width) {
  std::cout<<'[';
  bool first=true;
  for(const auto& vector:vectors) {
    if(!first)std::cout<<',';
    first=false;
    std::cout<<'[';
    bool inner_first=true;
    for(int i=0;i<width;++i)if((vector[i/64]>>(i%64))&1) {
      if(!inner_first)std::cout<<',';
      inner_first=false;
      std::cout<<i;
    }
    std::cout<<']';
  }
  std::cout<<']';
}

void check_residue_identities(const Bits& values,const Partition& label,
                             const std::vector<int>& hooks,int h) {
  for(int modulus=2;modulus<=h;modulus*=2) {
    std::vector<unsigned> actual(modulus),charge(modulus);
    for(int arm=0;arm<2*h;++arm) {
      const int i=hooks[arm];
      actual[arm%modulus]^=(values[i/64]>>(i%64))&1;
    }
    for(int arm:label.arms)charge[arm%modulus]^=1;
    for(int leg:label.legs)charge[modulus-1-leg%modulus]^=1;
    for(int residue=0;residue<modulus;++residue)
      if(actual[residue]!=(charge[residue]&(values[0]&1)))
        throw std::runtime_error("hook-residue identity failed");
  }
}

Partition parse(const std::string& input) {
  if(input=="-")return Partition();
  std::istringstream stream(input);std::string value;std::vector<int> p;
  while(std::getline(stream,value,','))p.push_back(std::stoi(value));
  return Partition(p);
}

int main(int argc,char** argv) {
  if(argc>1 && std::string(argv[1])=="--queries") {
    std::string p,o,i;
    while(std::cin>>p>>o>>i)std::cout<<mode2(parse(p),parse(o),parse(i))<<'\n';
    return 0;
  }
  const int h=argc>1?std::stoi(argv[1]):8;
  if(h<=0 || h>16 || (h&(h-1)))throw std::runtime_error("this finite scan supports powers of two through 16");
  const auto external=partitions(h), middle=partitions(2*h);
  const int d=external.size(), k=middle.size()+1, total=d*d;
  const int limit=argc>2?std::min(total,std::stoi(argv[2])):total;
  int multiplier=65537;
  while(std::gcd(multiplier,total)!=1)++multiplier;
  const int offset=1729%total;
  Basis left(k),right(k);Partition vacuum;
  std::vector<int> hooks(2*h);
  for(int j=0;j<(int)middle.size();++j)
    if(middle[j].rank()==1)hooks[middle[j].arms[0]]=j+1;
  auto start=std::chrono::steady_clock::now();
  for(int step=0;step<limit;++step) {
    const int pair=(int)((1LL*step*multiplier+offset)%total),x=pair/d,y=pair%d;
    Bits l(left.words),r(right.words);
    l[0]=mode2(external[x],vacuum,external[y]);
    r[0]=mode2(external[x],external[y],vacuum);
    for(int j=0;j<(int)middle.size();++j) {
      if(mode2(external[x],middle[j],external[y]))l[(j+1)/64]|=uint64_t(1)<<((j+1)%64);
      if(mode2(external[x],external[y],middle[j]))r[(j+1)/64]|=uint64_t(1)<<((j+1)%64);
    }
    check_residue_identities(l,external[y],hooks,h);
    check_residue_identities(r,external[y],hooks,h);
    left.insert(std::move(l));right.insert(std::move(r));
    if((step+1)%256==0)std::cerr<<"pairs "<<step+1<<" ranks "<<left.rank<<','<<right.rank<<'\n';
  }
  int rank=pairing_rank(left,right);
  if(h<=8 && rank!=pairing_rank_dense(left,right))
    throw std::runtime_error("independent channel-pairing rank methods disagree");
  double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
  std::cout<<"{\"h\":"<<h<<",\"dimension\":"<<d<<",\"channels\":"<<k
    <<",\"pairs\":"<<limit<<",\"all_pairs\":"<<(limit==total?"true":"false")
    <<",\"pair_multiplier\":"<<multiplier<<",\"pair_offset\":"<<offset
    <<",\"left_rank\":"<<left.rank<<",\"right_rank\":"<<right.rank
    <<",\"hook_residue_identities_checked\":"<<(2LL*(2*h-2)*limit)
    <<",\"central_rank"<<(limit==total?"":"_lower_bound")<<"\":"<<rank
    <<",\"seconds\":"<<seconds;
  if(left.width-left.rank<=8 && right.width-right.rank<=8) {
    std::cout<<",\"left_perpendicular_supports\":";
    print_supports(perpendicular(left),k);
    std::cout<<",\"right_perpendicular_supports\":";
    print_supports(perpendicular(right),k);
  }
  std::cout<<"}\n";
}
