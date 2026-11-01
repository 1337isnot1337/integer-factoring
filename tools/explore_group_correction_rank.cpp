// Exact native screen of central-shear crossing ranks over F_2.
// See research/noncommutative-crossing-orbits.md for the compressed matrix.
// This screens one necessary rank bound, not a complete factoring seed.
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

struct Screen {
  int k, q, n, stride;
  bool lookup;
  std::vector<int> mult;
  std::vector<std::pair<int,int>> terms;
  std::vector<uint64_t> rows;
  std::vector<int> pivot;

  Screen(int degree, const std::string& family): k(degree), q(1<<k), n(q*q*q),
      stride((n+63)/64), lookup(family=="lookup"), mult(q*q), rows(size_t(n)*stride), pivot(n,-1) {
    int modulus = k==1 ? 3 : k==2 ? 7 : k==3 ? 11 : k==4 ? 19 : 0;
    if (!modulus) throw std::string("k must lie between 1 and 4");
    for (int a=0;a<q;++a) for(int b=0;b<q;++b) {
      int x=a,y=b,value=0;
      while(y) { if(y&1)value^=x; y>>=1; x<<=1; if(x&q)x^=modulus; }
      mult[a*q+b]=value;
    }
    if (family=="bilinear") {
      for(int i=0;i<k;++i)for(int j=0;j<k;++j)terms.push_back({i,k+j});
    } else if(family=="quadratic") {
      terms.push_back({-1,-1});
      for(int i=0;i<2*k;++i)terms.push_back({i,-1});
      for(int i=0;i<2*k;++i)for(int j=i+1;j<2*k;++j)terms.push_back({i,j});
    } else if(lookup) {
      for(int a=0;a<q;++a)for(int b=0;b<q;++b)terms.push_back({a,b});
    } else throw std::string("family must be bilinear, quadratic, or lookup");
  }

  std::vector<int> correction(const std::vector<int>& coefficients) const {
    if(lookup)return coefficients;
    std::vector<int> f(q*q);
    for(int a=0;a<q;++a)for(int b=0;b<q;++b) {
      int bits=a|(b<<k),value=0;
      for(size_t j=0;j<terms.size();++j) {
        auto [u,v]=terms[j];
        if((u<0 || ((bits>>u)&1)) && (v<0 || ((bits>>v)&1)))value^=coefficients[j];
      }
      f[a*q+b]=value;
    }
    return f;
  }

  int rank(const std::vector<int>& f, int stop_at=0) {
    std::fill(rows.begin(),rows.end(),uint64_t(0));
    std::fill(pivot.begin(),pivot.end(),-1);
    for(int s=0;s<q;++s)for(int t=0;t<q;++t) {
      int shift=f[s*q+t], st=mult[s*q+t];
      for(int a=0;a<q;++a)for(int b=0;b<q;++b) {
        int delta=mult[a*q+t]^mult[b*q+s]^shift;
        int epsilon=delta^st;
        int row=(a*q+b)*q+delta;
        int col=((a^s)*q+(b^t))*q+epsilon;
        rows[size_t(row)*stride+(col>>6)] |= uint64_t(1)<<(col&63);
      }
    }
    int answer=0;
    for(int i=0;i<n;++i) {
      uint64_t* row=rows.data()+size_t(i)*stride;
      int block=stride-1;
      while(block>=0) {
        while(block>=0 && row[block]==0)--block;
        if(block<0)break;
        int col=block*64+63-__builtin_clzll(row[block]);
        if(pivot[col]<0) {
          pivot[col]=i;
          ++answer;
          if(stop_at && answer>=stop_at)return answer;
          break;
        }
        const uint64_t* other=rows.data()+size_t(pivot[col])*stride;
        for(int j=0;j<=block;++j)row[j]^=other[j];
      }
    }
    return answer;
  }
};

uint64_t random_word(uint64_t& state) {
  state ^= state >> 12; state ^= state << 25; state ^= state >> 27;
  return state*uint64_t(2685821657736338717ULL);
}

int main(int argc,char**argv) {
  try {
    int k=2;std::string family="quadratic";uint64_t start=0,count=1,seed=20261030;
    bool random=false,exact=false,normalized=false;
    for(int i=1;i<argc;++i) {
      std::string arg=argv[i];
      if(arg=="--k" && i+1<argc)k=std::stoi(argv[++i]);
      else if(arg=="--family" && i+1<argc)family=argv[++i];
      else if(arg=="--start" && i+1<argc)start=std::stoull(argv[++i]);
      else if(arg=="--count" && i+1<argc)count=std::stoull(argv[++i]);
      else if(arg=="--seed" && i+1<argc)seed=std::stoull(argv[++i]);
      else if(arg=="--random")random=true;
      else if(arg=="--exact")exact=true;
      else if(arg=="--normalized")normalized=true;
      else throw std::string("unrecognized argument: ")+arg;
    }
    Screen screen(k,family);
    int raw_bits=int(screen.terms.size())*k;
    std::vector<std::pair<int,int>> free_bits;
    for(size_t t=0;t<screen.terms.size();++t)for(int b=0;b<k;++b) {
      auto [u,v]=screen.terms[t];
      bool excluded=false;
      if(normalized && screen.lookup) {
        excluded=(u==0 && v==0) || (u==1 && v==0) || (u==0 && v==1);
        excluded=excluded || (u==1 && v==1 && b==0);
      } else if(normalized) {
        excluded=family=="quadratic" &&
            ((u==-1 && v==-1) || (v==-1 && (u==0 || u==k)));
        excluded=excluded || (u==0 && v==k && b==0);
      }
      if(!excluded)free_bits.push_back({int(t),b});
    }
    int bits=int(free_bits.size());
    if(!random && bits>63)throw std::string("exhaustive index needs at most 63 bits; use --random");
    if(!random && (start>=(uint64_t(1)<<bits) || count>(uint64_t(1)<<bits)-start))
      throw std::string("enumeration exceeds the specified family");
    std::vector<int> coefficients(screen.terms.size(),0),best_coefficients=coefficients;
    int best=screen.rank(screen.correction(coefficients));
    std::vector<uint64_t> histogram(screen.n+1,0);
    uint64_t digest=14695981039346656037ULL;
    auto began=std::chrono::steady_clock::now();
    std::cout<<"{\"event\":\"start\",\"k\":"<<k<<",\"q\":"<<screen.q
             <<",\"matrix_order\":"<<screen.n<<",\"family\":\""<<family
             <<"\",\"coefficient_bits\":"<<raw_bits<<",\"free_bits\":"<<bits
             <<",\"normalized\":"<<(normalized?"true":"false")<<",\"baseline_rank\":"<<best
             <<",\"random\":"<<(random?"true":"false")<<",\"exact_ranks\":"<<(exact?"true":"false")
             <<",\"start\":"<<start<<",\"count\":"<<count<<",\"seed\":"<<seed<<"}"<<std::endl;
    uint64_t random_state=seed?seed:1;
    for(uint64_t j=0;j<count;++j) {
      uint64_t code=start+j;
      std::fill(coefficients.begin(),coefficients.end(),0);
      uint64_t random_bits=random?random_word(random_state):0;
      for(size_t bit=0;bit<free_bits.size();++bit) {
        if(random && bit && bit%64==0)random_bits=random_word(random_state);
        int value=random ? ((random_bits>>(bit%64))&1) : ((code>>bit)&1);
        auto [term,place]=free_bits[bit];
        coefficients[term]|=value<<place;
      }
      int r=screen.rank(screen.correction(coefficients),exact?0:best);
      ++histogram[r];
      digest^=j+start;digest*=1099511628211ULL;digest^=uint64_t(r);digest*=1099511628211ULL;
      if(r<best) {
        best=r;best_coefficients=coefficients;
        std::cout<<"{\"event\":\"improvement\",\"index\":"<<start+j<<",\"rank\":"<<best<<",\"coefficients\":[";
        for(size_t t=0;t<coefficients.size();++t)std::cout<<(t?",":"")<<coefficients[t];
        std::cout<<"]}"<<std::endl;
      }
      if((j+1)%100000==0) {
        auto ms=std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now()-began).count();
        std::cout<<"{\"event\":\"progress\",\"completed\":"<<j+1<<",\"best_rank\":"<<best<<",\"milliseconds\":"<<ms<<"}"<<std::endl;
      }
    }
    auto ms=std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now()-began).count();
    std::cout<<"{\"event\":\"complete\",\"tested\":"<<count<<",\"best_rank\":"<<best<<",\"digest_fnv64\":\""<<digest<<"\",\"milliseconds\":"<<ms<<",\"histogram\":[";
    bool comma=false;for(int r=0;r<=screen.n;++r)if(histogram[r]){
      std::cout<<(comma?",":"")<<"["<<r<<","<<histogram[r]<<"]";comma=true;
    }
    std::cout<<"],\"best_coefficients\":[";
    for(size_t t=0;t<best_coefficients.size();++t)std::cout<<(t?",":"")<<best_coefficients[t];
    std::cout<<"]}"<<std::endl;
    return 0;
  } catch(const std::string&error) {std::cerr<<error<<std::endl;return 2;}
}
