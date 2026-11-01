#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <string>
#include <unordered_map>
#include <vector>

// Exact characteristic-two ranks for a nearest-neighbor grading of
// R/H/C tensor words. Initial screens use t<=6 and one packed word.
struct Candidate {
  std::array<int,3> vertex{};
  std::array<int,9> edge{};
  std::vector<int> row, col;
  bool polynomial_edges=false;
  bool auxiliary_matrices=false;
  int boundary=0;
  uint64_t maximum_product=0;
  long double rank_sum=0;
};

struct Ranker {
  int t, n;
  std::array<uint64_t,6> lower{};
  std::unordered_map<uint64_t,uint8_t> cache;
  uint64_t calls=0, hits=0, gaussian=0;
  explicit Ranker(int t_):t(t_),n(1<<t_) {
    for(int j=0;j<t;j++) for(int h=0;h<n;h++)
      if(!(h&(1<<j))) lower[j]|=uint64_t(1)<<h;
    cache.reserve(1<<18);
  }
  uint64_t translate(uint64_t word,int shift) const {
    for(int j=0;j<t;j++) if(shift&(1<<j)) {
      int delta=1<<j;
      word=((word&lower[j])<<delta)|((word&~lower[j])>>delta);
    }
    return word;
  }
  int direct(uint64_t kernel) const {
    uint64_t basis[64]={};
    int answer=0;
    for(int x=0;x<n;x++) {
      uint64_t word=translate(kernel,x);
      while(word) {
        int pivot=63-__builtin_clzll(word);
        if(!basis[pivot]) {basis[pivot]=word;answer++;break;}
        word^=basis[pivot];
      }
    }
    return answer;
  }
  int rank(uint64_t kernel) {
    calls++;
    if(!kernel) return 0;
    if(__builtin_popcountll(kernel)&1) return n;
    // A nonzero linear term in the square-zero group algebra gives
    // a free variable, hence rank n/2.
    for(int j=0;j<t;j++)
      if(__builtin_popcountll(kernel&~lower[j])&1) return n/2;
    auto found=cache.find(kernel);
    if(found!=cache.end()) {hits++;return found->second;}
    int answer=direct(kernel);
    gaussian++;
    if(cache.size()>(1u<<21)) cache.clear();
    cache.emplace(kernel,uint8_t(answer));
    return answer;
  }
};

Candidate evaluate(int t,const std::array<int,3>& vertex,
                   const std::array<int,9>& edge,Ranker& ranker,
                   bool polynomial_edges=false) {
  Candidate result;
  result.vertex=vertex;result.edge=edge;
  result.polynomial_edges=polynomial_edges;
  int n=1<<t;
  int maximum_edge=*std::max_element(edge.begin(),edge.end());
  if(polynomial_edges) maximum_edge=31-__builtin_clz(maximum_edge);
  int maximum_degree=t*(*std::max_element(vertex.begin(),vertex.end()))+(t-1)*maximum_edge;
  assert(maximum_degree<64);
  std::vector<std::vector<uint64_t>> rows(maximum_degree+1,std::vector<uint64_t>(n));
  std::vector<std::vector<uint64_t>> cols(maximum_degree+1,std::vector<uint64_t>(n));
  for(int g=0;g<n;g++) for(int h=0;h<n;h++) {
    uint64_t dp[3]={};
    dp[0]=(g&1)?0:uint64_t(1)<<vertex[0];
    dp[1]=uint64_t(1)<<vertex[1];
    dp[2]=(h&1)?0:uint64_t(1)<<vertex[2];
    for(int position=1;position<t;position++) {
      uint64_t next[3]={};
      for(int symbol=0;symbol<3;symbol++) {
        if((symbol==0&&(g&(1<<position))) ||
           (symbol==2&&(h&(1<<position)))) continue;
        for(int previous=0;previous<3;previous++) {
          if(polynomial_edges) {
            unsigned polynomial=edge[3*previous+symbol];
            while(polynomial) {
              int degree=__builtin_ctz(polynomial);polynomial&=polynomial-1;
              next[symbol]^=dp[previous]<<(vertex[symbol]+degree);
            }
          } else next[symbol]^=dp[previous]<<(vertex[symbol]+edge[3*previous+symbol]);
        }
      }
      std::copy(next,next+3,dp);
    }
    uint64_t polynomial=dp[0]^dp[1]^dp[2];
    assert((__builtin_popcountll(polynomial)&1)==(g==h));
    while(polynomial) {
      int degree=__builtin_ctzll(polynomial);
      polynomial&=polynomial-1;
      rows[degree][g]|=uint64_t(1)<<h;
      cols[degree][h]|=uint64_t(1)<<g;
    }
  }
  result.row.resize(maximum_degree+1);
  result.col.resize(maximum_degree+1);
  for(int degree=0;degree<=maximum_degree;degree++) {
    for(int g=0;g<n;g++) {
      result.row[degree]+=ranker.rank(rows[degree][g]);
      result.col[degree]+=ranker.rank(cols[degree][g]);
    }
    uint64_t product=uint64_t(result.row[degree])*result.col[degree];
    result.maximum_product=std::max(result.maximum_product,product);
    result.rank_sum+=std::sqrt((long double)product);
    assert((result.row[degree]==0)==(result.col[degree]==0));
  }
  return result;
}

Candidate evaluate_aux(int t,const std::array<int,3>& vertex,
                       const std::array<int,3>& matrices,int boundary,Ranker& ranker) {
  Candidate result;result.vertex=vertex;result.auxiliary_matrices=true;result.boundary=boundary;
  for(int s=0;s<3;s++) result.edge[s]=matrices[s];
  int n=1<<t,maximum_degree=t*(1+*std::max_element(vertex.begin(),vertex.end()));
  assert(maximum_degree<64);
  uint64_t local[4][4]={};
  for(int g=0;g<2;g++) for(int h=0;h<2;h++)
    for(int s=0;s<3;s++) {
      if((s==0&&g)||(s==2&&h)) continue;
      for(int i=0;i<2;i++) for(int j=0;j<2;j++) {
        uint64_t& value=local[2*g+h][2*i+j];
        if(i==j) value^=uint64_t(1)<<vertex[s];
        if(matrices[s]&(1<<(2*i+j)))
          value^=(uint64_t(1)<<vertex[s])^(uint64_t(1)<<(vertex[s]+1));
      }
    }
  std::vector<std::vector<uint64_t>> rows(maximum_degree+1,std::vector<uint64_t>(n));
  std::vector<std::vector<uint64_t>> cols(maximum_degree+1,std::vector<uint64_t>(n));
  for(int g=0;g<n;g++) for(int h=0;h<n;h++) {
    uint64_t dp[4]={1,0,0,1};
    for(int position=0;position<t;position++) {
      uint64_t next[4]={};
      const auto& factor=local[2*((g>>position)&1)+((h>>position)&1)];
      for(int i=0;i<(boundary?2:1);i++) for(int j=0;j<2;j++)
        for(int k=0;k<2;k++) {
          uint64_t polynomial=factor[2*k+j];
          while(polynomial) {
            int shift=__builtin_ctzll(polynomial);polynomial&=polynomial-1;
            next[2*i+j]^=dp[2*i+k]<<shift;
          }
        }
      std::copy(next,next+4,dp);
    }
    uint64_t polynomial=boundary?(dp[1]^dp[2]^dp[3]):dp[0];
    assert((__builtin_parityll(polynomial))==(g==h));
    while(polynomial) {
      int degree=__builtin_ctzll(polynomial);polynomial&=polynomial-1;
      rows[degree][g]|=uint64_t(1)<<h;
      cols[degree][h]|=uint64_t(1)<<g;
    }
  }
  result.row.resize(maximum_degree+1);result.col.resize(maximum_degree+1);
  for(int degree=0;degree<=maximum_degree;degree++) {
    for(int g=0;g<n;g++) {
      result.row[degree]+=ranker.rank(rows[degree][g]);
      result.col[degree]+=ranker.rank(cols[degree][g]);
    }
    uint64_t p=uint64_t(result.row[degree])*result.col[degree];
    result.maximum_product=std::max(result.maximum_product,p);
    result.rank_sum+=std::sqrt((long double)p);
    assert((result.row[degree]==0)==(result.col[degree]==0));
  }
  return result;
}

void output(std::ostream& stream,const Candidate& c,int t) {
  stream<<"{\"t\":"<<t<<",\"polynomial_edges\":"<<(c.polynomial_edges?"true":"false")
        <<",\"auxiliary_matrices\":"<<(c.auxiliary_matrices?"true":"false")
        <<",\"boundary\":"<<c.boundary<<",\"vertex\":[";
  for(int j=0;j<3;j++) stream<<(j?",":"")<<c.vertex[j];
  stream<<"],\"edge\":[";
  for(int j=0;j<9;j++) stream<<(j?",":"")<<c.edge[j];
  stream<<"],\"row_ranks\":[";
  for(size_t j=0;j<c.row.size();j++) stream<<(j?",":"")<<c.row[j];
  stream<<"],\"column_ranks\":[";
  for(size_t j=0;j<c.col.size();j++) stream<<(j?",":"")<<c.col[j];
  stream<<"],\"maximum_rank_product\":"<<c.maximum_product
        <<",\"rank_sum_estimate\":"<<std::setprecision(18)<<c.rank_sum
        <<",\"maximum_product_rho\":"
        <<std::log2((long double)c.maximum_product)/(4*t)
        <<",\"rank_sum_rho\":"<<std::log2(c.rank_sum)/(2*t)<<"}\n";
}

int main(int argc,char** argv) {
  if(argc<4) {
    std::cerr<<"usage: markov_pairing_search screen t output.jsonl [vertex_max]\n"
             <<"   or: markov_pairing_search profile t output.jsonl aR aH aC bRR ... bCC\n";
    return 2;
  }
  std::string mode=argv[1];
  int t=std::stoi(argv[2]);assert(t>=1&&t<=6);
  Ranker ranker(t);
  std::ofstream stream(argv[3]);assert(stream);
  if(mode=="selftest") {
    assert(t<=4);
    uint64_t count=uint64_t(1)<<ranker.n;
    for(uint64_t kernel=0;kernel<count;kernel++)
      assert(ranker.rank(kernel)==ranker.direct(kernel));
    stream<<"{\"status\":\"passed\",\"word_bits\":"<<t<<",\"all_kernels\":"<<count<<"}\n";
    std::cout<<"checked all "<<count<<" kernels"<<std::endl;
    return 0;
  }
  if(mode=="auxscreen") {
    assert(argc==8);
    std::array<int,3> vertex{};
    for(int j=0;j<3;j++) vertex[j]=std::stoi(argv[4+j]);
    int boundary=std::stoi(argv[7]);assert(boundary==0||boundary==1);
    std::vector<Candidate> best;
    for(int code=0;code<4096;code++) {
      std::array<int,3> matrices{code&15,(code>>4)&15,(code>>8)&15};
      auto candidate=evaluate_aux(t,vertex,matrices,boundary,ranker);
      output(stream,candidate,t);best.push_back(candidate);
      std::sort(best.begin(),best.end(),[](const auto& x,const auto& y) {
        if(x.maximum_product!=y.maximum_product) return x.maximum_product<y.maximum_product;
        return x.rank_sum<y.rank_sum;
      });
      if(best.size()>20) best.resize(20);
      if(code%512==511) std::cerr<<"auxiliary checked "<<code+1
        <<" boundary "<<boundary<<" best product "<<best[0].maximum_product
        <<" sum "<<double(best[0].rank_sum)<<std::endl;
    }
    for(const auto& candidate:best) output(std::cout,candidate,t);
    return 0;
  }
  if(mode=="auxprofile") {
    assert(argc==11);
    std::array<int,3> vertex{},matrices{};
    for(int j=0;j<3;j++) vertex[j]=std::stoi(argv[4+j]);
    for(int j=0;j<3;j++) matrices[j]=std::stoi(argv[7+j]);
    auto candidate=evaluate_aux(t,vertex,matrices,std::stoi(argv[10]),ranker);
    output(stream,candidate,t);output(std::cout,candidate,t);return 0;
  }
  if(mode=="profile"||mode=="polyprofile") {
    assert(argc==16);
    std::array<int,3> vertex{};std::array<int,9> edge{};
    for(int j=0;j<3;j++) vertex[j]=std::stoi(argv[4+j]);
    for(int j=0;j<9;j++) edge[j]=std::stoi(argv[7+j]);
    auto candidate=evaluate(t,vertex,edge,ranker,mode=="polyprofile");
    output(stream,candidate,t);output(std::cout,candidate,t);
    return 0;
  }
  if(mode=="polyscreen") {
    assert(argc==7);
    std::array<int,3> vertex{};
    for(int j=0;j<3;j++) vertex[j]=std::stoi(argv[4+j]);
    std::vector<Candidate> best;
    const int alphabet[3]={1,2,7};
    for(int code=0;code<19683;code++) {
      int remaining=code;std::array<int,9> edge{};
      for(int j=0;j<9;j++) {edge[j]=alphabet[remaining%3];remaining/=3;}
      auto candidate=evaluate(t,vertex,edge,ranker,true);
      output(stream,candidate,t);best.push_back(candidate);
      std::sort(best.begin(),best.end(),[](const auto& x,const auto& y) {
        if(x.maximum_product!=y.maximum_product) return x.maximum_product<y.maximum_product;
        return x.rank_sum<y.rank_sum;
      });
      if(best.size()>20) best.resize(20);
      if(code%2000==1999) std::cerr<<"polynomials checked "<<code+1
        <<" best product "<<best[0].maximum_product<<" sum "<<double(best[0].rank_sum)<<std::endl;
    }
    for(const auto& candidate:best) output(std::cout,candidate,t);
    return 0;
  }
  assert(mode=="screen");
  int vertex_max=argc>4?std::stoi(argv[4]):2;
  std::vector<Candidate> best;
  size_t checked=0;
  for(int a=0;a<=vertex_max;a++)
    for(int b=0;b<=vertex_max;b++)
      for(int c=a;c<=vertex_max;c++) {
        if(std::min({a,b,c})!=0) continue;
        for(int mask=0;mask<512;mask++) {
          std::array<int,9> edge{};
          for(int j=0;j<9;j++) edge[j]=(mask>>j)&1;
          auto candidate=evaluate(t,{a,b,c},edge,ranker);
          output(stream,candidate,t);
          best.push_back(candidate);
          std::sort(best.begin(),best.end(),[](const auto& x,const auto& y) {
            if(x.maximum_product!=y.maximum_product) return x.maximum_product<y.maximum_product;
            return x.rank_sum<y.rank_sum;
          });
          if(best.size()>20) best.resize(20);
          checked++;
        }
        std::cerr<<"checked "<<checked<<" vertex "<<a<<","<<b<<","<<c
                 <<" best product "<<best[0].maximum_product
                 <<" sum "<<double(best[0].rank_sum)<<std::endl;
      }
  std::cout<<"{\"checked\":"<<checked<<",\"rank_calls\":"<<ranker.calls
           <<",\"cache_hits\":"<<ranker.hits<<",\"gaussian_ranks\":"<<ranker.gaussian<<"}\n";
  for(const auto& candidate:best) output(std::cout,candidate,t);
}
