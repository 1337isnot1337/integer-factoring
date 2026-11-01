// Exact coefficient ranks for graded restricted-Heisenberg Hopf crossings.
// Generator substitutions are independently checked against the literal
// coproduct/antipode construction in explore_graded_hopf_crossings.py.
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using Term=std::pair<uint32_t,uint64_t>;
using Polynomial=std::vector<Term>;

uint64_t multiply_polynomial(uint64_t a,uint64_t b) {
  uint64_t out=0;
  while(b) {if(b&1)out^=a;b>>=1;a<<=1;}
  return out;
}

uint64_t mix64(uint64_t value) {
  value=(value^(value>>30))*0xbf58476d1ce4e5b9ULL;
  value=(value^(value>>27))*0x94d049bb133111ebULL;
  return value^(value>>31);
}

struct BlockLayout {
  std::vector<int> block,position,sizes,strides;
  std::vector<size_t> offsets;
  size_t words=0;
  BlockLayout(int d,const std::vector<int>&weight,int maximum,bool difference) {
    int count=2*maximum+1;
    block.resize(d*d);position.resize(d*d);sizes.assign(count,0);
    for(int a=0;a<d;++a)for(int b=0;b<d;++b) {
      int id=difference?weight[a]-weight[b]+maximum:weight[a]+weight[b];
      int pair=a*d+b;block[pair]=id;position[pair]=sizes[id]++;
    }
    for(int size:sizes){offsets.push_back(words);strides.push_back((size+63)/64);words+=size_t(size)*((size+63)/64);}
  }
};

int binary_rank(uint64_t*rows,int n,int stride) {
  std::vector<int> pivot(n,-1);int rank=0;
  for(int i=0;i<n;++i) {
    uint64_t*row=rows+size_t(i)*stride;int block=stride-1;
    while(block>=0) {
      while(block>=0&&!row[block])--block;
      if(block<0)break;
      int lead=64*block+63-__builtin_clzll(row[block]);
      if(pivot[lead]<0){pivot[lead]=i;++rank;break;}
      const uint64_t*other=rows+size_t(pivot[lead])*stride;
      for(int j=0;j<=block;++j)row[j]^=other[j];
    }
  }
  return rank;
}

struct Screen {
  int m,t,width,d,n,maximum;
  std::vector<int> columns,weight;
  std::vector<Polynomial> images;
  std::vector<uint64_t> accumulator;
  std::vector<uint32_t> marks;
  uint32_t stamp=0;
  BlockLayout row_layout,column_layout;
  std::vector<std::vector<uint64_t>> row_matrices,column_matrices;
  uint64_t digest_xor=0,digest_sum=0,nonzero_entries=0,finished_rows=0;
  std::vector<uint64_t> coefficient_pivots;

  static std::vector<int> weights(int width,int horizontal) {
    std::vector<int> result(1<<width);
    int mask=(1<<(2*horizontal))-1;
    for(int x=0;x<(1<<width);++x)result[x]=__builtin_popcount(unsigned(x&mask))+2*__builtin_popcount(unsigned(x>>(2*horizontal)));
    return result;
  }

  Screen(int horizontal,int center,std::vector<int> code,size_t memory_mib):
      m(horizontal),t(center),width(2*m+t),d(1<<width),n(d*d),maximum(2*m+2*t),
      columns(std::move(code)),weight(weights(width,m)),accumulator(n),marks(n,0),
      row_layout(d,weight,maximum,false),column_layout(d,weight,maximum,true) {
    coefficient_pivots.assign(maximum+1,0);
    size_t bytes=8*size_t(maximum+1)*(row_layout.words+column_layout.words);
    if(bytes>memory_mib*1024*1024)throw std::runtime_error("matrix memory exceeds --memory-mib guard");
    row_matrices.assign(maximum+1,std::vector<uint64_t>(row_layout.words,0));
    column_matrices.assign(maximum+1,std::vector<uint64_t>(column_layout.words,0));
    std::cout<<"{\"event\":\"start\",\"horizontal_pairs\":"<<m<<",\"center_bits\":"<<t
             <<",\"labels\":"<<d<<",\"label_bits\":"<<width<<",\"matrix_bytes\":"<<bytes<<",\"columns\":[";
    for(size_t i=0;i<columns.size();++i)std::cout<<(i?",":"")<<columns[i];
    std::cout<<"]}"<<std::endl;
    for(int side=0;side<2;++side)for(int variable=0;variable<width;++variable) {
      Polynomial image;
      uint32_t same=uint32_t(1)<<(side*width+variable),other=uint32_t(1)<<((1-side)*width+variable);
      if(variable<2*m)image={{same,3},{other,2}};
      else {
        image={{same,5},{other,4}};int center_bit=variable-2*m;
        for(int i=0;i<m;++i)if((columns[i]>>center_bit)&1) {
          uint32_t xL=uint32_t(1)<<i,yL=uint32_t(1)<<(m+i),xR=xL<<width,yR=yL<<width;
          if(side==0){image.push_back({xL|yR,3});image.push_back({yL|xR,5});image.push_back({xR|yR,6});}
          else {image.push_back({xL|yL,3});image.push_back({yL|xR,6});image.push_back({xR|yR,5});}
        }
      }
      images.push_back(std::move(image));
    }
  }

  Polynomial multiply(const Polynomial&a,const Polynomial&b) {
    if(++stamp==0){std::fill(marks.begin(),marks.end(),0);stamp=1;}
    std::vector<uint32_t> touched;
    for(auto [x,p]:a)for(auto [y,q]:b)if(!(x&y)) {
      uint32_t monomial=x|y;
      if(marks[monomial]!=stamp){marks[monomial]=stamp;accumulator[monomial]=0;touched.push_back(monomial);}
      accumulator[monomial]^=multiply_polynomial(p,q);
    }
    Polynomial out;out.reserve(touched.size());
    for(uint32_t monomial:touched)if(accumulator[monomial])out.push_back({monomial,accumulator[monomial]});
    return out;
  }

  void record(uint32_t output,const Polynomial&polynomial) {
    int x=output&(d-1),y=output>>width,row=x*d+y;
    int row_block=row_layout.block[row],row_stride=row_layout.strides[row_block];
    size_t row_offset=row_layout.offsets[row_block]+size_t(row_layout.position[row])*row_stride;
    bool crossing=false;
    for(auto [input,p]:polynomial) {
      int u=input&(d-1),v=input>>width,col=u*d+v;
      if(row_layout.block[col]!=row_block)throw std::runtime_error("weight conservation failed");
      bool endpoint=__builtin_parityll(p);
      if(endpoint!=(x==v&&y==u))throw std::runtime_error("crossing endpoint failed");
      crossing|=endpoint;
      if(p>>(maximum+1))throw std::runtime_error("unexpected coefficient degree");
      uint64_t independent=p;
      while(independent) {
        int lead=63-__builtin_clzll(independent);
        if(!coefficient_pivots[lead]){coefficient_pivots[lead]=independent;break;}
        independent^=coefficient_pivots[lead];
      }
      uint64_t token=mix64(mix64(uint64_t(output)^(uint64_t(input)<<32))^p);
      digest_xor^=token;digest_sum+=token;++nonzero_entries;
      int column_row=x*d+u,column_col=y*d+v;
      int column_block=column_layout.block[column_row];
      if(column_layout.block[column_col]!=2*maximum-column_block)throw std::runtime_error("column energy difference failed");
      int column_stride=column_layout.strides[column_block];
      size_t column_offset=column_layout.offsets[column_block]+size_t(column_layout.position[column_row])*column_stride;
      int row_bit=row_layout.position[col],column_bit=column_layout.position[column_col];
      while(p) {
        int degree=__builtin_ctzll(p);p&=p-1;
        row_matrices[degree][row_offset+(row_bit>>6)]^=uint64_t(1)<<(row_bit&63);
        column_matrices[degree][column_offset+(column_bit>>6)]^=uint64_t(1)<<(column_bit&63);
      }
    }
    if(!crossing)throw std::runtime_error("crossing permutation entry absent");
    ++finished_rows;
  }

  void visit(int variable,uint32_t output,const Polynomial&polynomial) {
    if(variable==2*width){record(output,polynomial);return;}
    visit(variable+1,output,polynomial);
    Polynomial product=multiply(polynomial,images[variable]);
    visit(variable+1,output|(uint32_t(1)<<variable),product);
  }

  uint64_t rank_combination(const std::vector<std::vector<uint64_t>>&matrices,
                            const BlockLayout&layout,uint32_t mask) const {
    std::vector<uint64_t> combined(layout.words,0);
    for(int degree=0;degree<=maximum;++degree)if((mask>>degree)&1)
      for(size_t j=0;j<combined.size();++j)combined[j]^=matrices[degree][j];
    uint64_t rank=0;
    for(size_t b=0;b<layout.sizes.size();++b)
      rank+=binary_rank(combined.data()+layout.offsets[b],layout.sizes[b],layout.strides[b]);
    return rank;
  }

  void run(uint32_t span_start,uint32_t span_count) {
    auto began=std::chrono::steady_clock::now();
    visit(0,0,{{0,1}});
    if(finished_rows!=uint64_t(n))throw std::runtime_error("incomplete output basis");
    for(uint64_t pivot:coefficient_pivots)if(!pivot)throw std::runtime_error("coefficient family is linearly dependent");
    uint64_t maximum_product=0;
    std::vector<uint64_t> basis_row(maximum+1),basis_column(maximum+1);
    for(int degree=0;degree<=maximum;++degree) {
      uint64_t r=rank_combination(row_matrices,row_layout,uint32_t(1)<<degree);
      uint64_t c=rank_combination(column_matrices,column_layout,uint32_t(1)<<degree);
      basis_row[degree]=r;basis_column[degree]=c;
      maximum_product=std::max(maximum_product,r*c);
      std::cout<<"{\"event\":\"coefficient\",\"degree\":"<<degree<<",\"row_rank\":"<<r
               <<",\"column_rank\":"<<c<<",\"product\":"<<r*c<<"}"<<std::endl;
    }
    if(span_count) {
      uint32_t limit=uint32_t(1)<<(maximum+1);
      if(!span_start||span_start>=limit||span_count>limit-span_start)throw std::runtime_error("invalid coefficient-span interval");
      for(uint32_t mask=span_start;mask<span_start+span_count;++mask) {
        bool single=!(mask&(mask-1));int bit=__builtin_ctz(mask);
        uint64_t r=single?basis_row[bit]:rank_combination(row_matrices,row_layout,mask);
        uint64_t c=single?basis_column[bit]:rank_combination(column_matrices,column_layout,mask);
        std::cout<<"{\"event\":\"recombination\",\"mask\":"<<mask<<",\"row_rank\":"<<r
                 <<",\"column_rank\":"<<c<<",\"product\":"<<r*c<<"}"<<std::endl;
      }
    }
    auto ms=std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now()-began).count();
    std::cout<<"{\"event\":\"complete\",\"nonzero_polynomial_entries\":"<<nonzero_entries
             <<",\"matrix_digest_xor\":\""<<digest_xor<<"\",\"matrix_digest_sum\":\""<<digest_sum
             <<"\",\"maximum_rank_product\":"<<maximum_product<<",\"span_start\":"<<span_start
             <<",\"span_count\":"<<span_count<<",\"milliseconds\":"<<ms<<"}"<<std::endl;
  }
};

int main(int argc,char**argv) {
  try {
    int m=1,t=1;std::string code="1";size_t memory=2048;uint32_t span_start=1,span_count=0;bool full_span=false;
    for(int i=1;i<argc;++i) {
      std::string arg=argv[i];
      if(arg=="--horizontal"&&i+1<argc)m=std::stoi(argv[++i]);
      else if(arg=="--center"&&i+1<argc)t=std::stoi(argv[++i]);
      else if(arg=="--columns"&&i+1<argc)code=argv[++i];
      else if(arg=="--memory-mib"&&i+1<argc)memory=std::stoull(argv[++i]);
      else if(arg=="--span")full_span=true;
      else if(arg=="--span-start"&&i+1<argc)span_start=std::stoul(argv[++i]);
      else if(arg=="--span-count"&&i+1<argc)span_count=std::stoul(argv[++i]);
      else throw std::runtime_error("unrecognized argument: "+arg);
    }
    if(m<1||t<1||2*m+t>8||2*m+2*t>30)throw std::runtime_error("require positive dimensions and label width at most eight");
    std::stringstream stream(code);std::string item;std::vector<int> columns;
    while(std::getline(stream,item,','))columns.push_back(std::stoi(item));
    if(int(columns.size())!=m)throw std::runtime_error("one column is required per horizontal pair");
    for(int value:columns)if(value<=0||value>=(1<<t))throw std::runtime_error("invalid binary code column");
    if(full_span)span_count=(uint32_t(1)<<(2*m+2*t+1))-span_start;
    Screen screen(m,t,columns,memory);screen.run(span_start,span_count);return 0;
  }catch(const std::exception&error){std::cerr<<error.what()<<std::endl;return 2;}
}
