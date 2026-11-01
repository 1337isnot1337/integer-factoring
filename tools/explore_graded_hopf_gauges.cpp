// Parameter-dependent, per-leg central gauges of the graded Hopf crossing.
// Every gauge vanishes at t=0, so the z=1 crossing endpoint remains exact.
#define main graded_hopf_base_main
#include "explore_graded_hopf_crossings.cpp"
#undef main
#include <cmath>
#include <iomanip>

struct GaugedScreen:Screen {
  int degree_bound,actual_degree=0,gauge_code;
  std::vector<int> gauge;
  GaugedScreen(int m,int t,const std::vector<int>&columns,int code,size_t memory):
      Screen(m,t,columns,memory),degree_bound(2*maximum+4*t),gauge_code(code),gauge(4) {
    for(int i=0;i<4;++i)gauge[i]=(code>>(2*i))&3;
    const uint64_t choices[4]={0,3,5,6};
    uint64_t outL=choices[gauge[0]],outR=choices[gauge[1]],inL=choices[gauge[2]],inR=choices[gauge[3]];
    size_t bytes=8*size_t(degree_bound+1)*(row_layout.words+column_layout.words);
    if(bytes>memory*1024*1024)throw std::runtime_error("gauge matrix memory exceeds --memory-mib");
    row_matrices.clear();column_matrices.clear();
    row_matrices.assign(degree_bound+1,std::vector<uint64_t>(row_layout.words,0));
    column_matrices.assign(degree_bound+1,std::vector<uint64_t>(column_layout.words,0));
    coefficient_pivots.assign(degree_bound+1,0);
    for(int side=0;side<2;++side)for(int j=0;j<t;++j) {
      Polynomial image={{uint32_t(1)<<(side*width+2*m+j),5},
                        {uint32_t(1)<<((1-side)*width+2*m+j),4}};
      uint64_t qL,qR,crossXY,crossYX;
      if(side==0) {
        qL=multiply_polynomial(outL,5)^multiply_polynomial(inL,5);
        qR=6^multiply_polynomial(outL,4)^multiply_polynomial(inR,4);
        crossXY=3^multiply_polynomial(outL,6);
        crossYX=5^multiply_polynomial(outL,6);
      }else {
        qL=3^multiply_polynomial(outR,4)^multiply_polynomial(inL,4);
        qR=5^multiply_polynomial(outR,5)^multiply_polynomial(inR,5);
        crossXY=multiply_polynomial(outR,6);
        crossYX=6^multiply_polynomial(outR,6);
      }
      for(int i=0;i<m;++i)if((columns[i]>>j)&1) {
        uint32_t xL=uint32_t(1)<<i,yL=uint32_t(1)<<(m+i),xR=xL<<width,yR=yL<<width;
        if(qL)image.push_back({xL|yL,qL});
        if(qR)image.push_back({xR|yR,qR});
        if(crossXY)image.push_back({xL|yR,crossXY});
        if(crossYX)image.push_back({yL|xR,crossYX});
      }
      images[side*width+2*m+j]=std::move(image);
    }
    std::cout<<"{\"event\":\"gauge_start\",\"gauge_code\":"<<code<<",\"gauges\":[";
    for(int i=0;i<4;++i)std::cout<<(i?",":"")<<gauge[i];
    std::cout<<"],\"degree_bound\":"<<degree_bound<<",\"matrix_bytes\":"<<bytes<<"}"<<std::endl;
  }

  void record_gauged(uint32_t output,const Polynomial&polynomial) {
    int x=output&(d-1),y=output>>width,row=x*d+y;
    int rb=row_layout.block[row],rs=row_layout.strides[rb];
    size_t ro=row_layout.offsets[rb]+size_t(row_layout.position[row])*rs;
    bool crossing=false;
    for(auto [input,p]:polynomial) {
      int u=input&(d-1),v=input>>width,col=u*d+v;
      if(row_layout.block[col]!=rb)throw std::runtime_error("gauge weight conservation failed");
      bool endpoint=__builtin_parityll(p);
      if(endpoint!=(x==v&&y==u))throw std::runtime_error("gauge crossing failed");
      crossing|=endpoint;
      if(p>>(degree_bound+1))throw std::runtime_error("gauge degree bound failed");
      actual_degree=std::max(actual_degree,63-__builtin_clzll(p));
      uint64_t independent=p;
      while(independent){int lead=63-__builtin_clzll(independent);if(!coefficient_pivots[lead]){coefficient_pivots[lead]=independent;break;}independent^=coefficient_pivots[lead];}
      uint64_t token=mix64(mix64(uint64_t(output)^(uint64_t(input)<<32))^p);
      digest_xor^=token;digest_sum+=token;++nonzero_entries;
      int cr=x*d+u,cc=y*d+v,cb=column_layout.block[cr];
      if(column_layout.block[cc]!=2*maximum-cb)throw std::runtime_error("gauge column block failed");
      size_t co=column_layout.offsets[cb]+size_t(column_layout.position[cr])*column_layout.strides[cb];
      int bitr=row_layout.position[col],bitc=column_layout.position[cc];
      while(p){int i=__builtin_ctzll(p);p&=p-1;row_matrices[i][ro+(bitr>>6)]^=uint64_t(1)<<(bitr&63);column_matrices[i][co+(bitc>>6)]^=uint64_t(1)<<(bitc&63);}
    }
    if(!crossing)throw std::runtime_error("gauge crossing permutation missing");
    ++finished_rows;
  }

  void visit_gauged(int variable,uint32_t output,const Polynomial&polynomial) {
    if(variable==2*width){record_gauged(output,polynomial);return;}
    visit_gauged(variable+1,output,polynomial);
    Polynomial product=multiply(polynomial,images[variable]);
    visit_gauged(variable+1,output|(uint32_t(1)<<variable),product);
  }

  void run_gauged() {
    auto began=std::chrono::steady_clock::now();visit_gauged(0,0,{{0,1}});
    if(finished_rows!=uint64_t(n))throw std::runtime_error("gauge output basis incomplete");
    uint64_t largest=0;double cost=0;int dimension=0;
    for(uint64_t value:coefficient_pivots)dimension+=value!=0;
    std::cout<<"{\"event\":\"gauge_profile\",\"gauge_code\":"<<gauge_code<<",\"coefficient_dimension\":"<<dimension<<",\"profile\":[";
    bool comma=false;
    for(int i=0;i<=actual_degree;++i) {
      uint64_t r=0,c=0;
      for(size_t b=0;b<row_layout.sizes.size();++b)r+=binary_rank(row_matrices[i].data()+row_layout.offsets[b],row_layout.sizes[b],row_layout.strides[b]);
      for(size_t b=0;b<column_layout.sizes.size();++b)c+=binary_rank(column_matrices[i].data()+column_layout.offsets[b],column_layout.sizes[b],column_layout.strides[b]);
      if(r){std::cout<<(comma?",":"")<<"["<<i<<","<<r<<","<<c<<"]";comma=true;}
      largest=std::max(largest,r*c);cost+=std::sqrt(double(r*c));
    }
    auto ms=std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now()-began).count();
    std::cout<<"],\"maximum_rank_product\":"<<largest<<",\"rank_sum_approx\":"<<std::setprecision(17)<<cost
             <<",\"matrix_digest_xor\":\""<<digest_xor<<"\",\"matrix_digest_sum\":\""<<digest_sum
             <<"\",\"nonzero_polynomial_entries\":"<<nonzero_entries<<",\"milliseconds\":"<<ms<<"}"<<std::endl;
  }
};

int main(int argc,char**argv) {
  try {
    int m=2,t=1,code=0;std::string columns_text="1,1";bool sweep=false;size_t memory=2048;
    for(int i=1;i<argc;++i){std::string arg=argv[i];
      if(arg=="--horizontal"&&i+1<argc)m=std::stoi(argv[++i]);
      else if(arg=="--center"&&i+1<argc)t=std::stoi(argv[++i]);
      else if(arg=="--columns"&&i+1<argc)columns_text=argv[++i];
      else if(arg=="--gauge"&&i+1<argc)code=std::stoi(argv[++i]);
      else if(arg=="--sweep")sweep=true;
      else if(arg=="--memory-mib"&&i+1<argc)memory=std::stoull(argv[++i]);
      else throw std::runtime_error("unrecognized argument: "+arg);
    }
    if(m<1||t<1||2*m+t>7||code<0||code>255)throw std::runtime_error("require positive dimensions, width<=7, and gauge code0..255");
    std::vector<int>columns;std::stringstream ss(columns_text);std::string item;
    while(std::getline(ss,item,','))columns.push_back(std::stoi(item));
    if(int(columns.size())!=m)throw std::runtime_error("incorrect code length");
    for(int value:columns)if(value<=0||value>=(1<<t))throw std::runtime_error("invalid code column");
    for(int gauge=sweep?0:code;gauge<=(sweep?255:code);++gauge){GaugedScreen screen(m,t,columns,gauge,memory);screen.run_gauged();}
    return 0;
  }catch(const std::exception&error){std::cerr<<error.what()<<std::endl;return 2;}
}
