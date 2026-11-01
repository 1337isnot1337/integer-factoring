// Exact local Smith spectra for polynomial crossing MPOs over F_2.
#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <string>
#include <vector>

using Polynomial=__uint128_t;
static int trailing_zeros(Polynomial a) {
  assert(a);
  uint64_t low=static_cast<uint64_t>(a);
  return low?__builtin_ctzll(low):64+__builtin_ctzll(static_cast<uint64_t>(a>>64));
}
static int polynomial_parity(Polynomial a) {
  return __builtin_parityll(static_cast<uint64_t>(a))^__builtin_parityll(static_cast<uint64_t>(a>>64));
}
using Spectrum=std::array<int,4>;
using Binary8=std::array<unsigned,8>;
static Binary8 binary_product(const Binary8& a,const Binary8& b) {
  Binary8 out{};
  for(int i=0;i<8;i++)for(int j=0;j<8;j++)if((a[i]>>j)&1)out[i]^=b[j];
  return out;
}
static int binary_rank4(unsigned matrix) {
  std::array<unsigned,4> rows;
  for(int i=0;i<4;i++)rows[i]=(matrix>>(4*i))&15;
  int out=0;
  for(int j=0;j<4;j++) {
    int pivot=out;
    while(pivot<4&&!((rows[pivot]>>j)&1))pivot++;
    if(pivot==4)continue;
    std::swap(rows[out],rows[pivot]);
    for(int i=out+1;i<4;i++)if((rows[i]>>j)&1)rows[i]^=rows[out];
    out++;
  }
  return out;
}
static Polynomial multiply(Polynomial a,Polynomial b,int precision) {
  Polynomial out=0,mask=(Polynomial(1)<<precision)-1;
  while(b) {int j=trailing_zeros(b);out^=a<<j;b&=b-1;}
  return out&mask;
}
static Polynomial inverse_unit(Polynomial a,int precision) {
  assert(a&1);Polynomial b=1;
  for(int i=1;i<precision;i++)
    if((multiply(a,b,precision)>>i)&1)b^=Polynomial(1)<<i;
  assert(multiply(a,b,precision)==1);return b;
}
static std::pair<bool,Spectrum> smith_at_precision(std::vector<std::array<Polynomial,4>> a,int precision) {
  assert(precision<128);
  Polynomial mask=(Polynomial(1)<<precision)-1;
  for(auto& row:a)for(auto& value:row)value&=mask;
  Spectrum out{};
  for(int k=0;k<4;k++) {
    int bi=-1,bj=-1,bv=precision;
    for(size_t i=k;i<a.size();i++)for(int j=k;j<4;j++)if(a[i][j]) {
      int v=trailing_zeros(a[i][j]);
      if(v<bv){bv=v;bi=i;bj=j;}
    }
    if(bi<0)return {false,out};
    out[k]=bv;
    std::swap(a[k],a[bi]);
    for(auto& row:a)std::swap(row[k],row[bj]);
    Polynomial inverse=inverse_unit(a[k][k]>>bv,precision);
    for(int j=k;j<4;j++)a[k][j]=multiply(a[k][j],inverse,precision);
    assert(a[k][k]==(Polynomial(1)<<bv));
    for(size_t i=k+1;i<a.size();i++) {
      assert(!a[i][k]||trailing_zeros(a[i][k])>=bv);
      Polynomial factor=a[i][k]>>bv;
      for(int j=k+1;j<4;j++)a[i][j]^=multiply(factor,a[k][j],precision);
      a[i][k]=0;
    }
  }
  return {true,out};
}
static Spectrum smith(const std::vector<std::array<Polynomial,4>>& a,int degree) {
  int maximum=4*degree+1;
  assert(maximum<128);
  int precision=degree+1;
  while(true) {
    auto result=smith_at_precision(a,precision);
    if(result.first)return result.second;
    assert(precision<maximum);
    precision=std::min(2*precision,maximum);
  }
}
static Polynomial reverse(Polynomial a,int degree) {
  Polynomial out=0;
  for(int j=0;j<=degree;j++)if((a>>j)&1)out|=Polynomial(1)<<(degree-j);
  return out;
}
static std::array<Spectrum,8> profiles(const std::array<std::array<Polynomial,4>,16>& local,int degree,int height=0,int ramification=1,std::array<int,4> physical={}) {
  std::array<Spectrum,8> out;
  degree*=ramification;
  int offset=std::abs(height);
  for(int h:physical)offset+=std::abs(h);
  int bound_degree=degree+2*offset;
  assert(4*bound_degree+1<128);
  for(int cut=0;cut<2;cut++)for(int side=0;side<2;side++) {
    std::vector<std::array<Polynomial,4>> matrix(16),infinity_matrix(16);
    for(int x=0;x<2;x++)for(int y=0;y<2;y++)for(int u=0;u<2;u++)for(int v=0;v<2;v++) {
      int p=cut?2*x+u:2*x+y,q=cut?2*y+v:2*u+v;
      for(int ab=0;ab<4;ab++) {
        int a=ab/2,b=ab%2,shift=offset+(a-b)*height;
        shift+=(2*x-1)*physical[0]+(2*y-1)*physical[1]+(2*u-1)*physical[2]+(2*v-1)*physical[3];
        auto raw=local[8*x+4*y+2*u+v][ab];
        Polynomial entry=0;
        for(int k=0;k<=degree/ramification;k++)if((raw>>k)&1)entry|=Polynomial(1)<<(k*ramification);
        matrix[4*(side?q:p)+ab][side?p:q]=entry<<shift;
        infinity_matrix[4*(side?q:p)+ab][side?p:q]=reverse(entry,degree)<<shift;
      }
    }
    Spectrum zero=smith(matrix,bound_degree);
    for(auto& value:zero)value-=offset;
    out[4*cut+2*side]=zero;
    Spectrum infinity=smith(infinity_matrix,bound_degree);
    for(auto& value:infinity)value=degree+offset-value;
    out[4*cut+2*side+1]=infinity;
  }
  return out;
}
static std::array<Spectrum,8> lattice_profiles(const std::array<std::array<Polynomial,4>,16>& local,int degree,std::array<int,3> lattice,int ramification=1,std::array<int,4> physical={}) {
  int a=lattice[0],b=lattice[1],h=a+b;
  Polynomial p=lattice[2],p_squared=multiply(p,p,127);
  degree*=ramification;
  int offset=0;for(int value:physical)offset+=std::abs(value);
  int bound_degree=degree+2*h+2*offset;
  assert(4*bound_degree+1<128);
  std::array<std::array<std::array<Polynomial,4>,16>,2> changed{};
  for(int labels=0;labels<16;labels++) {
    std::array<Polynomial,4> raw{};
    for(int ab=0;ab<4;ab++)for(int k=0;k<=degree/ramification;k++)
      if((local[labels][ab]>>k)&1)raw[ab]|=Polynomial(1)<<(k*ramification);
    for(int infinity=0;infinity<2;infinity++) {
      auto entry=raw;
      if(infinity)for(auto& value:entry)value=reverse(value,degree);
      Polynomial coupling=multiply(p,entry[2],127)<<b;
      changed[infinity][labels][0]=(entry[0]<<h)^coupling;
      changed[infinity][labels][1]=(multiply(p,entry[0]^entry[3],127)<<a)^multiply(p_squared,entry[2],127)^(entry[1]<<(2*a));
      changed[infinity][labels][2]=entry[2]<<(2*b);
      changed[infinity][labels][3]=(entry[3]<<h)^coupling;
    }
  }
  std::array<Spectrum,8> out{};
  for(int cut=0;cut<2;cut++)for(int side=0;side<2;side++)for(int infinity=0;infinity<2;infinity++) {
    std::vector<std::array<Polynomial,4>> matrix(16);
    for(int x=0;x<2;x++)for(int y=0;y<2;y++)for(int u=0;u<2;u++)for(int v=0;v<2;v++) {
      int left=cut?2*x+u:2*x+y,right=cut?2*y+v:2*u+v;
      int shift=offset+(2*x-1)*physical[0]+(2*y-1)*physical[1]+(2*u-1)*physical[2]+(2*v-1)*physical[3];
      for(int ab=0;ab<4;ab++)matrix[4*(side?right:left)+ab][side?left:right]=changed[infinity][8*x+4*y+2*u+v][ab]<<shift;
    }
    auto values=smith(matrix,bound_degree);
    for(auto& value:values)value=infinity?degree+h+offset-value:value-h-offset;
    out[4*cut+2*side+infinity]=values;
  }
  return out;
}
struct Value {double bound=1e100,s=1;int row=-1,column=-1;};
static std::pair<double,double> spectral_value(const Spectrum& degrees,const Spectrum& orders) {
  static std::map<std::array<int,8>,std::pair<double,double>> cache;
  std::array<int,8> key;
  for(int i=0;i<4;i++){key[i]=degrees[i];key[i+4]=orders[i];}
  auto old=cache.find(key);
  if(old!=cache.end())return old->second;
  std::array<int,16> exponents;
  for(int i=0;i<4;i++)for(int j=0;j<4;j++)exponents[4*i+j]=degrees[i]-orders[j];
  auto derivative=[&](double argument) {double out=0;for(int e:exponents)out+=e*std::pow(argument,e);return out;};
  double bests=1;
  if(derivative(1)>=0)bests=1;
  else if(derivative(16)<=0)bests=1048576;
  else {
    double low=1,high=16;
    for(int iteration=0;iteration<70;iteration++) {
      double middle=(low+high)/2;
      if(derivative(middle)<0)low=middle;else high=middle;
    }
    bests=std::round((low+high)*524288)/1048576.;
  }
  double best=0;for(int e:exponents)best+=std::pow(bests,e);
  return cache[key]={best,bests};
}
static Value optimize(const std::array<Spectrum,8>& spectrum) {
  Value result;
  for(int orientation=0;orientation<2;orientation++)
  for(int side_r=0;side_r<2;side_r++)for(int side_c=0;side_c<2;side_c++) {
    int r=4*orientation+2*side_r+1,c=4*(1-orientation)+2*side_c;
    auto value=spectral_value(spectrum[r],spectrum[c]);
    if(value.first<result.bound)result={value.first,value.second,r,c};
  }
  return result;
}
struct GaugedValue {
  Value value;
  std::array<Spectrum,8> spectrum{};
  int infinity_height=0,zero_height=0;
  std::array<int,4> infinity_physical{},zero_physical{};
  int infinity_basis=9,zero_basis=9;
  std::array<int,3> infinity_lattice{},zero_lattice{};
};
struct GaugeChoice {
  Spectrum spectrum;
  int height,index;
  std::array<int,4> physical;
  int basis=9;
  std::array<int,3> lattice{};
};
static std::array<std::array<Polynomial,4>,16> auxiliary_conjugate(const std::array<std::array<Polynomial,4>,16>& local,int basis) {
  int inverse=((basis>>3)&1)|(((basis>>1)&1)<<1)|(((basis>>2)&1)<<2)|((basis&1)<<3);
  assert(((basis&1)*((basis>>3)&1))^(((basis>>1)&1)*((basis>>2)&1)));
  std::array<std::array<Polynomial,4>,16> out{};
  for(int labels=0;labels<16;labels++)for(int a=0;a<2;a++)for(int b=0;b<2;b++)
  for(int c=0;c<2;c++)for(int d=0;d<2;d++)
    if(((basis>>(2*a+c))&1)&&((inverse>>(2*d+b))&1))out[labels][2*a+b]^=local[labels][2*c+d];
  return out;
}
static std::vector<GaugeChoice> pareto(const std::map<Spectrum,GaugeChoice>& source,bool infinity) {
  std::vector<GaugeChoice> result;
  for(const auto& entry:source) {
    bool dominated=false;
    for(const auto& other:source) {
      if(entry.first==other.first)continue;
      bool better=true;
      for(int i=0;i<4;i++)if(infinity?other.first[i]>entry.first[i]:other.first[i]<entry.first[i])better=false;
      if(better){dominated=true;break;}
    }
    if(!dominated)result.push_back(entry.second);
  }
  return result;
}
static GaugedValue gauge_optimize(const std::array<std::array<Polynomial,4>,16>& local,int degree,int maximum_height,int ramification=1,int maximum_physical=0,bool independent_bases=false) {
  std::array<std::map<Spectrum,GaugeChoice>,2> infinities,zeros;
  int radius=maximum_physical;
  std::vector<int> bases=independent_bases?std::vector<int>{6,7,9,11,13,14}:std::vector<int>{9};
  for(int basis:bases) {
  auto changed=auxiliary_conjugate(local,basis);
  for(int a=-radius;a<=radius;a++)for(int b=-radius;b<=radius;b++)
  for(int c=-radius;c<=radius;c++)for(int d=-radius;d<=radius;d++)
  for(int h=-maximum_height;h<=maximum_height;h++) {
    std::array<int,4> physical{a,b,c,d};
    auto current=profiles(changed,degree,h,ramification,physical);
    for(int cut=0;cut<2;cut++)for(int side=0;side<2;side++) {
      int zero=4*cut+2*side,infinity=zero+1;
      zeros[cut].emplace(current[zero],GaugeChoice{current[zero],h,zero,physical,basis});
      infinities[cut].emplace(current[infinity],GaugeChoice{current[infinity],h,infinity,physical,basis});
    }
  }
  }
  GaugedValue best;
  for(int cut=0;cut<2;cut++) {
    auto high=pareto(infinities[cut],true),low=pareto(zeros[1-cut],false);
    for(const auto& first:high)for(const auto& second:low) {
      auto value=spectral_value(first.spectrum,second.spectrum);
      if(value.first<best.value.bound) {
        best.value={value.first,value.second,first.index,second.index};
        best.spectrum[first.index]=first.spectrum;best.spectrum[second.index]=second.spectrum;
        best.infinity_height=first.height;best.zero_height=second.height;
        best.infinity_physical=first.physical;best.zero_physical=second.physical;
        best.infinity_basis=first.basis;best.zero_basis=second.basis;
      }
    }
  }
  return best;
}
static GaugedValue lattice_optimize(const std::array<std::array<Polynomial,4>,16>& local,int degree,int maximum_depth,int ramification=1,int maximum_physical=0) {
  std::array<std::map<Spectrum,GaugeChoice>,2> infinities,zeros;
  int radius=maximum_physical;
  for(int depth=0;depth<=maximum_depth;depth++)for(int a=0;a<=depth;a++) {
    int b=depth-a;
    for(int p=0;p<(1<<b);p++) {
      if(a&&b&&!(p&1))continue;
      std::array<int,3> lattice{a,b,p};
      for(int hx=-radius;hx<=radius;hx++)for(int hy=-radius;hy<=radius;hy++)
      for(int hu=-radius;hu<=radius;hu++)for(int hv=-radius;hv<=radius;hv++) {
        std::array<int,4> physical{hx,hy,hu,hv};
        auto current=lattice_profiles(local,degree,lattice,ramification,physical);
        for(int cut=0;cut<2;cut++)for(int side=0;side<2;side++) {
          int zero=4*cut+2*side,infinity=zero+1;
          zeros[cut].emplace(current[zero],GaugeChoice{current[zero],0,zero,physical,9,lattice});
          infinities[cut].emplace(current[infinity],GaugeChoice{current[infinity],0,infinity,physical,9,lattice});
        }
      }
    }
  }
  GaugedValue best;
  for(int cut=0;cut<2;cut++) {
    auto high=pareto(infinities[cut],true),low=pareto(zeros[1-cut],false);
    for(const auto& first:high)for(const auto& second:low) {
      auto value=spectral_value(first.spectrum,second.spectrum);
      if(value.first<best.value.bound) {
        best.value={value.first,value.second,first.index,second.index};
        best.spectrum[first.index]=first.spectrum;best.spectrum[second.index]=second.spectrum;
        best.infinity_physical=first.physical;best.zero_physical=second.physical;
        best.infinity_lattice=first.lattice;best.zero_lattice=second.lattice;
      }
    }
  }
  return best;
}
static std::string gauge_physical_record(const GaugedValue& selected) {
  std::string out=",\"infinity_physical\":[";
  for(int i=0;i<4;i++)out+=(i?",":"")+std::to_string(selected.infinity_physical[i]);
  out+="],\"zero_physical\":[";
  for(int i=0;i<4;i++)out+=(i?",":"")+std::to_string(selected.zero_physical[i]);
  out+="],\"infinity_basis\":"+std::to_string(selected.infinity_basis)+",\"zero_basis\":"+std::to_string(selected.zero_basis)+",\"infinity_lattice\":[";
  for(int i=0;i<3;i++)out+=(i?",":"")+std::to_string(selected.infinity_lattice[i]);
  out+="],\"zero_lattice\":[";
  for(int i=0;i<3;i++)out+=(i?",":"")+std::to_string(selected.zero_lattice[i]);
  return out+"]";
}
int main(int argc,char** argv) {
  assert(argc>=2&&argc<=9);double best=1e100;int checked=0;std::vector<std::string> records;
  std::string mode=argc>=3?argv[2]:"pairing";
  bool lattice_mode=mode.rfind("lattice-",0)==0;
  bool physical_mode=lattice_mode||mode.rfind("physical-",0)==0;
  bool gauged=physical_mode||mode=="gauged-pairing"||mode=="gauged-loops"||mode=="gauged-resolvents"||mode=="gauged-linear-factor"||mode=="gauged-ket-column";
  int maximum_height=argc>=4?std::stoi(argv[3]):4;
  int ramification=argc>=5?std::stoi(argv[4]):1;
  assert(ramification>=1&&ramification<=8);
  int maximum_physical=physical_mode?(argc>=6?std::stoi(argv[5]):1):0;
  bool independent_bases=argc>=7?std::stoi(argv[6]):false;
  auto optimize_local=[&](const auto& local,int degree) {
    return lattice_mode?lattice_optimize(local,degree,maximum_height,ramification,maximum_physical):gauge_optimize(local,degree,maximum_height,ramification,maximum_physical,independent_bases);
  };
  assert(mode=="pairing"||mode=="loops"||mode=="resolvents"||mode=="linear-factor"||gauged);
  if(mode=="physical-file"||mode=="lattice-file") {
    assert(argc==9);
    int degree=std::stoi(argv[7]);
    std::ifstream input(argv[8]);assert(input);
    std::array<std::array<Polynomial,4>,16> local{};
    for(auto& entries:local)for(auto& value:entries) {
      uint64_t raw;input>>raw;assert(input);value=raw;
    }
    auto selected=optimize_local(local,degree);auto spectrum=selected.spectrum;auto value=selected.value;
    checked=1;best=value.bound;
    std::string record="{\"degree\":"+std::to_string(degree)+gauge_physical_record(selected)+",\"infinity_height\":"+std::to_string(selected.infinity_height)+",\"zero_height\":"+std::to_string(selected.zero_height)+",\"s_numerator\":"+std::to_string(std::llround(value.s*1048576))+",\"s_denominator\":1048576,\"infinity_profile_index\":"+std::to_string(value.row)+",\"zero_profile_index\":"+std::to_string(value.column)+",\"degrees_at_infinity\":[";
    for(int i=0;i<4;i++)record+=(i?",":"")+std::to_string(spectrum[value.row][i]);
    record+="],\"orders_at_zero\":[";
    for(int i=0;i<4;i++)record+=(i?",":"")+std::to_string(spectrum[value.column][i]);
    record+="]}";records.push_back(record);
  } else if((mode=="physical-quadratic-bridge"||mode=="lattice-quadratic-bridge")) {
    for(unsigned v=1;v<16;v++) {
      unsigned swapped=(v&9)|((v&2)<<1)|((v&4)>>1);
      std::vector<unsigned> kernel;
      for(unsigned r=0;r<16;r++)if(!__builtin_parity(r&v)&&!__builtin_parity(r&swapped))kernel.push_back(r);
      for(unsigned r0:kernel)for(unsigned r1:kernel)for(unsigned r2:kernel)for(unsigned r3:kernel) {
        std::array<unsigned,4> r{r0,r1,r2,r3};
        std::array<std::array<Polynomial,4>,16> local{};
        for(int x=0;x<2;x++)for(int y=0;y<2;y++)for(int u=0;u<2;u++)for(int vv=0;vv<2;vv++) {
          unsigned identity=x==u&&y==vv,endpoint=x==vv&&y==u;
          for(int a=0;a<2;a++)for(int b=0;b<2;b++) {
            Polynomial p=a==b?identity|((identity^endpoint)<<1):0;
            if(((r[2*a+b]>>(2*x+y))&1)&&((v>>(2*u+vv))&1))p^=6;
            local[8*x+4*y+2*u+vv][2*a+b]=p;
            assert(polynomial_parity(p)==unsigned(endpoint&&(a==b)));
          }
        }
        auto selected=optimize_local(local,2);
        auto spectrum=selected.spectrum;auto value=selected.value;checked++;
        if(value.bound<best-1e-9){best=value.bound;records.clear();}
        if(value.bound<best+1e-9) {
          std::string record="{\"v\":"+std::to_string(v)+",\"r\":["+std::to_string(r0)+","+std::to_string(r1)+","+std::to_string(r2)+","+std::to_string(r3)+"]"+gauge_physical_record(selected)+",\"infinity_height\":"+std::to_string(selected.infinity_height)+",\"zero_height\":"+std::to_string(selected.zero_height)+",\"s_numerator\":"+std::to_string(std::llround(value.s*1048576))+",\"s_denominator\":1048576,\"infinity_profile_index\":"+std::to_string(value.row)+",\"zero_profile_index\":"+std::to_string(value.column)+",\"degrees_at_infinity\":[";
          for(int i=0;i<4;i++)record+=(i?",":"")+std::to_string(spectrum[value.row][i]);
          record+="],\"orders_at_zero\":[";
          for(int i=0;i<4;i++)record+=(i?",":"")+std::to_string(spectrum[value.column][i]);
          record+="]}";records.push_back(record);
        }
      }
    }
    assert(checked==30720);
  } else if((mode=="physical-ket-column"||mode=="lattice-ket-column")||mode=="gauged-ket-column") {
    for(unsigned h=0;h<16;h++)for(unsigned g=0;g<16;g++)for(unsigned d=0;d<16;d++) {
      std::array<std::array<Polynomial,4>,16> local{};
      for(int x=0;x<2;x++)for(int y=0;y<2;y++)for(int u=0;u<2;u++)for(int v=0;v<2;v++) {
        int index=8*x+4*y+2*u+v;
        unsigned identity=x==u&&y==v,endpoint=x==v&&y==u;
        local[index][0]=identity|((identity^endpoint)<<1);
        local[index][1]=(y==v)&&((h>>(2*x+u))&1);
        local[index][2]=((y==v)&&((g>>(2*x+u))&1))?3:0;
        local[index][3]=(y==v)&&((d>>(2*x+u))&1);
        assert(polynomial_parity(local[index][0])==endpoint);
        assert(!polynomial_parity(local[index][2]));
      }
      auto selected=optimize_local(local,1);
      auto spectrum=selected.spectrum;auto value=selected.value;checked++;
      if(value.bound<best-1e-9){best=value.bound;records.clear();}
      if(value.bound<best+1e-9) {
        std::string record="{\"h\":"+std::to_string(h)+",\"g\":"+std::to_string(g)+",\"d\":"+std::to_string(d)+gauge_physical_record(selected)+",\"infinity_height\":"+std::to_string(selected.infinity_height)+",\"zero_height\":"+std::to_string(selected.zero_height)+",\"s_numerator\":"+std::to_string(std::llround(value.s*1048576))+",\"s_denominator\":1048576,\"infinity_profile_index\":"+std::to_string(value.row)+",\"zero_profile_index\":"+std::to_string(value.column)+",\"degrees_at_infinity\":[";
        for(int i=0;i<4;i++)record+=(i?",":"")+std::to_string(spectrum[value.row][i]);
        record+="],\"orders_at_zero\":[";
        for(int i=0;i<4;i++)record+=(i?",":"")+std::to_string(spectrum[value.column][i]);
        record+="]}";records.push_back(record);
      }
    }
  } else if(mode=="linear-factor"||mode=="gauged-linear-factor"||(mode=="physical-linear-factor"||mode=="lattice-linear-factor")) {
    for(unsigned A=0;A<65536;A++)for(int factor_rank=1;factor_rank<=2;factor_rank++) {
      std::array<std::array<Polynomial,4>,16> local{};
      for(int x=0;x<2;x++)for(int y=0;y<2;y++)for(int u=0;u<2;u++)for(int v=0;v<2;v++)
      for(int a=0;a<2;a++)for(int b=0;b<2;b++) {
        unsigned initial=((factor_rank==2)?y==v:!y&&!v)&&((A>>(4*(2*x+a)+(2*u+b)))&1);
        unsigned endpoint=x==v&&y==u&&a==b;
        local[8*x+4*y+2*u+v][2*a+b]=initial|((initial^endpoint)<<1);
      }
      auto selected=gauged?optimize_local(local,1):GaugedValue{};
      auto spectrum=gauged?selected.spectrum:profiles(local,1,0,ramification);
      auto value=gauged?selected.value:optimize(spectrum);checked++;
      if(value.bound<best-1e-9){best=value.bound;records.clear();}
      if(value.bound<best+1e-9) {
        std::string record="{\"A\":"+std::to_string(A)+",\"factor_rank\":"+std::to_string(factor_rank)+gauge_physical_record(selected)+",\"infinity_height\":"+std::to_string(selected.infinity_height)+",\"zero_height\":"+std::to_string(selected.zero_height)+",\"s_numerator\":"+std::to_string(std::llround(value.s*1048576))+",\"s_denominator\":1048576,\"infinity_profile_index\":"+std::to_string(value.row)+",\"zero_profile_index\":"+std::to_string(value.column)+",\"degrees_at_infinity\":[";
        for(int i=0;i<4;i++)record+=(i?",":"")+std::to_string(spectrum[value.row][i]);
        record+="],\"orders_at_zero\":[";
        for(int i=0;i<4;i++)record+=(i?",":"")+std::to_string(spectrum[value.column][i]);
        record+="]}";records.push_back(record);
      }
    }
  } else if(mode=="resolvents"||mode=="gauged-resolvents"||(mode=="physical-resolvents"||mode=="lattice-resolvents")) {
    int invertible=0;
    for(unsigned A=0;A<65536;A++) {
      if(binary_rank4(A)!=4)continue;
      invertible++;
      Binary8 initial{},nilpotent{},identity{};
      for(int i=0;i<8;i++)identity[i]=1u<<i;
      for(int x=0;x<2;x++)for(int y=0;y<2;y++)for(int u=0;u<2;u++)for(int v=0;v<2;v++)
      for(int a=0;a<2;a++)for(int b=0;b<2;b++) {
        int row=2*(2*x+y)+a,col=2*(2*u+v)+b;
        if(y==v&&((A>>(4*(2*x+a)+(2*u+b)))&1))initial[row]|=1u<<col;
        if(y==u&&((A>>(4*(2*x+a)+(2*v+b)))&1))nilpotent[row]|=1u<<col;
      }
      for(int i=0;i<8;i++)nilpotent[i]^=identity[i];
      Binary8 power=identity;
      int index=0;
      for(int k=1;k<=8;k++) {
        power=binary_product(power,nilpotent);
        if(power==Binary8{}){index=k;break;}
      }
      if(!index)continue;
      std::array<std::array<Polynomial,4>,16> local{};
      power=identity;
      for(int k=0;k<index;k++) {
        auto coefficient=binary_product(power,initial);
        for(int x=0;x<2;x++)for(int y=0;y<2;y++)for(int u=0;u<2;u++)for(int v=0;v<2;v++)
        for(int a=0;a<2;a++)for(int b=0;b<2;b++) {
          int row=2*(2*x+y)+a,col=2*(2*u+v)+b;
          if((coefficient[row]>>col)&1)local[8*x+4*y+2*u+v][2*a+b]|=Polynomial(1)<<k;
        }
        power=binary_product(power,nilpotent);
      }
      for(int x=0;x<2;x++)for(int y=0;y<2;y++)for(int u=0;u<2;u++)for(int v=0;v<2;v++)
      for(int a=0;a<2;a++)for(int b=0;b<2;b++)assert(polynomial_parity(local[8*x+4*y+2*u+v][2*a+b])==unsigned(x==v&&y==u&&a==b));
      int degree=index-1;
      auto selected=gauged?optimize_local(local,degree):GaugedValue{};
      auto spectrum=gauged?selected.spectrum:profiles(local,degree,0,ramification);
      auto value=gauged?selected.value:optimize(spectrum);checked++;
      if(value.bound<best-1e-9){best=value.bound;records.clear();}
      if(value.bound<best+1e-9) {
        std::string record="{\"A\":"+std::to_string(A)+",\"nilpotency_index\":"+std::to_string(index)+gauge_physical_record(selected)+",\"infinity_height\":"+std::to_string(selected.infinity_height)+",\"zero_height\":"+std::to_string(selected.zero_height)+",\"s_numerator\":"+std::to_string(std::llround(value.s*1048576))+",\"s_denominator\":1048576,\"infinity_profile_index\":"+std::to_string(value.row)+",\"zero_profile_index\":"+std::to_string(value.column)+",\"degrees_at_infinity\":[";
        for(int i=0;i<4;i++)record+=(i?",":"")+std::to_string(spectrum[value.row][i]);
        record+="],\"orders_at_zero\":[";
        for(int i=0;i<4;i++)record+=(i?",":"")+std::to_string(spectrum[value.column][i]);
        record+="]}";records.push_back(record);
      }
    }
    assert(invertible==20160);
  } else if(mode=="loops"||mode=="gauged-loops"||(mode=="physical-loops"||mode=="lattice-loops")) {
    for(unsigned first=1;first<256;first++)for(unsigned second=1;second<256;second++) {
      if(__builtin_parity(first&second))continue;
      std::array<std::array<Polynomial,4>,16> local{};
      for(int x=0;x<2;x++)for(int y=0;y<2;y++)for(int u=0;u<2;u++)for(int v=0;v<2;v++) {
        int identity=(x==u)&&(y==v),cross=(x==v)&&(y==u);
        for(int a=0;a<2;a++)for(int b=0;b<2;b++) {
          int row=2*(2*x+y)+a,col=2*(2*u+v)+b,swapcol=2*(2*v+u)+b;
          Polynomial p=(a==b)?unsigned(identity)|(unsigned(identity^cross)<<1):0;
          if((first>>row)&1) {
            if((second>>col)&1)p^=6;
            if(((second>>col)^(second>>swapcol))&1)p^=12;
          }
          local[8*x+4*y+2*u+v][2*a+b]=p;
          assert(polynomial_parity(p)==unsigned(cross&&(a==b)));
        }
      }
      auto selected=gauged?optimize_local(local,3):GaugedValue{};
      auto spectrum=gauged?selected.spectrum:profiles(local,3,0,ramification);
      auto value=gauged?selected.value:optimize(spectrum);checked++;
      if(value.bound<best-1e-9){best=value.bound;records.clear();}
      if(value.bound<best+1e-9) {
        std::string record="{\"u\":"+std::to_string(first)+",\"v\":"+std::to_string(second)+gauge_physical_record(selected)+",\"infinity_height\":"+std::to_string(selected.infinity_height)+",\"zero_height\":"+std::to_string(selected.zero_height)+",\"s_numerator\":"+std::to_string(std::llround(value.s*1048576))+",\"s_denominator\":1048576,\"infinity_profile_index\":"+std::to_string(value.row)+",\"zero_profile_index\":"+std::to_string(value.column)+",\"degrees_at_infinity\":[";
        for(int i=0;i<4;i++)record+=(i?",":"")+std::to_string(spectrum[value.row][i]);
        record+="],\"orders_at_zero\":[";
        for(int i=0;i<4;i++)record+=(i?",":"")+std::to_string(spectrum[value.column][i]);
        record+="]}";records.push_back(record);
      }
    }
  } else
  for(int grading=0;grading<2;grading++) {
    std::array<int,3> powers=grading?std::array<int,3>{0,1,2}:std::array<int,3>{0,0,1};
    int degree=grading?3:2;
    for(int br=0;br<16;br++)for(int bh=0;bh<16;bh++)for(int bc=0;bc<16;bc++) {
      std::array<int,3> bs{br,bh,bc};
      std::array<std::array<Polynomial,4>,16> local{};
      for(int x=0;x<2;x++)for(int y=0;y<2;y++)for(int u=0;u<2;u++)for(int v=0;v<2;v++) {
        int R=(x==y)&&(u==v),C=(x==u)&&(y==v),E=(x==v)&&(y==u),H=R^C^E;
        std::array<int,3> values{R,H,C};
        for(int s=0;s<3;s++)if(values[s])for(int a=0;a<2;a++)for(int b=0;b<2;b++) {
          Polynomial p=unsigned(a==b)^(((bs[s]>>(2*a+b))&1)?3:0);
          local[8*x+4*y+2*u+v][2*a+b]^=p<<powers[s];
        }
        for(int a=0;a<2;a++)for(int b=0;b<2;b++)assert(polynomial_parity(local[8*x+4*y+2*u+v][2*a+b])==unsigned(E&&(a==b)));
      }
      auto selected=gauged?optimize_local(local,degree):GaugedValue{};
      auto spectrum=gauged?selected.spectrum:profiles(local,degree,0,ramification);
      auto value=gauged?selected.value:optimize(spectrum);checked++;
      if(value.bound<best-1e-9){best=value.bound;records.clear();}
      if(value.bound<best+1e-9) {
        std::string record="{\"grading\":"+std::to_string(grading)+",\"matrices\":["+std::to_string(br)+","+std::to_string(bh)+","+std::to_string(bc)+"],\"infinity_height\":"+std::to_string(selected.infinity_height)+",\"zero_height\":"+std::to_string(selected.zero_height)+",\"s_numerator\":"+std::to_string(std::llround(value.s*1048576))+",\"s_denominator\":1048576,\"infinity_profile_index\":"+std::to_string(value.row)+",\"zero_profile_index\":"+std::to_string(value.column)+",\"degrees_at_infinity\":[";
        for(int i=0;i<4;i++)record+=(i?",":"")+std::to_string(spectrum[value.row][i]);
        record+="],\"orders_at_zero\":[";
        for(int i=0;i<4;i++)record+=(i?",":"")+std::to_string(spectrum[value.column][i]);
        record+="]}";records.push_back(record);
      }
    }
  }
  std::ofstream out(argv[1]);assert(out);
  out<<"{\"lattice_gauges\":"<<(lattice_mode?"true":"false")<<",\"mode\":\""<<mode<<"\",\"maximum_height\":"<<(gauged?maximum_height:0)<<",\"independent_auxiliary_bases\":"<<(independent_bases?"true":"false")<<",\"maximum_physical\":"<<maximum_physical<<",\"ramification\":"<<ramification<<",\"checked\":"<<checked<<",\"best_bound_estimate\":"<<std::setprecision(17)<<best<<",\"rho_estimate\":"<<std::log2(best)/4<<",\"minimizers\":[";
  for(size_t i=0;i<records.size();i++)out<<(i?",":"")<<records[i];
  out<<"]}\n";
  std::cout<<"models="<<checked<<" best="<<std::setprecision(17)<<best<<" rho="<<std::log2(best)/4<<" minimizers="<<records.size()<<"\n";
  return 0;
}
