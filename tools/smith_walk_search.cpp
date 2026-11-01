// Degree-sensitive auxiliary walk bounds for binary polynomial crossing MPOs.
// All algebra is exact over F_2; optimization is only a numerical screen.
#define main plain_smith_search_main
#include "smith_mpo_search.cpp"
#undef main
#include <functional>
#include <random>
#include <set>
#include <sstream>

using Local=std::array<std::array<Polynomial,4>,16>;
using Rectangular=std::vector<std::array<Polynomial,4>>;
static constexpr int infinity_order=1000000;
static bool tropical_enabled=false;
struct WalkFactor {
  Spectrum base{};
  std::array<std::array<int,4>,4> kappa{};
  std::array<std::array<int,2>,4> pivots{};
  int precision=0;
  std::vector<std::array<int,4>> tropical_next,tropical_exponent;
};

static void determinize(WalkFactor& factor) {
  std::array<std::array<int,4>,4> costs=factor.kappa;
  int maximum=0;
  for(int i=0;i<4;i++)for(int ab=0;ab<4;ab++) {
    // Replacing infinity by any finite lower bound only enlarges the
    // counted language. The precision bound is convenient and fixed.
    costs[i][ab]=std::min(costs[i][ab],factor.precision);
    maximum=std::max(maximum,costs[i][ab]);
  }
  std::map<int,int> index{{0,0}};std::vector<int> states{0};
  for(size_t state=0;state<states.size();state++) {
    int difference=states[state],a=std::max(0,-difference),b=std::max(0,difference);
    std::array<int,4> next{},exponent{};
    for(int channel=0;channel<4;channel++) {
      int x=std::min(a+costs[channel][0],b+costs[channel][2]);
      int y=std::min(a+costs[channel][1],b+costs[channel][3]);
      int removed=std::min(x,y),difference_next=y-x;
      auto found=index.find(difference_next);
      if(found==index.end()) {
        int id=states.size();index[difference_next]=id;states.push_back(difference_next);next[channel]=id;
      } else next[channel]=found->second;
      exponent[channel]=factor.base[channel]-removed;
    }
    factor.tropical_next.push_back(next);factor.tropical_exponent.push_back(exponent);
  }
  assert(states.size()<=size_t(2*maximum+1));
}

static WalkFactor schur_walk(Rectangular residual,int degree,int shift,int pivot_rule=0) {
  // The residual after k-1 pivots is a bordered k-minor divided by
  // a (k-1)-minor. Thus precision 4*degree+1 also detects exact zeros.
  int precision=4*degree+1;
  assert(precision<128);
  WalkFactor result;result.precision=precision;
  const Polynomial mask=(Polynomial(1)<<precision)-1;
  for(auto& row:residual)for(auto& value:row)value&=mask;
  for(int channel=0;channel<4;channel++) {
    int pivot_row=-1,pivot_column=-1,order=precision;
    for(int ii=0;ii<16;ii++)for(int jj=0;jj<4;jj++) {
      int i=(pivot_rule&1)?15-ii:ii;
      i=(i+2*(pivot_rule/4))%16;
      int j=(pivot_rule&2)?3-jj:jj;
      if(!residual[i][j])continue;
      int valuation=trailing_zeros(residual[i][j]);
      if(valuation<order){order=valuation;pivot_row=i;pivot_column=j;}
    }
    assert(pivot_row>=0);
    result.base[channel]=shift-order;
    result.pivots[channel]={pivot_row,pivot_column};
    result.kappa[channel].fill(infinity_order);
    std::array<Polynomial,16> left{};
    std::array<Polynomial,4> right{};
    int available=precision-order;
    Polynomial inverse=inverse_unit(residual[pivot_row][pivot_column]>>order,available);
    for(int i=0;i<16;i++) {
      left[i]=residual[i][pivot_column]>>order;
      if(left[i])result.kappa[channel][i%4]=std::min(result.kappa[channel][i%4],trailing_zeros(left[i]));
    }
    for(int j=0;j<4;j++)right[j]=multiply(residual[pivot_row][j]>>order,inverse,available);
    for(int i=0;i<16;i++)for(int j=0;j<4;j++)
      residual[i][j]^=multiply(left[i],right[j],available)<<order;
    for(int j=0;j<4;j++)assert(!residual[pivot_row][j]);
    for(int i=0;i<16;i++)assert(!residual[i][pivot_column]);
  }
  for(auto row:residual)for(auto value:row)assert(!value);
  if(tropical_enabled)determinize(result);
  return result;
}

static std::array<WalkFactor,8> local_walks(const Local& local,int original_degree,
    std::array<int,3> lattice,int ramification,std::array<int,4> physical,int pivot_rule=0) {
  int a=lattice[0],b=lattice[1],h=a+b,degree=original_degree*ramification;
  Polynomial p=lattice[2],p_squared=multiply(p,p,127);
  int offset=0;for(int value:physical)offset+=std::abs(value);
  int bound_degree=degree+2*h+2*offset;
  assert(4*bound_degree+1<128);
  std::array<std::array<std::array<Polynomial,4>,16>,2> changed{};
  for(int labels=0;labels<16;labels++) {
    std::array<Polynomial,4> raw{};
    for(int ab=0;ab<4;ab++)for(int k=0;k<=original_degree;k++)
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
  std::array<WalkFactor,8> out{};
  for(int cut=0;cut<2;cut++)for(int side=0;side<2;side++)for(int infinity=0;infinity<2;infinity++) {
    Rectangular matrix(16);
    for(int x=0;x<2;x++)for(int y=0;y<2;y++)for(int u=0;u<2;u++)for(int v=0;v<2;v++) {
      int left=cut?2*x+u:2*x+y,right=cut?2*y+v:2*u+v;
      int shift=offset+(2*x-1)*physical[0]+(2*y-1)*physical[1]+(2*u-1)*physical[2]+(2*v-1)*physical[3];
      for(int ab=0;ab<4;ab++)matrix[4*(side?right:left)+ab][side?left:right]=changed[infinity][8*x+4*y+2*u+v][ab]<<shift;
    }
    int restore=h+offset+(infinity?degree:0);
    out[4*cut+2*side+infinity]=schur_walk(matrix,bound_degree,restore,pivot_rule);
  }
  return out;
}

static double plain_cost(const WalkFactor& f,double s) {
  double out=0;for(int exponent:f.base)out+=std::pow(s,exponent);return out;
}
static std::array<double,4> walk_matrix(const WalkFactor& f,double s) {
  std::array<double,4> out{};
  for(int i=0;i<4;i++)for(int ab=0;ab<4;ab++)if(f.kappa[i][ab]<infinity_order)
    out[ab]+=std::pow(s,f.base[i]-f.kappa[i][ab]);
  return out;
}
static double walk_cost(const WalkFactor& f,double s) {
  auto w=walk_matrix(f,s);
  return (w[0]+w[3]+std::sqrt((w[0]-w[3])*(w[0]-w[3])+4*w[1]*w[2]))/2;
}
static double tropical_cost(const WalkFactor& f,double s) {
  int n=f.tropical_next.size();assert(n);
  std::vector<std::array<double,4>> weights(n);
  double lazy=0,minimum_row=1e100;
  for(int i=0;i<n;i++) {
    double row=0;
    for(int j=0;j<4;j++){weights[i][j]=std::pow(s,f.tropical_exponent[i][j]);row+=weights[i][j];}
    lazy=std::max(lazy,row);
    minimum_row=std::min(minimum_row,row);
  }
  if(lazy-minimum_row<1e-14*lazy)return lazy;
  std::vector<double> potential(n,1),changed(n);
  double upper=lazy;
  int stable=0;
  for(int iteration=0;iteration<128;iteration++) {
    double estimate=0,maximum=0;
    for(int i=0;i<n;i++) {
      double value=0;for(int j=0;j<4;j++)value+=weights[i][j]*potential[f.tropical_next[i][j]];
      estimate=std::max(estimate,value/potential[i]);
      changed[i]=value+lazy*potential[i];maximum=std::max(maximum,changed[i]);
    }
    double previous=upper;upper=std::min(upper,estimate);
    stable=previous-upper<1e-14*upper?stable+1:0;
    if(stable>=8)break;
    for(int i=0;i<n;i++)potential[i]=changed[i]/maximum;
  }
  return upper;
}
struct WalkChoice {
  WalkFactor factor{};
  std::array<int,3> lattice{};
  std::array<int,4> physical{};
  int index=0,basis=9,pivot_rule=0;
  bool walk=false;
  double cost=1e100;
  int method=0;
};
struct WalkValue {double bound=1e100,s=1;WalkChoice infinity{},zero{};};
static double selected_cost(const WalkChoice& choice,double s) {
  return choice.method==2?tropical_cost(choice.factor,s):(choice.walk?walk_cost(choice.factor,s):plain_cost(choice.factor,s));
}
static std::vector<std::array<int,3>> hermite_forms(int maximum_depth) {
  std::vector<std::array<int,3>> out;
  for(int depth=0;depth<=maximum_depth;depth++)for(int a=0;a<=depth;a++) {
    int b=depth-a;
    for(int p=0;p<(1<<b);p++)if(!a||!b||(p&1))out.push_back({a,b,p});
  }
  return out;
}
static const std::vector<double> search_parameters={1,1.0625,1.125,1.25,1.5,1.75,2,2.5,3,4,6,8,12,16};
static WalkValue optimize_walks(const Local& local,int degree,int depth,int ramification,
    int physical_radius,int basis_count,int pivot_rules) {
  const int size=search_parameters.size();
  std::vector<std::array<WalkChoice,4>> best(size);
  std::vector<int> bases=basis_count==6?std::vector<int>{6,7,9,11,13,14}:std::vector<int>{9};
  auto lattices=hermite_forms(depth);
  int r=physical_radius;
  for(int basis:bases) {
    Local conjugated=auxiliary_conjugate(local,basis);
    for(auto lattice:lattices)
    for(int hx=-r;hx<=r;hx++)for(int hy=-r;hy<=r;hy++)for(int hu=-r;hu<=r;hu++)for(int hv=-r;hv<=r;hv++) {
      std::array<int,4> physical{hx,hy,hu,hv};
      for(int pivot_rule=0;pivot_rule<pivot_rules;pivot_rule++) {
        auto factors=local_walks(conjugated,degree,lattice,ramification,physical,pivot_rule);
        for(int i=0;i<size;i++)for(int index=0;index<8;index++) {
          int slot=2*(index/4)+(index%2);
          double plain=plain_cost(factors[index],search_parameters[i]);
          double walk=walk_cost(factors[index],search_parameters[i]);
          int method=walk<plain?1:0;
          double cost=std::min(plain,walk);
          if(tropical_enabled) {
            double deterministic=tropical_cost(factors[index],search_parameters[i]);
            if(deterministic<cost){method=2;cost=deterministic;}
          }
          if(cost<best[i][slot].cost-1e-12)best[i][slot]={factors[index],lattice,physical,index,basis,pivot_rule,method!=0,cost,method};
        }
      }
    }
  }
  WalkValue result;
  std::set<std::vector<int>> refined;
  for(int i=0;i<size;i++)for(int orientation=0;orientation<2;orientation++) {
    auto first=best[i][2*orientation+1],second=best[i][2*(1-orientation)];
    double bound=first.cost*second.cost;
    if(bound<result.bound)result={bound,search_parameters[i],first,second};
    std::vector<int> signature;
    for(auto candidate:{first,second}) {
      signature.push_back(candidate.method);
      for(int value:candidate.factor.base)signature.push_back(value);
      if(candidate.walk)for(auto row:candidate.factor.kappa)for(int value:row)signature.push_back(value);
    }
    if(!refined.insert(signature).second)continue;
    // Refine each selected pair. The log of either fixed path or plain
    // cost is convex in log(s). Rounding gives a reproducible rational s.
    double low=0,high=std::log(64.0);
    for(int iteration=0;iteration<64;iteration++) {
      double a=(2*low+high)/3,b=(low+2*high)/3;
      double va=selected_cost(first,std::exp(a))*selected_cost(second,std::exp(a));
      double vb=selected_cost(first,std::exp(b))*selected_cost(second,std::exp(b));
      if(va<=vb)high=b;else low=a;
    }
    double s=std::round(std::exp((low+high)/2)*1048576)/1048576.;
    bound=selected_cost(first,s)*selected_cost(second,s);
    if(bound<result.bound)result={bound,s,first,second};
  }
  return result;
}
static void write_choice(std::ostream& out,const WalkChoice& c) {
  out<<"{\"profile_index\":"<<c.index<<",\"basis\":"<<c.basis<<",\"pivot_rule\":"<<c.pivot_rule<<",\"walk\":"<<(c.walk?"true":"false")<<",\"method\":"<<c.method<<",\"tropical_states\":"<<c.factor.tropical_next.size()<<",\"precision\":"<<c.factor.precision<<",\"lattice\":[";
  for(int i=0;i<3;i++)out<<(i?",":"")<<c.lattice[i];
  out<<"],\"physical\":[";for(int i=0;i<4;i++)out<<(i?",":"")<<c.physical[i];
  out<<"],\"base_exponents\":[";for(int i=0;i<4;i++)out<<(i?",":"")<<c.factor.base[i];
  out<<"],\"kappa\":[";
  for(int i=0;i<4;i++){out<<(i?",[":"[");for(int ab=0;ab<4;ab++){if(ab)out<<",";if(c.factor.kappa[i][ab]==infinity_order)out<<"null";else out<<c.factor.kappa[i][ab];}out<<"]";}
  out<<"],\"pivots\":[";for(int i=0;i<4;i++)out<<(i?",[":"[")<<c.factor.pivots[i][0]<<","<<c.factor.pivots[i][1]<<"]";out<<"]}";
}
static Local elementary_local() {
  Local out{};
  for(int x=0;x<2;x++)for(int y=0;y<2;y++)for(int u=0;u<2;u++)for(int v=0;v<2;v++) {
    unsigned I=x==u&&y==v,E=x==v&&y==u;
    out[8*x+4*y+2*u+v][0]=out[8*x+4*y+2*u+v][3]=I|((I^E)<<1);
  }
  return out;
}
static Local trace_shear_local() {
  Local out{};
  for(int x=0;x<2;x++)for(int y=0;y<2;y++)for(int u=0;u<2;u++)for(int v=0;v<2;v++) {
    unsigned R=x==y&&u==v,C=x==u&&y==v,E=x==v&&y==u;
    out[8*x+4*y+2*u+v]={R,(R^E)<<1,C^E,C<<1};
  }
  return out;
}
static Local s3_local(unsigned a0,unsigned a1) {
  constexpr std::array<std::array<int,3>,6> permutations={{{0,1,2},{0,2,1},{1,0,2},{1,2,0},{2,0,1},{2,1,0}}};
  Local out{};unsigned a2=4^a0^a1;
  for(int p=0;p<6;p++) {
    Polynomial polynomial=((a0>>p)&1)|(((a1>>p)&1)<<1)|(((a2>>p)&1)<<2);
    for(int u=0;u<2;u++)for(int v=0;v<2;v++)for(int b=0;b<2;b++) {
      std::array<int,3> inputs{u,v,b};
      int x=inputs[permutations[p][0]],y=inputs[permutations[p][1]],a=inputs[permutations[p][2]];
      out[8*x+4*y+2*u+v][2*a+b]^=polynomial;
    }
  }
  return out;
}

#ifndef SMITH_WALK_NO_MAIN
int main(int argc,char** argv) {
  assert(argc>=3&&argc<=11);
  std::string output=argv[1],mode=argv[2];
  if(mode.rfind("det-",0)==0){tropical_enabled=true;mode=mode.substr(4);}
  int depth=argc>3?std::stoi(argv[3]):0,ramification=argc>4?std::stoi(argv[4]):1;
  int physical_radius=argc>5?std::stoi(argv[5]):0,basis_count=argc>6?std::stoi(argv[6]):1;
  int pivot_rules=argc>7?std::stoi(argv[7]):1;
  WalkValue best;std::string best_model;Local best_local{};int best_degree=0,count=0;
  std::map<std::string,int> histogram;
  auto submit=[&](const Local& local,int degree,const std::string& model) {
    auto value=optimize_walks(local,degree,depth,ramification,physical_radius,basis_count,pivot_rules);count++;
    std::ostringstream key;key<<std::setprecision(8)<<value.bound;histogram[key.str()]++;
    if(value.bound<best.bound-1e-10) {
      best=value;best_model=model;best_local=local;best_degree=degree;
      std::cout<<"newbest models="<<count<<" B="<<std::setprecision(17)<<best.bound<<" rho="<<std::log2(best.bound)/4<<" model="<<model<<std::endl;
    }
    if(count%256==0)std::cout<<"models="<<count<<std::endl;
  };
  if(mode=="elementary")submit(elementary_local(),1,"{\"kind\":\"elementary\"}");
  else if(mode=="trace-shear")submit(trace_shear_local(),1,"{\"kind\":\"trace-shear\"}");
  else if(mode=="s3")for(unsigned a0=0;a0<64;a0++)for(unsigned a1=0;a1<64;a1++)
    submit(s3_local(a0,a1),2,"{\"layers\":["+std::to_string(a0)+","+std::to_string(a1)+"]}");
  else if(mode=="conjugations") {
    for(unsigned first=1;first<256;first++)for(unsigned second=1;second<256;second++) {
      if(__builtin_parity(first&second))continue;
      unsigned swapped_first=0;
      for(int x=0;x<2;x++)for(int y=0;y<2;y++)for(int a=0;a<2;a++)
        if((first>>(2*(2*y+x)+a))&1)swapped_first|=1u<<(2*(2*x+y)+a);
      unsigned coupling=__builtin_parity(second&swapped_first);
      Local local=elementary_local();
      for(int x=0;x<2;x++)for(int y=0;y<2;y++)for(int u=0;u<2;u++)for(int v=0;v<2;v++)
      for(int a=0;a<2;a++)for(int b=0;b<2;b++) {
        int row=2*(2*x+y)+a,col=2*(2*u+v)+b,swaprow=2*(2*y+x)+a,swapcol=2*(2*v+u)+b;
        unsigned mixed=(((first>>row)&1)&((second>>swapcol)&1))^(((first>>swaprow)&1)&((second>>col)&1));
        if(mixed)local[8*x+4*y+2*u+v][2*a+b]^=6;
        if(coupling&&((first>>row)&1)&&((second>>col)&1))local[8*x+4*y+2*u+v][2*a+b]^=10;
      }
      submit(local,coupling?3:2,"{\"u\":"+std::to_string(first)+",\"v\":"+std::to_string(second)+",\"coupling\":"+std::to_string(coupling)+"}");
    }
  }
  else if(mode=="loops") {
    for(unsigned first=1;first<256;first++)for(unsigned second=1;second<256;second++) {
      if(__builtin_parity(first&second))continue;
      Local local=elementary_local();
      for(int x=0;x<2;x++)for(int y=0;y<2;y++)for(int u=0;u<2;u++)for(int v=0;v<2;v++)
      for(int a=0;a<2;a++)for(int b=0;b<2;b++) {
        int row=2*(2*x+y)+a,col=2*(2*u+v)+b,swapcol=2*(2*v+u)+b;
        if((first>>row)&1) {
          if((second>>col)&1)local[8*x+4*y+2*u+v][2*a+b]^=6;
          if(((second>>col)^(second>>swapcol))&1)local[8*x+4*y+2*u+v][2*a+b]^=12;
        }
      }
      submit(local,3,"{\"u\":"+std::to_string(first)+",\"v\":"+std::to_string(second)+"}");
    }
  } else if(mode=="bridge") {
    for(unsigned v=1;v<16;v++) {
      unsigned swapped=(v&9)|((v&2)<<1)|((v&4)>>1);
      std::vector<unsigned> kernel;
      for(unsigned r=0;r<16;r++)if(!__builtin_parity(r&v)&&!__builtin_parity(r&swapped))kernel.push_back(r);
      for(unsigned r0:kernel)for(unsigned r1:kernel)for(unsigned r2:kernel)for(unsigned r3:kernel) {
        std::array<unsigned,4> r{r0,r1,r2,r3};Local local=elementary_local();
        for(int x=0;x<2;x++)for(int y=0;y<2;y++)for(int u=0;u<2;u++)for(int vv=0;vv<2;vv++)
        for(int ab=0;ab<4;ab++)if(((r[ab]>>(2*x+y))&1)&&((v>>(2*u+vv))&1))local[8*x+4*y+2*u+vv][ab]^=6;
        submit(local,2,"{\"v\":"+std::to_string(v)+",\"r\":["+std::to_string(r0)+","+std::to_string(r1)+","+std::to_string(r2)+","+std::to_string(r3)+"]}");
      }
    }
  } else if(mode=="resolvents") {
    for(unsigned A=0;A<65536;A++) {
      if(binary_rank4(A)!=4)continue;
      Binary8 initial{},nilpotent{},identity{};
      for(int i=0;i<8;i++)identity[i]=1u<<i;
      for(int x=0;x<2;x++)for(int y=0;y<2;y++)for(int u=0;u<2;u++)for(int v=0;v<2;v++)
      for(int a=0;a<2;a++)for(int b=0;b<2;b++) {
        int row=2*(2*x+y)+a,col=2*(2*u+v)+b;
        if(y==v&&((A>>(4*(2*x+a)+(2*u+b)))&1))initial[row]|=1u<<col;
        if(y==u&&((A>>(4*(2*x+a)+(2*v+b)))&1))nilpotent[row]|=1u<<col;
      }
      for(int i=0;i<8;i++)nilpotent[i]^=identity[i];
      Binary8 power=identity;int index=0;
      for(int k=1;k<=8;k++){power=binary_product(power,nilpotent);if(power==Binary8{}){index=k;break;}}
      if(!index)continue;
      Local local{};power=identity;
      for(int k=0;k<index;k++) {
        auto coefficient=binary_product(power,initial);
        for(int x=0;x<2;x++)for(int y=0;y<2;y++)for(int u=0;u<2;u++)for(int v=0;v<2;v++)
        for(int a=0;a<2;a++)for(int b=0;b<2;b++) {
          int row=2*(2*x+y)+a,col=2*(2*u+v)+b;
          if((coefficient[row]>>col)&1)local[8*x+4*y+2*u+v][2*a+b]|=Polynomial(1)<<k;
        }
        power=binary_product(power,nilpotent);
      }
      submit(local,index-1,"{\"A\":"+std::to_string(A)+",\"nilpotency_index\":"+std::to_string(index)+"}");
    }
  }
  else if(mode=="file") {
    assert(argc==10);int degree=std::stoi(argv[8]);std::ifstream input(argv[9]);assert(input);Local local{};
    for(auto& entries:local)for(auto& value:entries){uint64_t raw;input>>raw;assert(input);value=raw;}
    submit(local,degree,"{\"kind\":\"file\",\"filename\":\""+std::string(argv[9])+"\"}");
  } else assert(false);
  std::ofstream out(output);assert(out);
  out<<"{\"mode\":\""<<mode<<"\",\"depth\":"<<depth<<",\"ramification\":"<<ramification<<",\"physical_radius\":"<<physical_radius<<",\"basis_count\":"<<basis_count<<",\"pivot_rules\":"<<pivot_rules<<",\"models\":"<<count<<",\"best_bound_estimate\":"<<std::setprecision(17)<<best.bound<<",\"rho_estimate\":"<<std::log2(best.bound)/4<<",\"model\":"<<best_model<<",\"degree\":"<<best_degree<<",\"s_numerator\":"<<std::llround(best.s*1048576)<<",\"s_denominator\":1048576,\"infinity\":";
  write_choice(out,best.infinity);out<<",\"zero\":";write_choice(out,best.zero);
  out<<",\"local\":[";
  for(int labels=0;labels<16;labels++){out<<(labels?",[":"[");for(int ab=0;ab<4;ab++)out<<(ab?",":"")<<static_cast<uint64_t>(best_local[labels][ab]);out<<"]";}
  out<<"],\"bound_histogram\":{";bool comma=false;for(auto entry:histogram){out<<(comma?",":"")<<"\""<<entry.first<<"\":"<<entry.second;comma=true;}out<<"}}\n";
  std::cout<<"complete models="<<count<<" best="<<std::setprecision(17)<<best.bound<<" rho="<<std::log2(best.bound)/4<<std::endl;
}
#endif
