// Optional GF(4) weighted Smith search engine. Not an exact real-weight certificate.
// Exact field arithmetic uses omega^2=omega+1 and 512 coefficients per series.
// A missing pivot returns an error; no incomplete profile is reported.
// The Python wrapper builds this source in a temporary cache on first use.
// Exploratory weighted Smith oracle with a 512-coefficient series window.
// Polynomial arithmetic is exact; floating pivot choices require exact replay.
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <vector>
#if defined(__PCLMUL__)
#include <wmmintrin.h>
#endif

constexpr int words=8,precision=64*words;
struct Poly{
 std::array<std::uint64_t,words>w{};
 Poly(std::uint64_t x=0){w[0]=x;}
 explicit operator bool()const{for(auto x:w)if(x)return true;return false;}
 Poly& operator^=(const Poly&b){for(int i=0;i<words;++i)w[i]^=b.w[i];return *this;}
 Poly& operator&=(const Poly&b){for(int i=0;i<words;++i)w[i]&=b.w[i];return *this;}
 friend Poly operator^(Poly a,const Poly&b){return a^=b;}
 friend Poly operator&(Poly a,const Poly&b){return a&=b;}
 friend bool operator==(const Poly&a,const Poly&b){return a.w==b.w;}
 friend bool operator!=(const Poly&a,const Poly&b){return !(a==b);}
 Poly operator<<(int s)const{
  Poly r;if(s<0)throw std::runtime_error("negative shift");if(s>=precision)return r;
  int k=s/64,b=s%64;for(int i=words-1;i>=k;--i){r.w[i]=w[i-k]<<b;if(b&&i>k)r.w[i]|=w[i-k-1]>>(64-b);}return r;
 }
 Poly operator>>(int s)const{
  Poly r;if(s<0)throw std::runtime_error("negative shift");if(s>=precision)return r;
  int k=s/64,b=s%64;for(int i=0;i+k<words;++i){r.w[i]=w[i+k]>>b;if(b&&i+k+1<words)r.w[i]|=w[i+k+1]<<(64-b);}return r;
 }
};
int valuation(const Poly&a){for(int i=0;i<words;++i)if(a.w[i])return 64*i+__builtin_ctzll(a.w[i]);throw std::runtime_error("zero valuation");}
int degree(const Poly&a){for(int i=words-1;i>=0;--i)if(a.w[i])return 64*i+63-__builtin_clzll(a.w[i]);return -1;}
Poly low_mask(int bits){Poly p;for(int i=0;i<words;++i){int b=std::clamp(bits-64*i,0,64);p.w[i]=b==64?~std::uint64_t(0):(b?(std::uint64_t(1)<<b)-1:0);}return p;}
Poly multiply(const Poly&a,const Poly&b){
 Poly out;
 for(int i=0;i<words;++i)if(a.w[i])for(int j=0;i+j<words;++j)if(b.w[j]){
#if defined(__PCLMUL__)
  auto x=_mm_set_epi64x(0,a.w[i]),y=_mm_set_epi64x(0,b.w[j]);auto z=_mm_clmulepi64_si128(x,y,0);
  alignas(16) std::uint64_t v[2];_mm_store_si128((__m128i*)v,z);out.w[i+j]^=v[0];if(i+j+1<words)out.w[i+j+1]^=v[1];
#else
  std::uint64_t x=a.w[i];
  while(x){int bit=__builtin_ctzll(x);x&=x-1;out.w[i+j]^=b.w[j]<<bit;
   if(bit&&i+j+1<words)out.w[i+j+1]^=b.w[j]>>(64-bit);}
#endif
 }
 return out;
}
namespace f4 {
struct Poly {
 ::Poly lo,hi;
 Poly(unsigned coefficient=0):lo(coefficient&1),hi((coefficient>>1)&1){}
 explicit operator bool()const{return bool(lo)||bool(hi);}
 Poly& operator^=(const Poly&b){lo^=b.lo;hi^=b.hi;return *this;}
 Poly& operator&=(const ::Poly&b){lo&=b;hi&=b;return *this;}
 friend Poly operator^(Poly a,const Poly&b){return a^=b;}
 friend bool operator==(const Poly&a,const Poly&b){return a.lo==b.lo&&a.hi==b.hi;}
 friend bool operator!=(const Poly&a,const Poly&b){return !(a==b);}
 Poly operator<<(int s)const{Poly r;r.lo=lo<<s;r.hi=hi<<s;return r;}
 Poly operator>>(int s)const{Poly r;r.lo=lo>>s;r.hi=hi>>s;return r;}
 unsigned constant()const{return unsigned(lo.w[0]&1)|unsigned((hi.w[0]&1)<<1);}
};
int valuation(const Poly&a){return std::min(a.lo?::valuation(a.lo):precision,a.hi?::valuation(a.hi):precision);}
int degree(const Poly&a){return std::max(::degree(a.lo),::degree(a.hi));}
Poly multiply(const Poly&a,const Poly&b){
 if(!a.hi&&!b.hi){Poly out;out.lo=::multiply(a.lo,b.lo);return out;}
 auto p=::multiply(a.lo,b.lo),q=::multiply(a.hi,b.hi),r=::multiply(a.lo^a.hi,b.lo^b.hi);
 Poly out;out.lo=p^q;out.hi=r^p;return out;
}
struct Series{Poly bits;int base=0;};
struct Pivot{double value;int base;std::array<int,4>labels;};
using Matrix=std::vector<std::vector<Series>>;
int alphabet=0;
std::array<Matrix,2>matrix_by_cut;
std::array<std::vector<std::array<int,4>>,2>row_labels,column_labels;
struct Block{Matrix matrix;std::vector<std::array<int,4>>rows,cols;};
std::array<std::vector<Block>,2>blocks_by_cut;

Series normalize(Series a,double tag,double scale,double cap){
 if(!a.bits)return {};
 int keep=int(std::ceil((cap-tag)/scale-a.base-1e-8));if(keep<=0)return {};
 if(keep<precision)a.bits&=::low_mask(keep);
 if(a.bits){int v=valuation(a.bits);a.base+=v;a.bits=a.bits>>v;}else a.base=0;
 return a;
}
Series add(Series a,Series b,double tag,double scale,double cap){
 if(!a.bits)return normalize(b,tag,scale,cap);if(!b.bits)return normalize(a,tag,scale,cap);
 if(a.base>b.base)std::swap(a,b);a.bits^=b.bits<<(b.base-a.base);return normalize(a,tag,scale,cap);
}
double height(const double*g,const std::array<int,4>&labels){double h=0;for(int l=0;l<4;++l)if(labels[l])h+=g[(alphabet-1)*l+labels[l]-1];return h;}
std::vector<Pivot> profile_block(const double*g,const Matrix&matrix,
                               std::vector<std::array<int,4>>rows,
                               std::vector<std::array<int,4>>cols,bool infinity,double scale){
 auto a=matrix;int n=a.size();
 std::vector<double>r(n),c(n);for(int i=0;i<n;++i){r[i]=height(g,rows[i]);c[i]=height(g,cols[i]);}
 if(infinity){for(auto&row:a)for(auto&entry:row)if(entry.bits){
  int d=degree(entry.bits);Poly rev;
  for(int k=0;k<=d;++k){
   if(entry.bits.lo.w[k/64]>>(k%64)&1)rev.lo.w[(d-k)/64]|=std::uint64_t(1)<<((d-k)%64);
   if(entry.bits.hi.w[k/64]>>(k%64)&1)rev.hi.w[(d-k)/64]|=std::uint64_t(1)<<((d-k)%64);
  }
  entry.bits=rev;entry.base=-d;
 }}else{for(auto&row:a)for(auto&entry:row)if(entry.bits){int v=valuation(entry.bits);entry.bits=entry.bits>>v;entry.base=v;}}
 double minimum=1e100;for(int i=0;i<n;++i)for(int j=0;j<n;++j)if(a[i][j].bits)minimum=std::min(minimum,scale*a[i][j].base+r[i]+c[j]);
 double cap=minimum+precision*scale;
 for(int i=0;i<n;++i)for(int j=0;j<n;++j)a[i][j]=normalize(a[i][j],r[i]+c[j],scale,cap);
 std::vector<Pivot>out;
 for(int k=0;k<n;++k){
  double best=1e100;int pi=-1,pj=-1;
  for(int i=k;i<n;++i)for(int j=k;j<n;++j)if(a[i][j].bits){double x=scale*a[i][j].base+r[i]+c[j];if(x<best){best=x;pi=i;pj=j;}}
  if(pi<0)return {};
  std::swap(a[k],a[pi]);std::swap(r[k],r[pi]);std::swap(rows[k],rows[pi]);
  for(auto&row:a)std::swap(row[k],row[pj]);std::swap(c[k],c[pj]);std::swap(cols[k],cols[pj]);
  std::array<int,4>labels;for(int l=0;l<4;++l)labels[l]=rows[k][l]+cols[k][l];
  int pivotbase=a[k][k].base;out.push_back({best,pivotbase,labels});
  Poly unit=a[k][k].bits;unsigned inv_table[4]={0,1,3,2};Poly inv(inv_table[unit.constant()]);
  if(!inv)throw std::runtime_error("zero leading coefficient");
  if(unit!=Poly(1))for(int p=1;p<precision;p*=2)inv=multiply(unit,multiply(inv,inv));
  for(int j=k;j<n;++j)if(a[k][j].bits){a[k][j].bits=multiply(a[k][j].bits,inv);a[k][j]=normalize(a[k][j],r[k]+c[j],scale,cap);}
  if(a[k][k].bits!=Poly(1)||a[k][k].base!=pivotbase)throw std::runtime_error("pivot");
  std::vector<int>support;for(int j=k+1;j<n;++j)if(a[k][j].bits)support.push_back(j);
  for(int i=k+1;i<n;++i)if(a[i][k].bits){auto factor=a[i][k];factor.base-=pivotbase;
   for(int j:support){Series p{multiply(factor.bits,a[k][j].bits),factor.base+a[k][j].base};a[i][j]=add(a[i][j],p,r[i]+c[j],scale,cap);}a[i][k]={};
  }
 }
 return out;
}
std::vector<Pivot> profile(const double*g,int cut,bool infinity,double scale){
 if(blocks_by_cut[cut].empty())return profile_block(g,matrix_by_cut[cut],row_labels[cut],column_labels[cut],infinity,scale);
 std::vector<Pivot>out;
 for(const auto&block:blocks_by_cut[cut]){
  auto part=profile_block(g,block.matrix,block.rows,block.cols,infinity,scale);
  if(part.size()!=block.matrix.size())return {};
  out.insert(out.end(),part.begin(),part.end());
 }
 std::stable_sort(out.begin(),out.end(),[](const auto&a,const auto&b){return a.value<b.value;});
 return out;
}
void split_support_components(){
 int n=alphabet*alphabet;
 for(int cut=0;cut<2;++cut){
  std::vector<std::vector<int>>adj(2*n);std::vector<bool>seen(2*n,false);
  for(int i=0;i<n;++i)for(int j=0;j<n;++j)if(matrix_by_cut[cut][i][j].bits){adj[i].push_back(n+j);adj[n+j].push_back(i);}
  std::vector<std::pair<std::vector<int>,std::vector<int>>>components;
  bool balanced=true;
  for(int root=0;root<n;++root)if(!seen[root]){
   std::vector<int>todo{root},rr,cc;seen[root]=true;
   while(!todo.empty()){int a=todo.back();todo.pop_back();if(a<n)rr.push_back(a);else cc.push_back(a-n);
    for(int b:adj[a])if(!seen[b]){seen[b]=true;todo.push_back(b);}
   }
   std::sort(rr.begin(),rr.end());std::sort(cc.begin(),cc.end());
   if(rr.size()!=cc.size()||cc.empty()){balanced=false;break;}
   components.push_back({rr,cc});
  }
  if(!balanced||components.size()<2)continue;
  for(const auto&component:components){
   Block block;int size=component.first.size();block.matrix=Matrix(size,std::vector<Series>(size));
   for(int i=0;i<size;++i){block.rows.push_back(row_labels[cut][component.first[i]]);block.cols.push_back(column_labels[cut][component.second[i]]);
    for(int j=0;j<size;++j)block.matrix[i][j]=matrix_by_cut[cut][component.first[i]][component.second[j]];
   }
   blocks_by_cut[cut].push_back(std::move(block));
  }
  matrix_by_cut[cut].clear();row_labels[cut].clear();column_labels[cut].clear();
 }
}
}

extern "C" int gf4_set_seed(int q,const std::uint64_t* coefficients,int words_per_value,int endpoint_check){
 try{
  using namespace f4;
  if(q<2||q>32||words_per_value<1||words_per_value>words)return 1;
  f4::alphabet=q;int n=q*q;f4::blocks_by_cut={};
  for(int cut=0;cut<2;++cut){f4::matrix_by_cut[cut]=f4::Matrix(n,std::vector<f4::Series>(n));f4::row_labels[cut].clear();f4::column_labels[cut].clear();
   for(int a=0;a<q;++a)for(int b=0;b<q;++b){f4::row_labels[cut].push_back(cut==0?std::array<int,4>{a,b,0,0}:std::array<int,4>{a,0,b,0});f4::column_labels[cut].push_back(cut==0?std::array<int,4>{0,0,a,b}:std::array<int,4>{0,a,0,b});}
  }
  for(int v=0;v<q;++v)for(int u=0;u<q;++u)for(int y=0;y<q;++y)for(int x=0;x<q;++x){
   int index=x+q*(y+q*(u+q*v));f4::Poly p;unsigned parity=0;
   for(int k=0;k<words_per_value;++k){p.lo.w[k]=coefficients[(index*2)*words_per_value+k];p.hi.w[k]=coefficients[(index*2+1)*words_per_value+k];parity^=__builtin_parityll(p.lo.w[k]);parity^=2*__builtin_parityll(p.hi.w[k]);}
   if(endpoint_check&&parity!=unsigned(x==v&&y==u))return 2;
   f4::matrix_by_cut[0][q*x+y][q*u+v]={p,0};f4::matrix_by_cut[1][q*x+u][q*y+v]={p,0};
  }
  f4::split_support_components();return 0;
 }catch(...){return 3;}
}

extern "C" int gf4_set_paired_seed(int bulk_modes,const std::uint64_t*coefficients,int words_per_value){
 try{
  if(bulk_modes<1||bulk_modes>2||words_per_value<1||words_per_value>words)return 1;
  int one_q=1<<(bulk_modes+1),q=1<<(2*bulk_modes+1),n=q*q,bulk_mask=(1<<bulk_modes)-1;
  std::vector<f4::Poly>one(one_q*one_q*one_q*one_q);
  for(int i=0;i<int(one.size());++i)for(int k=0;k<words_per_value;++k){one[i].lo.w[k]=coefficients[(i*2)*words_per_value+k];one[i].hi.w[k]=coefficients[(i*2+1)*words_per_value+k];}
  f4::alphabet=q;f4::matrix_by_cut={};f4::row_labels={};f4::column_labels={};f4::blocks_by_cut={};
  auto charge=[&](int a){int core=a>>(2*bulk_modes)&1;return ((__builtin_parity(a&bulk_mask)^core)|((__builtin_parity((a>>bulk_modes)&bulk_mask)^core)<<1));};
  std::array<std::vector<int>,4>indices;std::vector<int>position(n),signature(n);
  for(int a=0;a<q;++a)for(int b=0;b<q;++b){int i=q*a+b,s=charge(a)^charge(b);signature[i]=s;position[i]=indices[s].size();indices[s].push_back(i);}
  for(int cut=0;cut<2;++cut){
   f4::blocks_by_cut[cut].resize(4);
   for(int s=0;s<4;++s){auto&block=f4::blocks_by_cut[cut][s];int size=indices[s].size();block.matrix=f4::Matrix(size,std::vector<f4::Series>(size));
    for(int i:indices[s]){int a=i/q,b=i%q;block.rows.push_back(cut==0?std::array<int,4>{a,b,0,0}:std::array<int,4>{a,0,b,0});block.cols.push_back(cut==0?std::array<int,4>{0,0,a,b}:std::array<int,4>{0,a,0,b});}
   }
  }
  std::vector<int>first(q),second(q);
  for(int a=0;a<q;++a){int core=a>>(2*bulk_modes);first[a]=(a&bulk_mask)|(core<<bulk_modes);second[a]=((a>>bulk_modes)&bulk_mask)|(core<<bulk_modes);}
  for(int v=0;v<q;++v)for(int u=0;u<q;++u)for(int y=0;y<q;++y)for(int x=0;x<q;++x){
   auto a=one[first[x]+one_q*(first[y]+one_q*(first[u]+one_q*first[v]))];
   auto b=one[second[x]+one_q*(second[y]+one_q*(second[u]+one_q*second[v]))];b.lo^=b.hi;
   auto p=f4::multiply(a,b);if(!p)continue;
   for(int cut=0;cut<2;++cut){int r=cut?q*x+u:q*x+y,c=cut?q*y+v:q*u+v,s=signature[r];
    if(signature[c]!=s)return 2;f4::blocks_by_cut[cut][s].matrix[position[r]][position[c]]={p,0};
   }
  }
  return 0;
 }catch(...){return 3;}
}

extern "C" int gf4_paired_entry(int x,int y,int u,int v,std::uint64_t*output){
 try{
  int q=f4::alphabet;if(f4::blocks_by_cut[0].empty())return 1;int r=q*x+y,c=q*u+v;
  for(const auto&block:f4::blocks_by_cut[0]){
   int ri=-1,ci=-1;
   for(int i=0;i<int(block.rows.size());++i){if(q*block.rows[i][0]+block.rows[i][1]==r)ri=i;if(q*block.cols[i][2]+block.cols[i][3]==c)ci=i;}
   if(ri>=0&&ci>=0){auto p=block.matrix[ri][ci].bits;for(int k=0;k<words;++k){output[k]=p.lo.w[k];output[words+k]=p.hi.w[k];}return 0;}
  }
  std::fill(output,output+2*words,0);return 0;
 }catch(...){return 2;}
}

extern "C" int gf4_cost_gradient(const double*params,const double*pi,int cut,double*cost,double*gradient,double*spectra){
 try{
  if(!f4::alphabet||cut<0||cut>1)return 1;int q=f4::alphabet,stride=4*(q-1),scale_index=2*stride,n=q*q;double scale=params[scale_index];if(scale<=0)return 2;
  std::fill(gradient,gradient+scale_index+1,0);*cost=0;
  for(int side=0;side<2;++side){auto v=f4::profile(params+stride*side,side==0?cut:1-cut,side==0,scale);if(int(v.size())!=n)return 3;
   double smallest=1e100;for(auto p:v)smallest=std::min(smallest,p.value);
   double total=0;for(auto p:v)total+=std::exp(-(p.value-smallest));*cost+=-smallest+std::log(total);
   for(int i=0;i<n;++i){auto p=v[i];if(spectra)spectra[n*side+i]=p.value;double weight=std::exp(-(p.value-smallest))/total;gradient[scale_index]-=weight*p.base;
    for(int l=0;l<4;++l)if(p.labels[l])gradient[side*stride+(q-1)*l+p.labels[l]-1]-=weight;}
   for(int l=0;l<4;++l)for(int a=1;a<q;++a){int i=side*stride+(q-1)*l+a-1;*cost+=pi[a]*params[i];gradient[i]+=pi[a];}
  }
  return 0;
 }catch(...){return 4;}
}

extern "C" int gf4_prefixes(const double*params,int cut,double*orders,std::int64_t*bases,std::uint64_t*row_masks,std::uint64_t*column_masks){
 try{
  int q=f4::alphabet;if(q<2||q>8||cut<0||cut>1)return 1;int stride=4*(q-1),n=q*q;double scale=params[2*stride];if(scale<=0)return 2;
  for(int side=0;side<2;++side){int current=side==0?cut:1-cut;auto pivots=f4::profile(params+stride*side,current,side==0,scale);if(int(pivots.size())!=n)return 3;
   std::uint64_t rows=0,cols=0;std::int64_t base=0;
   for(int i=0;i<n;++i){auto p=pivots[i];int r=current?q*p.labels[0]+p.labels[2]:q*p.labels[0]+p.labels[1];int c=current?q*p.labels[1]+p.labels[3]:q*p.labels[2]+p.labels[3];
    if((rows>>r&1)||(cols>>c&1))return 4;rows|=std::uint64_t(1)<<r;cols|=std::uint64_t(1)<<c;base+=p.base;int index=side*n+i;
    orders[index]=p.value;bases[index]=base;row_masks[index]=rows;column_masks[index]=cols;
   }
  }return 0;
 }catch(...){return 5;}
}
