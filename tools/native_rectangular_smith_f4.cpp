// Optional rectangular GF(4) weighted Smith search engine. Not an exact real-weight certificate.
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

#ifndef F4_SERIES_WORDS
#define F4_SERIES_WORDS 8
#endif
constexpr int words=F4_SERIES_WORDS,precision=64*words;
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
std::array<int,4>alphabets{},offsets{};
int parameter_stride=0;
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
double height(const double*g,const std::array<int,4>&labels){double h=0;for(int l=0;l<4;++l)if(labels[l])h+=g[offsets[l]+labels[l]-1];return h;}
std::vector<Pivot> profile_block(const double*g,const Matrix&matrix,
                               std::vector<std::array<int,4>>rows,
                               std::vector<std::array<int,4>>cols,bool infinity,double scale,
                               int wanted_rank=-1,const std::vector<double>*given_rows=nullptr,
                               const std::vector<double>*given_cols=nullptr){
 auto a=matrix;int n=a.size(),rank=wanted_rank<0?n:wanted_rank;
 if(rank<0||rank>n)throw std::runtime_error("rank");
 std::vector<double>r(n),c(n);for(int i=0;i<n;++i){r[i]=given_rows?(*given_rows)[i]:height(g,rows[i]);c[i]=given_cols?(*given_cols)[i]:height(g,cols[i]);}
 if(infinity){for(auto&row:a)for(auto&entry:row)if(entry.bits){
  int d=degree(entry.bits);Poly rev;
  for(int k=0;k<=d;++k){
   if(entry.bits.lo.w[k/64]>>(k%64)&1)rev.lo.w[(d-k)/64]|=std::uint64_t(1)<<((d-k)%64);
   if(entry.bits.hi.w[k/64]>>(k%64)&1)rev.hi.w[(d-k)/64]|=std::uint64_t(1)<<((d-k)%64);
  }
  entry.bits=rev;entry.base=-d-entry.base;
 }}else{for(auto&row:a)for(auto&entry:row)if(entry.bits){int v=valuation(entry.bits);entry.bits=entry.bits>>v;entry.base+=v;}}
 double minimum=1e100;for(int i=0;i<n;++i)for(int j=0;j<n;++j)if(a[i][j].bits)minimum=std::min(minimum,scale*a[i][j].base+r[i]+c[j]);
 double cap=minimum+precision*scale;
 for(int i=0;i<n;++i)for(int j=0;j<n;++j)a[i][j]=normalize(a[i][j],r[i]+c[j],scale,cap);
 std::vector<Pivot>out;
 for(int k=0;k<rank;++k){
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
 int n=alphabets[0]*alphabets[1];
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

extern "C" int rectangular_f4_set_seed(int A,int B,const std::uint64_t*coefficients,int words_per_value,int endpoint_check){
 try{
  using namespace f4;
  if(A<2||B<2||A>64||B>64||words_per_value<1||words_per_value>words)return 1;
  f4::alphabet=A;alphabets={A,B,B,A};parameter_stride=0;blocks_by_cut={};
  for(int l=0;l<4;++l){offsets[l]=parameter_stride;parameter_stride+=alphabets[l]-1;}
  int n=A*B;
  for(int cut=0;cut<2;++cut){matrix_by_cut[cut]=Matrix(n,std::vector<Series>(n));row_labels[cut].clear();column_labels[cut].clear();
   for(int x=0;x<A;++x)for(int b=0;b<B;++b)row_labels[cut].push_back(cut==0?std::array<int,4>{x,b,0,0}:std::array<int,4>{x,0,b,0});
   for(int b=0;b<B;++b)for(int v=0;v<A;++v)column_labels[cut].push_back(cut==0?std::array<int,4>{0,0,b,v}:std::array<int,4>{0,b,0,v});
  }
  for(int v=0;v<A;++v)for(int u=0;u<B;++u)for(int y=0;y<B;++y)for(int x=0;x<A;++x){
   int index=x+A*(y+B*(u+B*v));f4::Poly p;unsigned parity=0;
   for(int k=0;k<words_per_value;++k){p.lo.w[k]=coefficients[(index*2)*words_per_value+k];p.hi.w[k]=coefficients[(index*2+1)*words_per_value+k];parity^=__builtin_parityll(p.lo.w[k]);parity^=2*__builtin_parityll(p.hi.w[k]);}
   if(endpoint_check&&parity!=unsigned(x==v&&y==u))return 2;
   matrix_by_cut[0][B*x+y][A*u+v]={p,0};matrix_by_cut[1][B*x+u][A*y+v]={p,0};
  }
  split_support_components();return 0;
 }catch(...){return 3;}
}
extern "C" int rectangular_f4_cost_gradient(const double*params,const double*piX,const double*piY,int cut,double*cost,double*gradient,double*spectra){
 try{
  using namespace f4;
  if(!alphabet||cut<0||cut>1)return 1;int stride=parameter_stride,scale_index=2*stride,n=alphabets[0]*alphabets[1];double scale=params[scale_index];if(scale<=0)return 2;
  std::fill(gradient,gradient+scale_index+1,0);*cost=0;
  for(int side=0;side<2;++side){auto v=profile(params+stride*side,side==0?cut:1-cut,side==0,scale);if(int(v.size())!=n)return 3;
   double smallest=1e100;for(auto p:v)smallest=std::min(smallest,p.value);
   double total=0;for(auto p:v)total+=std::exp(-(p.value-smallest));*cost+=-smallest+std::log(total);
   for(int i=0;i<n;++i){auto p=v[i];if(spectra)spectra[n*side+i]=p.value;double weight=std::exp(-(p.value-smallest))/total;gradient[scale_index]-=weight*p.base;
    for(int l=0;l<4;++l)if(p.labels[l])gradient[side*stride+offsets[l]+p.labels[l]-1]-=weight;}
   for(int l=0;l<4;++l){const double*pi=(l==0||l==3)?piX:piY;
    for(int a=1;a<alphabets[l];++a){int i=side*stride+offsets[l]+a-1;*cost+=pi[a]*params[i];gradient[i]+=pi[a];}}
  }
  return 0;
 }catch(...){return 4;}
}

// Selected row/column indices certify nested attaining minors. The finite
// field arithmetic is exact; the real-weight selection still needs replay.
extern "C" int rectangular_f4_prefixes(const double*params,int cut,double*orders,
 std::int64_t*bases,std::uint32_t*rows,std::uint32_t*cols){
 try{
  using namespace f4;
  if(!alphabet||cut<0||cut>1)return 1;
  int A=alphabets[0],B=alphabets[1],n=A*B,stride=parameter_stride;
  double scale=params[2*stride];if(scale<=0)return 2;
  for(int side=0;side<2;++side){
   int selected=side==0?cut:1-cut;
   auto pivots=profile(params+side*stride,selected,side==0,scale);
   if(int(pivots.size())!=n)return 3;
   std::int64_t cumulative=0;std::vector<bool>seen_r(n),seen_c(n);
   for(int i=0;i<n;++i){auto p=pivots[i];
    int r=B*p.labels[0]+p.labels[selected==0?1:2];
    int c=A*p.labels[selected==0?2:1]+p.labels[3];
    if(seen_r[r]||seen_c[c])return 4;
    seen_r[r]=seen_c[c]=true;cumulative+=p.base;
    int index=side*n+i;orders[index]=p.value;bases[index]=cumulative;
    rows[index]=r;cols[index]=c;
   }
  }
  return 0;
 }catch(...){return 5;}
}
