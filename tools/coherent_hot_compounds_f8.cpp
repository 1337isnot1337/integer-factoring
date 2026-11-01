// F8-weighted coherent maps with three binary polynomial coefficient planes.
// Derived from coherent_hot_compounds_fast.cpp; same exact degree bound.
// K_ij=(1+t)t^(i xor j) within a side and K_ij=[i=j] across sides.
// Every B-minor is t^e times its K-minor, and deg(det K)<=e.
#include "selected_gaussian_coefficients.cpp"

Poly hot_determinant(std::vector<std::vector<Poly>> a,int bound){
 if(bound>=precision)throw std::runtime_error("hot determinant degree exceeds precision");
 const int n=int(a.size());
 Poly determinant(1);
 for(int k=0;k<n;++k){
  int pi=-1,pj=-1,v=precision;
  for(int i=k;i<n;++i)for(int j=k;j<n;++j)if(a[i][j]){
   int candidate=valuation(a[i][j]);if(candidate<v){v=candidate;pi=i;pj=j;}
  }
  if(pi<0)return Poly();
  std::swap(a[k],a[pi]);for(auto&row:a)std::swap(row[k],row[pj]);
  determinant=multiply(determinant,a[k][k]);
  Poly unit=a[k][k]>>v,inverse(1);
  if(unit!=Poly(1))for(int width=1;width<precision;width*=2)
   inverse=multiply(unit,multiply(inverse,inverse));
  for(int i=k+1;i<n;++i)if(a[i][k]){
   Poly factor=multiply(a[i][k]>>v,inverse);
   for(int j=k+1;j<n;++j)if(a[k][j])a[i][j]^=multiply(factor,a[k][j]);
  }
 }
 if(degree(determinant)>bound)throw std::runtime_error("hot determinant exceeds proved degree bound");
 return determinant;
}


namespace f8 {
struct Poly {
 ::Poly lo,hi,top;
 Poly(unsigned c=0):lo(c&1),hi((c>>1)&1),top((c>>2)&1){}
 explicit operator bool()const{return bool(lo)||bool(hi)||bool(top);}
 Poly& operator^=(const Poly&b){lo^=b.lo;hi^=b.hi;top^=b.top;return *this;}
 Poly& operator&=(const ::Poly&b){lo&=b;hi&=b;top&=b;return *this;}
 friend Poly operator^(Poly a,const Poly&b){return a^=b;}
 friend bool operator==(const Poly&a,const Poly&b){return a.lo==b.lo&&a.hi==b.hi&&a.top==b.top;}
 friend bool operator!=(const Poly&a,const Poly&b){return !(a==b);}
 Poly operator<<(int s)const{Poly r;r.lo=lo<<s;r.hi=hi<<s;r.top=top<<s;return r;}
 Poly operator>>(int s)const{Poly r;r.lo=lo>>s;r.hi=hi>>s;r.top=top>>s;return r;}
 unsigned constant()const{return unsigned(lo.w[0]&1)|unsigned((hi.w[0]&1)<<1)|unsigned((top.w[0]&1)<<2);}
};
int valuation(const Poly&a){return std::min({a.lo?::valuation(a.lo):precision,a.hi?::valuation(a.hi):precision,a.top?::valuation(a.top):precision});}
int degree(const Poly&a){return std::max({::degree(a.lo),::degree(a.hi),::degree(a.top)});}
Poly multiply(const Poly&a,const Poly&b){
 Poly out;
 if(!a.hi&&!a.top&&!b.hi&&!b.top){out.lo=::multiply(a.lo,b.lo);return out;}
 auto d0=::multiply(a.lo,b.lo),d1=::multiply(a.lo,b.hi)^::multiply(a.hi,b.lo);
 auto d2=::multiply(a.lo,b.top)^::multiply(a.hi,b.hi)^::multiply(a.top,b.lo);
 auto d3=::multiply(a.hi,b.top)^::multiply(a.top,b.hi),d4=::multiply(a.top,b.top);
 out.lo=d0^d3;out.hi=d1^d3^d4;out.top=d2^d4;return out;
}
struct Series{Poly bits;int base=0;};
using Matrix=std::vector<std::vector<Series>>;
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
}

#ifndef HOT_PROFILE_CUTOFF
#define HOT_PROFILE_CUTOFF 48
#endif
std::vector<int> hot_partial_profile(f8::Matrix a){
 const int n=int(a.size()),cutoff=HOT_PROFILE_CUTOFF;const double cap=cutoff+1;
 if(cutoff>=precision)throw std::runtime_error("partial precision");
 for(auto&row:a)for(auto&entry:row)entry=f8::normalize(entry,0,1,cap);
 std::vector<int>out;
 for(int k=0;k<n;++k){
  int best=cutoff+1,pi=-1,pj=-1;
  for(int i=k;i<n;++i)for(int j=k;j<n;++j)
   if(a[i][j].bits&&a[i][j].base<best){best=a[i][j].base;pi=i;pj=j;}
  if(pi<0)break;
  std::swap(a[k],a[pi]);for(auto&row:a)std::swap(row[k],row[pj]);
  out.push_back(best);auto unit=a[k][k].bits;unsigned inverses[8]={0,1,5,6,7,2,3,4};f8::Poly inv(inverses[unit.constant()]);
  if(unit!=f8::Poly(1))for(int p=1;p<precision;p*=2)inv=f8::multiply(unit,f8::multiply(inv,inv));
  for(int j=k;j<n;++j)if(a[k][j].bits){a[k][j].bits=f8::multiply(a[k][j].bits,inv);a[k][j]=f8::normalize(a[k][j],0,1,cap);}
  std::vector<int>support;for(int j=k+1;j<n;++j)if(a[k][j].bits)support.push_back(j);
  for(int i=k+1;i<n;++i)if(a[i][k].bits){auto factor=a[i][k];factor.base-=best;
   for(int j:support){f8::Series p{f8::multiply(factor.bits,a[k][j].bits),factor.base+a[k][j].base};a[i][j]=f8::add(a[i][j],p,0,1,cap);}a[i][k]={};
  }
 }
 return out;
}

extern "C" int coherent_hot_compound_f8(
 int modes,int count,
 const std::uint32_t*row_offsets,const std::uint32_t*row_masks,const std::uint32_t*row_scalars,
 const std::uint32_t*column_offsets,const std::uint32_t*column_masks,const std::uint32_t*column_scalars,
 std::int64_t*orders,std::uint64_t*coefficients){
 try{
  if(modes<1||modes>16||count<1||count>4096)return 1;
  auto scalar_mul=[](unsigned a,unsigned b){unsigned z=0;for(int i=0;i<3;++i){if(b&1)z^=a;b>>=1;a<<=1;if(a&8)a^=11;}return z;};
  const int n=2*modes;
  std::vector<std::vector<Poly>>one(n,std::vector<Poly>(n));
  for(int i=0;i<n;++i)for(int j=0;j<n;++j){
   if(i/modes==j/modes)one[i][j]=Poly(3)<<((i%modes)^(j%modes));
   else if(i%modes==j%modes)one[i][j]=Poly(1);
  }
  const std::uint32_t low_mask=(std::uint32_t(1)<<modes)-1;
  int charge=-1;
  std::vector<std::vector<std::vector<int>>>rr(count),cc(count);
  std::vector<int>row_index_sum(count,-1),column_index_sum(count,-1);
  for(int i=0;i<count;++i){
   if(row_offsets[i]>=row_offsets[i+1]||column_offsets[i]>=column_offsets[i+1])return 2;
   for(int side=0;side<2;++side){
    auto offsets=side?column_offsets:row_offsets;
    auto masks=side?column_masks:row_masks;
    auto&selected_vectors=side?cc[i]:rr[i];
    auto&index_sum=side?column_index_sum[i]:row_index_sum[i];
    for(auto a=offsets[i];a<offsets[i+1];++a){
     std::vector<int>selected;int total=0;
     for(int j=0;j<n;++j)if((masks[a]>>j)&1){selected.push_back(j);total+=j%modes;}
     if(charge<0)charge=int(selected.size());
     if(int(selected.size())!=charge)return 2;
     if(index_sum<0)index_sum=total;
     if(index_sum!=total)return 2;
     selected_vectors.push_back(std::move(selected));
    }
   }
  }
  f8::Matrix matrix(count,std::vector<f8::Series>(count));
  std::vector<std::array<int,4>>row_labels(count),column_labels(count);
  std::vector<double>weights(count,0.0);
  for(int i=0;i<count;++i){
   row_labels[i]={i,0,0,0};column_labels[i]={0,0,i,0};
   for(int j=0;j<count;++j){
    const int exponent=row_index_sum[i]+column_index_sum[j]+charge;
    if(exponent>HOT_PROFILE_CUTOFF)continue;
    f8::Poly value;
    for(std::size_t ri=0;ri<rr[i].size();++ri)for(std::size_t ci=0;ci<cc[j].size();++ci){
     const auto&r=rr[i][ri];const auto&c=cc[j][ci];
     auto rc=row_scalars[row_offsets[i]+ri],ccv=column_scalars[column_offsets[j]+ci];
     if(rc<1||rc>7||ccv<1||ccv>7)return 10;
     auto factor=scalar_mul(rc,ccv);
     std::vector<std::vector<Poly>>minor(charge,std::vector<Poly>(charge));
     for(int a=0;a<charge;++a)for(int b=0;b<charge;++b)minor[a][b]=one[r[a]][c[b]];
     auto term=hot_determinant(std::move(minor),exponent);
     if(factor&1)value.lo^=term;if(factor&2)value.hi^=term;if(factor&4)value.top^=term;
    }
    int endpoint=0;
    for(int a=0;a<words;++a)endpoint^=__builtin_parityll(value.lo.w[a])|(__builtin_parityll(value.hi.w[a])<<1)|(__builtin_parityll(value.top.w[a])<<2);
    if(endpoint!=int(i==j))return 3;
    int expected=0;
    for(auto a=row_offsets[i];a<row_offsets[i+1];++a)
     for(auto b=column_offsets[j];b<column_offsets[j+1];++b){
      auto c=column_masks[b];auto swapped=(c>>modes)|((c&low_mask)<<modes);
      if(row_masks[a]==swapped)expected^=scalar_mul(row_scalars[a],column_scalars[b]);
     }
    if(expected!=int(i==j))return 4;
    matrix[i][j].bits=value;matrix[i][j].base=exponent;
    if(coefficients){
     if(degree(value)+exponent>=precision)return 9;
     auto full=value<<exponent;
     for(int a=0;a<words;++a){coefficients[3*(std::size_t(i)*count+j)*words+a]=full.lo.w[a];
      coefficients[(3*(std::size_t(i)*count+j)+1)*words+a]=full.hi.w[a];coefficients[(3*(std::size_t(i)*count+j)+2)*words+a]=full.top.w[a];}
    }
   }
  }
  auto retained=hot_partial_profile(std::move(matrix));
  for(int i=0;i<count;++i)orders[i]=i<int(retained.size())?retained[i]:HOT_PROFILE_CUTOFF+1;
  return 0;
 }catch(...){return 8;}
}
