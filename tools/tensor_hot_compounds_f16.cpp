// Exact coherent compounds over F16=F4[z]/(z^2+z+omega).
// Binary minors and their degree bounds reuse the established F4 constructor.
// The unweighted F16 Smith elimination keeps all four binary planes.
#include "coherent_hot_compounds_f4.cpp"

namespace tensor_f16 {
unsigned scalar_mul(unsigned a,unsigned b){
 auto four=[](unsigned x,unsigned y){unsigned out=0;for(int i=0;i<2;++i){if(y&1)out^=x;y>>=1;x<<=1;if(x&4)x^=7;}return out;};
 unsigned p=four(a&3,b&3),q=four(a>>2,b>>2),r=four((a&3)^(a>>2),(b&3)^(b>>2));
 return (p^four(2,q))|((r^p)<<2);
}
unsigned scalar_inverse(unsigned a){
 if(!a)throw std::runtime_error("zero F16 inverse");
 unsigned out=1,p=a;for(int e=14;e;e>>=1,p=scalar_mul(p,p))if(e&1)out=scalar_mul(out,p);return out;
}
f4::Poly times_omega(const f4::Poly&a){f4::Poly out;out.lo=a.hi;out.hi=a.lo^a.hi;return out;}
struct Poly{
 f4::Poly lo,hi;
 Poly(unsigned coefficient=0):lo(coefficient&3),hi(coefficient>>2){}
 explicit operator bool()const{return bool(lo)||bool(hi);}
 Poly&operator^=(const Poly&b){lo^=b.lo;hi^=b.hi;return *this;}
 Poly&operator&=(const ::Poly&b){lo&=b;hi&=b;return *this;}
 friend Poly operator^(Poly a,const Poly&b){return a^=b;}
 friend bool operator==(const Poly&a,const Poly&b){return a.lo==b.lo&&a.hi==b.hi;}
 friend bool operator!=(const Poly&a,const Poly&b){return !(a==b);}
 Poly operator<<(int s)const{Poly out;out.lo=lo<<s;out.hi=hi<<s;return out;}
 Poly operator>>(int s)const{Poly out;out.lo=lo>>s;out.hi=hi>>s;return out;}
 unsigned constant()const{return lo.constant()|(hi.constant()<<2);}
 ::Poly&plane(int i){return i==0?lo.lo:i==1?lo.hi:i==2?hi.lo:hi.hi;}
 const ::Poly&plane(int i)const{return i==0?lo.lo:i==1?lo.hi:i==2?hi.lo:hi.hi;}
};
int valuation(const Poly&a){return std::min(f4::valuation(a.lo),f4::valuation(a.hi));}
int degree(const Poly&a){return std::max(f4::degree(a.lo),f4::degree(a.hi));}
Poly multiply(const Poly&a,const Poly&b){
 if(!a.hi&&!b.hi){Poly out;out.lo=f4::multiply(a.lo,b.lo);return out;}
 auto p=f4::multiply(a.lo,b.lo),q=f4::multiply(a.hi,b.hi),r=f4::multiply(a.lo^a.hi,b.lo^b.hi);
 Poly out;out.lo=p^times_omega(q);out.hi=r^p;return out;
}
struct Series{Poly bits;int base=0;};
using Matrix=std::vector<std::vector<Series>>;
Series normalize(Series a,int cap){
 if(!a.bits)return {};
 int keep=cap-a.base;if(keep<=0)return {};
 if(keep<precision)a.bits&=::low_mask(keep);
 if(a.bits){int v=valuation(a.bits);a.base+=v;a.bits=a.bits>>v;}else a.base=0;
 return a;
}
Series add(Series a,Series b,int cap){
 if(!a.bits)return normalize(b,cap);if(!b.bits)return normalize(a,cap);
 if(a.base>b.base)std::swap(a,b);a.bits^=b.bits<<(b.base-a.base);return normalize(a,cap);
}
std::vector<int> profile(Matrix a,int order_cap=0){
 const int n=int(a.size());int minimum=precision;
 for(auto&row:a)for(auto&entry:row)if(entry.bits){int v=valuation(entry.bits);entry.bits=entry.bits>>v;entry.base+=v;minimum=std::min(minimum,entry.base);}
 const int cap=minimum+precision;
 for(auto&row:a)for(auto&entry:row)entry=normalize(entry,cap);
 std::vector<int>orders;
 for(int k=0;k<n;++k){
  int best=cap,pi=-1,pj=-1;
  for(int i=k;i<n;++i)for(int j=k;j<n;++j)if(a[i][j].bits&&a[i][j].base<best){best=a[i][j].base;pi=i;pj=j;}
  if(order_cap>0&&best>=order_cap&&order_cap<=cap){
   orders.insert(orders.end(),n-k,order_cap);return orders;
  }
  if(pi<0)return {};
  std::swap(a[k],a[pi]);for(auto&row:a)std::swap(row[k],row[pj]);orders.push_back(best);
  Poly unit=a[k][k].bits,inv(scalar_inverse(unit.constant()));
  if(unit!=Poly(1))for(int width=1;width<precision;width*=2)inv=multiply(unit,multiply(inv,inv));
  for(int j=k;j<n;++j)if(a[k][j].bits){a[k][j].bits=multiply(a[k][j].bits,inv);a[k][j]=normalize(a[k][j],cap);}
  if(a[k][k].bits!=Poly(1)||a[k][k].base!=best)throw std::runtime_error("F16 pivot normalization");
  std::vector<int>support;for(int j=k+1;j<n;++j)if(a[k][j].bits)support.push_back(j);
  for(int i=k+1;i<n;++i)if(a[i][k].bits){auto factor=a[i][k];factor.base-=best;
   for(int j:support){Series term{multiply(factor.bits,a[k][j].bits),factor.base+a[k][j].base};a[i][j]=add(a[i][j],term,cap);}a[i][k]={};
  }
 }
 return orders;
}
}

extern "C" int tensor_f16_scalar_control(unsigned a,unsigned b){
 if(a>15||b>15)return -1;return int(tensor_f16::scalar_mul(a,b));
}

int coherent_hot_compound_f16_impl(
 int modes,int count,
 const std::uint32_t*row_offsets,const std::uint32_t*row_masks,const std::uint32_t*row_scalars,
 const std::uint32_t*column_offsets,const std::uint32_t*column_masks,const std::uint32_t*column_scalars,
 std::int64_t*orders,std::uint64_t*coefficients,int order_cap){
 try{
  if(modes<1||modes>16||count<1||count>4096)return 1;
  const int n=2*modes;
  std::vector<std::vector<::Poly>>one(n,std::vector<::Poly>(n));
  for(int i=0;i<n;++i)for(int j=0;j<n;++j){
   if(i/modes==j/modes)one[i][j]=::Poly(3)<<((i%modes)^(j%modes));
   else if(i%modes==j%modes)one[i][j]=::Poly(1);
  }
  const std::uint32_t low_mask=(std::uint32_t(1)<<modes)-1;
  int charge=-1;
  std::vector<std::vector<std::vector<int>>>rr(count),cc(count);
  std::vector<int>row_index_sum(count,-1),column_index_sum(count,-1);
  for(int i=0;i<count;++i){
   if(row_offsets[i]>=row_offsets[i+1]||column_offsets[i]>=column_offsets[i+1])return 2;
   for(int side=0;side<2;++side){
    auto offsets=side?column_offsets:row_offsets;auto masks=side?column_masks:row_masks;
    auto&selected=side?cc[i]:rr[i];auto&index_sum=side?column_index_sum[i]:row_index_sum[i];
    for(auto a=offsets[i];a<offsets[i+1];++a){
     std::vector<int>indices;int total=0;
     for(int j=0;j<n;++j)if((masks[a]>>j)&1){indices.push_back(j);total+=j%modes;}
     if(charge<0)charge=int(indices.size());
     if(int(indices.size())!=charge)return 2;
     if(index_sum<0)index_sum=total;if(index_sum!=total)return 2;
     selected.push_back(std::move(indices));
    }
   }
  }
  tensor_f16::Matrix matrix(count,std::vector<tensor_f16::Series>(count));
  for(int i=0;i<count;++i)for(int j=0;j<count;++j){
   const int exponent=row_index_sum[i]+column_index_sum[j]+charge;
   tensor_f16::Poly value;
   for(std::size_t ri=0;ri<rr[i].size();++ri)for(std::size_t ci=0;ci<cc[j].size();++ci){
    auto rc=row_scalars[row_offsets[i]+ri],ccv=column_scalars[column_offsets[j]+ci];
    if(rc<1||rc>15||ccv<1||ccv>15)return 10;
    auto factor=tensor_f16::scalar_mul(rc,ccv);
    const auto&r=rr[i][ri];const auto&c=cc[j][ci];
    std::vector<std::vector<::Poly>>minor(charge,std::vector<::Poly>(charge));
    for(int a=0;a<charge;++a)for(int b=0;b<charge;++b)minor[a][b]=one[r[a]][c[b]];
    auto term=hot_determinant(std::move(minor),exponent);
    for(int plane=0;plane<4;++plane)if((factor>>plane)&1)value.plane(plane)^=term;
   }
   unsigned endpoint=0;
   for(int plane=0;plane<4;++plane)for(int k=0;k<words;++k)endpoint^=unsigned(__builtin_parityll(value.plane(plane).w[k]))<<plane;
   if(endpoint!=unsigned(i==j))return 3;
   unsigned expected=0;
   for(auto a=row_offsets[i];a<row_offsets[i+1];++a)for(auto b=column_offsets[j];b<column_offsets[j+1];++b){
    auto c=column_masks[b],swapped=(c>>modes)|((c&low_mask)<<modes);
    if(row_masks[a]==swapped)expected^=tensor_f16::scalar_mul(row_scalars[a],column_scalars[b]);
   }
   if(expected!=unsigned(i==j))return 4;
   matrix[i][j].bits=value;matrix[i][j].base=exponent;
   if(coefficients){
    if(tensor_f16::degree(value)+exponent>=precision)return 9;
    auto full=value<<exponent;
    for(int plane=0;plane<4;++plane)for(int k=0;k<words;++k)
     coefficients[(4*(std::size_t(i)*count+j)+plane)*words+k]=full.plane(plane).w[k];
   }
  }
  auto result=tensor_f16::profile(std::move(matrix),order_cap);
  if(int(result.size())!=count)return 6;
  for(int i=0;i<count;++i){if(result[i]<0)return 7;orders[i]=result[i];}
  return 0;
 }catch(...){return 8;}
}

extern "C" int coherent_hot_compound_f16(
 int modes,int count,
 const std::uint32_t*row_offsets,const std::uint32_t*row_masks,const std::uint32_t*row_scalars,
 const std::uint32_t*column_offsets,const std::uint32_t*column_masks,const std::uint32_t*column_scalars,
 std::int64_t*orders,std::uint64_t*coefficients){
 return coherent_hot_compound_f16_impl(modes,count,row_offsets,row_masks,row_scalars,
   column_offsets,column_masks,column_scalars,orders,coefficients,0);
}

extern "C" int coherent_hot_compound_f16_capped(
 int modes,int count,
 const std::uint32_t*row_offsets,const std::uint32_t*row_masks,const std::uint32_t*row_scalars,
 const std::uint32_t*column_offsets,const std::uint32_t*column_masks,const std::uint32_t*column_scalars,
 std::int64_t*orders,std::uint64_t*coefficients,int order_cap){
 if(order_cap<1)return 11;
 return coherent_hot_compound_f16_impl(modes,count,row_offsets,row_masks,row_scalars,
   column_offsets,column_masks,column_scalars,orders,coefficients,order_cap);
}
