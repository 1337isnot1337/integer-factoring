// Binary subset-determinant acceleration for exact F16 coherent compounds.
// Existing F16 scalar arithmetic and Smith elimination are reused unchanged.
#include "tensor_hot_compounds_f16.cpp"

namespace tensor_f16 {
// Capped Smith data depend only on the jet modulo t^order_cap.
std::vector<int> fast_profile(Matrix a,int order_cap=0){
 const int n=int(a.size());int minimum=precision;
 for(auto&row:a)for(auto&entry:row)if(entry.bits){int v=valuation(entry.bits);entry.bits=entry.bits>>v;entry.base+=v;minimum=std::min(minimum,entry.base);}
 const int cap=order_cap>0?std::min(minimum+precision,order_cap):minimum+precision;
 for(auto&row:a)for(auto&entry:row)entry=normalize(entry,cap);
 std::vector<int>orders;
 for(int k=0;k<n;++k){
  int best=cap,pi=-1,pj=-1;
  const int lower=orders.empty()?minimum:orders.back();
  bool attained=false;
  for(int i=k;i<n&&!attained;++i)for(int j=k;j<n;++j)if(a[i][j].bits&&a[i][j].base<best){
   best=a[i][j].base;pi=i;pj=j;
   // Minimum-valuation elimination cannot create a smaller valuation.
   if(best==lower){attained=true;break;}
  }
  if(order_cap>0&&best>=order_cap&&order_cap<=cap){
   orders.insert(orders.end(),n-k,order_cap);return orders;
  }
  if(pi<0)return {};
  std::swap(a[k],a[pi]);for(auto&row:a)std::swap(row[k],row[pj]);orders.push_back(best);
  Poly unit=a[k][k].bits,inv(scalar_inverse(unit.constant()));
  if(unit!=Poly(1))for(int width=1;width<cap-best;width*=2){inv=multiply(unit,multiply(inv,inv));inv&=::low_mask(cap-best);}
  for(int j=k;j<n;++j)if(a[k][j].bits){a[k][j].bits=multiply(a[k][j].bits,inv);a[k][j]=normalize(a[k][j],cap);}
  if(a[k][k].bits!=Poly(1)||a[k][k].base!=best)throw std::runtime_error("F16 pivot normalization");
  std::vector<int>support;for(int j=k+1;j<n;++j)if(a[k][j].bits)support.push_back(j);
  for(int i=k+1;i<n;++i)if(a[i][k].bits){auto factor=a[i][k];factor.base-=best;
   for(int j:support){const int base=factor.base+a[k][j].base;if(base>=cap)continue;
    Series term{multiply(factor.bits,a[k][j].bits),base};a[i][j]=add(a[i][j],term,cap);}a[i][k]={};
  }
 }
 return orders;
}
}

::Poly tensor_fast_binary_minor(const std::vector<int>&row,
 const std::vector<int>&column,int modes,int bound,
 const std::vector<std::vector<::Poly>>&one){
 const int n=int(row.size());
 if(bound>=64||n>12){
  std::vector<std::vector<::Poly>>minor(n,std::vector<::Poly>(n));
  for(int i=0;i<n;++i)for(int j=0;j<n;++j)minor[i][j]=one[row[i]][column[j]];
  return hot_determinant(std::move(minor),bound);
 }
 // In characteristic two determinant signs disappear. DP[mask] sums
 // products for assignments of the first popcount(mask) rows to mask.
 // Every partial product has degree <= bound, because each selected
 // entry has degree <= row_index+column_index+1 and all indices >=0.
 std::array<std::uint64_t,1<<12>dp;
 std::fill(dp.begin(),dp.begin()+(1u<<n),0);
 dp[0]=1;
 const unsigned full=(1u<<n)-1;
 for(unsigned mask=0;mask<full;++mask){
  const auto value=dp[mask];if(!value)continue;
  const int i=__builtin_popcount(mask),r=row[i];
  for(unsigned remaining=full^mask;remaining;remaining&=remaining-1){
   const unsigned bit=remaining&-remaining;const int j=__builtin_ctz(bit),c=column[j];
   if(r/modes==c/modes){
    const int shift=(r%modes)^(c%modes);
    dp[mask|bit]^=(value<<shift)^(value<<(shift+1));
   }else if(r%modes==c%modes)dp[mask|bit]^=value;
  }
 }
 ::Poly result;result.w[0]=dp[full];
 if(degree(result)>bound)throw std::runtime_error("F16 subset determinant degree bound");
 return result;
}


extern "C" int tensor_f16_fast_minor_control(int modes,std::uint32_t rmask,
 std::uint32_t cmask,std::uint64_t*out){
 try{
  if(modes<1||modes>16)return 1;
  std::vector<int>r,c;int bound=0;
  for(int i=0;i<2*modes;++i){
   if((rmask>>i)&1){r.push_back(i);bound+=i%modes;}
   if((cmask>>i)&1){c.push_back(i);bound+=i%modes;}
  }
  if(r.size()!=c.size())return 2;bound+=int(r.size());
  const int n=2*modes;
  std::vector<std::vector<::Poly>>one(n,std::vector<::Poly>(n));
  for(int i=0;i<n;++i)for(int j=0;j<n;++j){
   if(i/modes==j/modes)one[i][j]=::Poly(3)<<((i%modes)^(j%modes));
   else if(i%modes==j%modes)one[i][j]=::Poly(1);
  }
  auto result=tensor_fast_binary_minor(r,c,modes,bound,one);
  for(int i=0;i<words;++i)out[i]=result.w[i];return 0;
 }catch(...){return 8;}
}

int coherent_hot_compound_f16_fast_impl(
 int modes,int count,
 const std::uint32_t*row_offsets,const std::uint32_t*row_masks,const std::uint32_t*row_scalars,
 const std::uint32_t*column_offsets,const std::uint32_t*column_masks,const std::uint32_t*column_scalars,
 std::int64_t*orders,std::uint64_t*coefficients,int order_cap){
 try{
  if(modes<1||modes>16||count<1||count>4608)return 1;
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
    auto term=tensor_fast_binary_minor(r,c,modes,exponent,one);
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
  auto result=tensor_f16::fast_profile(std::move(matrix),order_cap);
  if(int(result.size())!=count)return 6;
  for(int i=0;i<count;++i){if(result[i]<0)return 7;orders[i]=result[i];}
  return 0;
 }catch(...){return 8;}
}

extern "C" int coherent_hot_compound_f16_fast(
 int modes,int count,
 const std::uint32_t*row_offsets,const std::uint32_t*row_masks,const std::uint32_t*row_scalars,
 const std::uint32_t*column_offsets,const std::uint32_t*column_masks,const std::uint32_t*column_scalars,
 std::int64_t*orders,std::uint64_t*coefficients){
 return coherent_hot_compound_f16_fast_impl(modes,count,row_offsets,row_masks,row_scalars,
   column_offsets,column_masks,column_scalars,orders,coefficients,0);
}

extern "C" int coherent_hot_compound_f16_fast_capped(
 int modes,int count,
 const std::uint32_t*row_offsets,const std::uint32_t*row_masks,const std::uint32_t*row_scalars,
 const std::uint32_t*column_offsets,const std::uint32_t*column_masks,const std::uint32_t*column_scalars,
 std::int64_t*orders,std::uint64_t*coefficients,int order_cap){
 if(order_cap<1)return 11;
 return coherent_hot_compound_f16_fast_impl(modes,count,row_offsets,row_masks,row_scalars,
   column_offsets,column_masks,column_scalars,orders,coefficients,order_cap);
}
