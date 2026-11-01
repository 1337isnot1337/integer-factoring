// Height-homogeneous F16 compounds allowing mixed occupation.
// Terms of unequal occupation pair to zero; equal terms use the exact
// binary minor, with exponent one half the common row/column height sum.
#include "tensor_hot_compounds_f16_fast.cpp"

int coherent_hot_compound_f16_mixed_impl(
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
  std::vector<std::vector<std::vector<int>>>rr(count),cc(count);
  std::vector<int>row_height(count,-1),column_height(count,-1);
  for(int i=0;i<count;++i){
   if(row_offsets[i]>=row_offsets[i+1]||column_offsets[i]>=column_offsets[i+1])return 2;
   for(int side=0;side<2;++side){
    auto offsets=side?column_offsets:row_offsets;auto masks=side?column_masks:row_masks;
    auto&selected=side?cc[i]:rr[i];auto&vector_height=side?column_height[i]:row_height[i];
    for(auto a=offsets[i];a<offsets[i+1];++a){
     std::vector<int>indices;int total=0;
     for(int j=0;j<n;++j)if((masks[a]>>j)&1){indices.push_back(j);total+=j%modes;}
     const int height=2*total+int(indices.size());
     if(vector_height<0)vector_height=height;if(vector_height!=height)return 2;
     selected.push_back(std::move(indices));
    }
   }
  }
  tensor_f16::Matrix matrix(count,std::vector<tensor_f16::Series>(count));
  for(int i=0;i<count;++i)for(int j=0;j<count;++j){
   const int height_sum=row_height[i]+column_height[j];
   if(height_sum&1)return 2;const int exponent=height_sum/2;
   tensor_f16::Poly value;
   for(std::size_t ri=0;ri<rr[i].size();++ri)for(std::size_t ci=0;ci<cc[j].size();++ci){
    auto rc=row_scalars[row_offsets[i]+ri],ccv=column_scalars[column_offsets[j]+ci];
    if(rc<1||rc>15||ccv<1||ccv>15)return 10;
    auto factor=tensor_f16::scalar_mul(rc,ccv);
    const auto&r=rr[i][ri];const auto&c=cc[j][ci];
    if(r.size()!=c.size())continue;
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

extern "C" int coherent_hot_compound_f16_mixed(
 int modes,int count,
 const std::uint32_t*row_offsets,const std::uint32_t*row_masks,const std::uint32_t*row_scalars,
 const std::uint32_t*column_offsets,const std::uint32_t*column_masks,const std::uint32_t*column_scalars,
 std::int64_t*orders,std::uint64_t*coefficients){
 return coherent_hot_compound_f16_mixed_impl(modes,count,row_offsets,row_masks,row_scalars,
   column_offsets,column_masks,column_scalars,orders,coefficients,0);
}

extern "C" int coherent_hot_compound_f16_mixed_capped(
 int modes,int count,
 const std::uint32_t*row_offsets,const std::uint32_t*row_masks,const std::uint32_t*row_scalars,
 const std::uint32_t*column_offsets,const std::uint32_t*column_masks,const std::uint32_t*column_scalars,
 std::int64_t*orders,std::uint64_t*coefficients,int order_cap){
 if(order_cap<1)return 11;
 return coherent_hot_compound_f16_mixed_impl(modes,count,row_offsets,row_masks,row_scalars,
   column_offsets,column_masks,column_scalars,orders,coefficients,order_cap);
}
