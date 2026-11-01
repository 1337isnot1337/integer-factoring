// Literal Pfaffian contractions with termwise Laurent monomial gauges.
// A reported common clearing power is retained in every returned profile.
#include "tensor_literal_hot_pfaffians_f16.cpp"

extern "C" int tensor_literal_weighted_block_f16(
 int modes,int count,
 const std::uint32_t*row_offsets,const std::uint32_t*row_masks,const std::uint32_t*row_scalars,
 const std::int32_t*row_shifts,
 const std::uint32_t*column_offsets,const std::uint32_t*column_masks,const std::uint32_t*column_scalars,
 const std::int32_t*column_shifts,
 int clearing_power,std::int64_t*orders,std::uint64_t*coefficients,
 int order_cap,std::uint64_t budget){
 try{
  if(modes<1||modes>16||count<1||count>4608||order_cap<0||clearing_power<0)return 1;
  tensor_literal_pf::initialize(modes,std::size_t(budget));
  const std::uint32_t low=(std::uint32_t(1)<<modes)-1;
  int parity=-1,maxheight[2]={-1000000,-1000000},minshift[2]={1000000,1000000};
  for(int i=0;i<count;++i)for(int side=0;side<2;++side){
   auto offsets=side?column_offsets:row_offsets;auto masks=side?column_masks:row_masks;
   auto scalars=side?column_scalars:row_scalars;auto shifts=side?column_shifts:row_shifts;
   if(offsets[i]>=offsets[i+1])return 2;
   for(auto p=offsets[i];p<offsets[i+1];++p){
    if(scalars[p]<1||scalars[p]>15)return 10;
    int current=0;for(int bit=0;bit<2*modes;++bit)if((masks[p]>>bit)&1)current+=2*(bit%modes)+1;
    if(parity<0)parity=current&1;if(parity!=(current&1))return 2;
    maxheight[side]=std::max(maxheight[side],current+shifts[p]);
    minshift[side]=std::min(minshift[side],int(shifts[p]));
   }
  }
  if(minshift[0]+minshift[1]+clearing_power<0)return 13;
  if(coefficients&&maxheight[0]+maxheight[1]+clearing_power>=precision)return 9;
  tensor_f16::Matrix matrix(count,std::vector<tensor_f16::Series>(count));
  for(int i=0;i<count;++i)for(int j=0;j<count;++j){
   tensor_f16::Poly value;unsigned endpoint=0;
   for(auto a=row_offsets[i];a<row_offsets[i+1];++a)for(auto b=column_offsets[j];b<column_offsets[j+1];++b){
    auto factor=tensor_f16::scalar_mul(row_scalars[a],column_scalars[b]);
    int shift=row_shifts[a]+column_shifts[b]+clearing_power;
    if(shift<0)return 13;
    if(shift<precision){
     std::uint64_t mask=std::uint64_t(row_masks[a])|(std::uint64_t(column_masks[b])<<(2*modes));
     auto term=tensor_literal_pf::pfaffian(mask)<<shift;
     for(int plane=0;plane<4;++plane)if((factor>>plane)&1)value.plane(plane)^=term;
    }
    auto column=column_masks[b],swapped=(column>>modes)|((column&low)<<modes);
    if(row_masks[a]==swapped)endpoint^=factor;
   }
   if(endpoint!=unsigned(i==j))return 4;
   matrix[i][j].bits=value;matrix[i][j].base=0;
   if(coefficients)for(int plane=0;plane<4;++plane)for(int k=0;k<words;++k)
    coefficients[(4*(std::size_t(i)*count+j)+plane)*words+k]=value.plane(plane).w[k];
  }
  auto result=tensor_f16::fast_profile(std::move(matrix),order_cap);
  if(int(result.size())!=count)return 6;
  for(int i=0;i<count;++i)orders[i]=result[i];return 0;
 }catch(const std::length_error&){return 12;}catch(...){return 8;}
}
