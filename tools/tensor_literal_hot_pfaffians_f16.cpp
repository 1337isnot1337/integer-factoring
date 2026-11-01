// Literal regular hot-reflection Pfaffians, without deletion-closure assumptions.
// The binary Pfaffian cache stores exact jets modulo t^precision. Endpoint
// identity is checked independently from the exact kernel at t=1.
#include "tensor_hot_compounds_f16_fast.cpp"
#include <unordered_map>

namespace tensor_literal_pf {
struct Edge { int shift=-1; bool binomial=false; };
static int active_modes=0;
static std::size_t cache_budget=8000000;
static std::vector<std::vector<Edge>>edges;
static std::unordered_map<std::uint64_t,::Poly>cache;
static std::uint64_t requests=0,hits=0;

void initialize(int modes,std::size_t budget){
 if(modes<1||modes>16||budget<1)throw std::runtime_error("literal Pfaffian dimensions");
 cache_budget=budget;
 if(active_modes==modes)return;
 active_modes=modes;cache.clear();requests=hits=0;
 edges.assign(4*modes,std::vector<Edge>(4*modes));
 for(int i=0;i<4*modes;++i)for(int j=i+1;j<4*modes;++j){
  int l=i/modes,k=j/modes,a=i%modes,b=j%modes;Edge e;
  if((l==0&&k==1)||(l==2&&k==3))e={2*(a&b),true};
  else if((l==0&&k==2)||(l==1&&k==3))e={2*(a|b)+1,true};
  else if(((l==0&&k==3)||(l==1&&k==2))&&a==b)e={2*a+1,false};
  edges[i][j]=edges[j][i]=e;
 }
}

::Poly pfaffian(std::uint64_t mask){
 ++requests;
 if(!mask)return ::Poly(1);
 if(__builtin_popcountll(mask)&1)return ::Poly();
 auto found=cache.find(mask);if(found!=cache.end()){++hits;return found->second;}
 int first=__builtin_ctzll(mask);std::uint64_t rest=mask&(mask-1),remaining=rest;::Poly result;
 while(remaining){
  auto chosen=remaining&-remaining;remaining^=chosen;int j=__builtin_ctzll(chosen);
  auto e=edges[first][j];if(e.shift<0||e.shift>=precision)continue;
  auto child=pfaffian(rest^chosen);result^=child<<e.shift;
  if(e.binomial)result^=child<<(e.shift+1);
 }
 if(cache.size()>=cache_budget)throw std::length_error("literal Pfaffian cache budget");
 cache.emplace(mask,result);return result;
}
}

extern "C" int tensor_literal_pfaffian_control(int modes,std::uint64_t mask,
 std::uint64_t*out,std::uint64_t budget){
 try{
  tensor_literal_pf::initialize(modes,std::size_t(budget));
  if(modes<16&&(mask>>(4*modes)))return 1;
  auto value=tensor_literal_pf::pfaffian(mask);
  for(int i=0;i<words;++i)out[i]=value.w[i];return 0;
 }catch(const std::length_error&){return 12;}catch(...){return 8;}
}

extern "C" void tensor_literal_pfaffian_cache_info(std::uint64_t*out){
 out[0]=tensor_literal_pf::cache.size();out[1]=tensor_literal_pf::requests;out[2]=tensor_literal_pf::hits;
}

extern "C" int tensor_literal_hot_block_f16(
 int modes,int count,
 const std::uint32_t*row_offsets,const std::uint32_t*row_masks,const std::uint32_t*row_scalars,
 const std::uint32_t*column_offsets,const std::uint32_t*column_masks,const std::uint32_t*column_scalars,
 std::int64_t*orders,std::uint64_t*coefficients,int order_cap,std::uint64_t budget){
 try{
  if(modes<1||modes>16||count<1||count>4608||order_cap<0)return 1;
  tensor_literal_pf::initialize(modes,std::size_t(budget));
  const std::uint32_t low=(std::uint32_t(1)<<modes)-1;
  std::vector<int>rh(count,-1),ch(count,-1);
  int parity=-1;
  for(int i=0;i<count;++i)for(int side=0;side<2;++side){
   auto offsets=side?column_offsets:row_offsets;auto masks=side?column_masks:row_masks;
   auto scalars=side?column_scalars:row_scalars;int&height=side?ch[i]:rh[i];
   if(offsets[i]>=offsets[i+1])return 2;
   for(auto p=offsets[i];p<offsets[i+1];++p){
    if(scalars[p]<1||scalars[p]>15)return 10;
    int current=0;for(int bit=0;bit<2*modes;++bit)if((masks[p]>>bit)&1)current+=2*(bit%modes)+1;
    if(height<0)height=current;if(height!=current)return 2;
    if(parity<0)parity=current&1;if(parity!=(current&1))return 2;
   }
  }
  if(coefficients&&*std::max_element(rh.begin(),rh.end())+*std::max_element(ch.begin(),ch.end())>=precision)return 9;
  tensor_f16::Matrix matrix(count,std::vector<tensor_f16::Series>(count));
  for(int i=0;i<count;++i)for(int j=0;j<count;++j){
   tensor_f16::Poly value;unsigned endpoint=0;
   for(auto a=row_offsets[i];a<row_offsets[i+1];++a)for(auto b=column_offsets[j];b<column_offsets[j+1];++b){
    auto factor=tensor_f16::scalar_mul(row_scalars[a],column_scalars[b]);
    std::uint64_t mask=std::uint64_t(row_masks[a])|(std::uint64_t(column_masks[b])<<(2*modes));
    auto term=tensor_literal_pf::pfaffian(mask);
    for(int plane=0;plane<4;++plane)if((factor>>plane)&1)value.plane(plane)^=term;
    auto column=column_masks[b],swapped=(column>>modes)|((column&low)<<modes);
    if(row_masks[a]==swapped)endpoint^=factor;
   }
   // At t=1 the only kernel edges pair x with v and y with u.
   // This exact endpoint check is independent of jet truncation.
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
