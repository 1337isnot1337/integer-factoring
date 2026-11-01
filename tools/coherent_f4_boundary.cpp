// Exact F4 coefficient compression of homogeneous boundary retractions.
#include "coherent_hot_compounds_fast.cpp"
#include <unordered_map>

namespace coherent_boundary_f4 {
static std::unordered_map<std::uint64_t, Poly> cache;
static int active_modes = -1;
static std::vector<std::vector<Poly>> one;
unsigned multiply(unsigned a, unsigned b) {
 unsigned a0=a&1,a1=a>>1,b0=b&1,b1=b>>1;
 return ((a0*b0)^(a1*b1)) | (((a0*b1)^(a1*b0)^(a1*b1))<<1);
}
int height(std::uint32_t mask,int modes){
 int value=0;for(int j=0;j<2*modes;++j)if(mask>>j&1)value+=2*(j%modes)+1;return value;
}
Poly determinant(std::uint32_t r,std::uint32_t c,int modes,int bound){
 if(__builtin_popcount(r)!=__builtin_popcount(c))return Poly();
 std::uint64_t key=(std::uint64_t(r)<<32)|c;
 auto found=cache.find(key);if(found!=cache.end())return found->second;
 std::vector<int>rr,cc;
 for(int j=0;j<2*modes;++j){if(r>>j&1)rr.push_back(j);if(c>>j&1)cc.push_back(j);}
 const int charge=int(rr.size());
 std::vector<std::vector<Poly>>minor(charge,std::vector<Poly>(charge));
 for(int i=0;i<charge;++i)for(int j=0;j<charge;++j)minor[i][j]=one[rr[i]][cc[j]];
 Poly value=hot_determinant(std::move(minor),bound);cache.emplace(key,value);return value;
}
}

extern "C" int coherent_boundary_f4_block(
 int modes,int count,
 const std::uint32_t*row_offsets,const std::uint32_t*row_masks,const std::uint32_t*row_coeff,
 const std::uint32_t*column_offsets,const std::uint32_t*column_masks,const std::uint32_t*column_coeff,
 std::int64_t*orders){
 try{
  if(modes<1||modes>16||count<1||count>4096)return 1;
  if(coherent_boundary_f4::active_modes!=modes){
   coherent_boundary_f4::active_modes=modes;coherent_boundary_f4::cache.clear();
   int n=2*modes;coherent_boundary_f4::one.assign(n,std::vector<Poly>(n));
   for(int i=0;i<n;++i)for(int j=0;j<n;++j){
    if(i/modes==j/modes)coherent_boundary_f4::one[i][j]=Poly(3)<<((i%modes)^(j%modes));
    else if(i%modes==j%modes)coherent_boundary_f4::one[i][j]=Poly(1);
   }
  }
  std::vector<int>rh(count,-1),ch(count,-1);
  for(int i=0;i<count;++i)for(int side=0;side<2;++side){
   const auto offsets=side?column_offsets:row_offsets;
   const auto masks=side?column_masks:row_masks;
   const auto coeff=side?column_coeff:row_coeff;
   int&h=side?ch[i]:rh[i];
   for(auto j=offsets[i];j<offsets[i+1];++j){
    int value=coherent_boundary_f4::height(masks[j],modes);
    if(h<0)h=value;if(h!=value||coeff[j]<1||coeff[j]>3)return 2;
   }
  }
  f4::Matrix matrix(count,std::vector<f4::Series>(count));
  std::vector<std::array<int,4>>row_labels(count),column_labels(count);
  std::vector<double>weights(count,0.0);
  for(int i=0;i<count;++i){
   row_labels[i]={i,0,0,0};column_labels[i]={0,0,i,0};
   for(int j=0;j<count;++j){
    if((rh[i]+ch[j])&1)return 2;const int exponent=(rh[i]+ch[j])/2;
    Poly lo,hi;
    for(auto a=row_offsets[i];a<row_offsets[i+1];++a)
     for(auto b=column_offsets[j];b<column_offsets[j+1];++b){
      unsigned coefficient=coherent_boundary_f4::multiply(row_coeff[a],column_coeff[b]);
      Poly value=coherent_boundary_f4::determinant(row_masks[a],column_masks[b],modes,exponent);
      if(coefficient&1)lo^=value;if(coefficient&2)hi^=value;
     }
    unsigned endpoint=0;
    for(int a=0;a<words;++a)endpoint^=__builtin_parityll(lo.w[a])|(__builtin_parityll(hi.w[a])<<1);
    if(endpoint!=unsigned(i==j))return 3;
    matrix[i][j].bits.lo=lo;matrix[i][j].bits.hi=hi;matrix[i][j].base=exponent;
   }
  }
  auto pivots=f4::profile_block(nullptr,matrix,row_labels,column_labels,false,1.0,count,&weights,&weights);
  if(int(pivots.size())!=count)return 6;
  for(int i=0;i<count;++i){if(pivots[i].value!=double(pivots[i].base)||pivots[i].base<0)return 7;orders[i]=pivots[i].base;}
  return 0;
 }catch(...){return 8;}
}
