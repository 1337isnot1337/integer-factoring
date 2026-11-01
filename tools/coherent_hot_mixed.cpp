// Height-homogeneous coherent maps may mix occupation numbers.
// Only minors with equal row/column occupations contribute to a compound entry.
#include "coherent_hot_compounds_fast.cpp"

extern "C" int coherent_hot_mixed_block(
 int modes,int count,
 const std::uint32_t*row_offsets,const std::uint32_t*row_masks,
 const std::uint32_t*column_offsets,const std::uint32_t*column_masks,
 std::int64_t*orders,std::uint64_t*coefficients){
 try{
  if(modes<1||modes>16||count<1||count>4096)return 1;
  const int n=2*modes;
  std::vector<std::vector<Poly>>one(n,std::vector<Poly>(n));
  for(int i=0;i<n;++i)for(int j=0;j<n;++j){
   if(i/modes==j/modes)one[i][j]=Poly(3)<<((i%modes)^(j%modes));
   else if(i%modes==j%modes)one[i][j]=Poly(1);
  }
  const std::uint32_t low_mask=(std::uint32_t(1)<<modes)-1;
  std::vector<std::vector<std::vector<int>>>rr(count),cc(count);
  std::vector<int>row_height(count,-1),column_height(count,-1);
  for(int i=0;i<count;++i){
   if(row_offsets[i]>=row_offsets[i+1]||column_offsets[i]>=column_offsets[i+1])return 2;
   for(int side=0;side<2;++side){
    auto offsets=side?column_offsets:row_offsets;
    auto masks=side?column_masks:row_masks;
    auto&vectors=side?cc[i]:rr[i];
    auto&hs=side?column_height[i]:row_height[i];
    for(auto a=offsets[i];a<offsets[i+1];++a){
     std::vector<int>selected;int h=0;
     for(int j=0;j<n;++j)if((masks[a]>>j)&1){selected.push_back(j);h+=2*(j%modes)+1;}
     if(hs<0)hs=h;
     if(hs!=h)return 2;
     vectors.push_back(std::move(selected));
    }
   }
  }
  f4::Matrix matrix(count,std::vector<f4::Series>(count));
  std::vector<std::array<int,4>>row_labels(count),column_labels(count);
  std::vector<double>weights(count,0.0);
  for(int i=0;i<count;++i){
   row_labels[i]={i,0,0,0};column_labels[i]={0,0,i,0};
   for(int j=0;j<count;++j){
    int height_sum=row_height[i]+column_height[j];
    if(height_sum&1)return 2;
    const int exponent=height_sum/2;
    Poly value;
    for(const auto&r:rr[i])for(const auto&c:cc[j]){
     if(r.size()!=c.size())continue;
     const int charge=int(r.size());
     std::vector<std::vector<Poly>>minor(charge,std::vector<Poly>(charge));
     for(int a=0;a<charge;++a)for(int b=0;b<charge;++b)minor[a][b]=one[r[a]][c[b]];
     value^=hot_determinant(std::move(minor),exponent);
    }
    int endpoint=0;
    for(int a=0;a<words;++a)endpoint^=__builtin_parityll(value.w[a]);
    if(endpoint!=int(i==j))return 3;
    int expected=0;
    for(auto a=row_offsets[i];a<row_offsets[i+1];++a)
     for(auto b=column_offsets[j];b<column_offsets[j+1];++b){
      auto c=column_masks[b];auto swapped=(c>>modes)|((c&low_mask)<<modes);
      expected^=int(row_masks[a]==swapped);
     }
    if(expected!=int(i==j))return 4;
    matrix[i][j].bits.lo=value;matrix[i][j].base=exponent;
    if(coefficients){
     if(degree(value)+exponent>=precision)return 9;
     auto full=value<<exponent;
     for(int a=0;a<words;++a)coefficients[(std::size_t(i)*count+j)*words+a]=full.w[a];
    }
   }
  }
  auto pivots=f4::profile_block(nullptr,matrix,row_labels,column_labels,false,1.0,
                               count,&weights,&weights);
  if(int(pivots.size())!=count)return 6;
  for(int i=0;i<count;++i){
   const auto&p=pivots[i];
   if(p.value!=double(p.base)||p.base<0)return 7;
   orders[i]=p.base;
  }
  return 0;
 }catch(...){return 8;}
}
