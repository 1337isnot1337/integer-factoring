// Hot-kernel specialization with an exact height-based degree bound.
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

extern "C" int coherent_hot_compound_fast(
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
  f4::Matrix matrix(count,std::vector<f4::Series>(count));
  std::vector<std::array<int,4>>row_labels(count),column_labels(count);
  std::vector<double>weights(count,0.0);
  for(int i=0;i<count;++i){
   row_labels[i]={i,0,0,0};column_labels[i]={0,0,i,0};
   for(int j=0;j<count;++j){
    const int exponent=row_index_sum[i]+column_index_sum[j]+charge;
    Poly value;
    for(const auto&r:rr[i])for(const auto&c:cc[j]){
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
