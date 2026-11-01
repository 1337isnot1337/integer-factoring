// Complete binary compound blocks under constant coherent row/column maps.
// The determinant degree bound and DVR arithmetic are inherited unchanged.
#include "selected_gaussian_coefficients.cpp"

extern "C" int coherent_hot_compound_block(
 int modes,int count,const std::uint64_t*input,int input_words,
 const std::uint32_t*row_offsets,const std::uint32_t*row_masks,
 const std::uint32_t*column_offsets,const std::uint32_t*column_masks,
 std::int64_t*orders,std::uint64_t*coefficients){
 try{
  if(modes<1||modes>16||count<1||count>4096||input_words<1||input_words>words)return 1;
  const int n=2*modes;
  std::vector<std::vector<Poly>>one(n,std::vector<Poly>(n));
  for(int i=0;i<n;++i)for(int j=0;j<n;++j)for(int k=0;k<input_words;++k)
   one[i][j].w[k]=input[(i*n+j)*input_words+k];
  const std::uint32_t low_mask=(std::uint32_t(1)<<modes)-1;
  int charge=-1;
  std::vector<std::vector<std::vector<int>>>rr(count),cc(count);
  for(int i=0;i<count;++i){
   if(row_offsets[i]>=row_offsets[i+1]||column_offsets[i]>=column_offsets[i+1])return 2;
   for(auto a=row_offsets[i];a<row_offsets[i+1];++a){
    std::vector<int>selected;
    for(int j=0;j<n;++j)if((row_masks[a]>>j)&1)selected.push_back(j);
    if(charge<0)charge=int(selected.size());
    if(int(selected.size())!=charge)return 2;
    rr[i].push_back(std::move(selected));
   }
   for(auto a=column_offsets[i];a<column_offsets[i+1];++a){
    std::vector<int>selected;
    for(int j=0;j<n;++j)if((column_masks[a]>>j)&1)selected.push_back(j);
    if(int(selected.size())!=charge)return 2;
    cc[i].push_back(std::move(selected));
   }
  }
  f4::Matrix matrix(count,std::vector<f4::Series>(count));
  std::vector<std::array<int,4>>row_labels(count),column_labels(count);
  std::vector<double>weights(count,0.0);
  for(int i=0;i<count;++i){
   row_labels[i]={i,0,0,0};column_labels[i]={0,0,i,0};
   for(int j=0;j<count;++j){
    Poly value;
    for(const auto&r:rr[i])for(const auto&c:cc[j]){
     std::vector<std::vector<Poly>>minor(charge,std::vector<Poly>(charge));
     for(int a=0;a<charge;++a)for(int b=0;b<charge;++b)minor[a][b]=one[r[a]][c[b]];
     value^=exact_binary_determinant(std::move(minor));
    }
    int endpoint=0;
    for(int a=0;a<words;++a)endpoint^=__builtin_parityll(value.w[a]);
    if(endpoint!=int(i==j))return 3;
    // Independently evaluate the sparse-map endpoint before polynomial arithmetic.
    int expected=0;
    for(auto a=row_offsets[i];a<row_offsets[i+1];++a)
     for(auto b=column_offsets[j];b<column_offsets[j+1];++b){
      auto c=column_masks[b];
      auto swapped=(c>>modes)|((c&low_mask)<<modes);
      expected^=int(row_masks[a]==swapped);
     }
    if(expected!=int(i==j))return 4;
    matrix[i][j].bits.lo=value;
    if(coefficients)for(int a=0;a<words;++a)
     coefficients[(std::size_t(i)*count+j)*words+a]=value.w[a];
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
