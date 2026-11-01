// Exact selected compound blocks and their complete ordinary Smith profiles.
// Existing determinant and DVR elimination kernels are included unchanged.
#include "selected_gaussian_coefficients.cpp"

extern "C" int hot_code_compound_block(
 int modes,int count,const std::uint64_t*input,int input_words,
 const std::uint32_t*row_masks,const std::uint32_t*column_masks,
 std::int64_t*orders,std::uint64_t*coefficients){
 try{
  if(modes<1||modes>16||count<1||count>4096||input_words<1||input_words>words)return 1;
  const int n=2*modes;
  std::vector<std::vector<Poly>>one(n,std::vector<Poly>(n));
  for(int i=0;i<n;++i)for(int j=0;j<n;++j)for(int k=0;k<input_words;++k)
   one[i][j].w[k]=input[(i*n+j)*input_words+k];
  std::vector<std::vector<int>>rr(count),cc(count);
  for(int i=0;i<count;++i){
   for(int j=0;j<n;++j){
    if((row_masks[i]>>j)&1)rr[i].push_back(j);
    if((column_masks[i]>>j)&1)cc[i].push_back(j);
   }
   if(rr[i].size()!=rr[0].size()||cc[i].size()!=rr[0].size())return 2;
  }
  f4::Matrix matrix(count,std::vector<f4::Series>(count));
  std::vector<std::array<int,4>>row_labels(count),column_labels(count);
  std::vector<double>weights(count,0.0);
  const std::uint32_t mask=(std::uint32_t(1)<<modes)-1;
  std::vector<bool>matched(count,false);
  for(int i=0;i<count;++i){
   row_labels[i]={i,0,0,0};column_labels[i]={0,0,i,0};
   for(int j=0;j<count;++j){
    const int k=int(rr[i].size());
    std::vector<std::vector<Poly>>minor(k,std::vector<Poly>(k));
    for(int a=0;a<k;++a)for(int b=0;b<k;++b)minor[a][b]=one[rr[i][a]][cc[j][b]];
    Poly value=exact_binary_determinant(std::move(minor));
    int endpoint=0;for(int a=0;a<words;++a)endpoint^=__builtin_parityll(value.w[a]);
    const std::uint32_t swapped=(column_masks[j]>>modes)|((column_masks[j]&mask)<<modes);
    const bool expected=row_masks[i]==swapped;
    if(endpoint!=int(expected))return 3;
    if(expected){if(matched[i])return 4;matched[i]=true;}
    matrix[i][j].bits.lo=value;
    if(coefficients)for(int a=0;a<words;++a)
     coefficients[(std::size_t(i)*count+j)*words+a]=value.w[a];
   }
  }
  for(bool found:matched)if(!found)return 5;
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
