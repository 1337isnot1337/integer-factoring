// Generic Laurent-polynomial matrix interface for the existing F4 engine.
// The caller supplies the generic rank; truncated zeros do not prove it.
#include "native_rectangular_smith_f4.cpp"

extern "C" int matrix_f4_profile(int n,int rank,const std::uint64_t*coefficients,
 int words_per_value,const std::int64_t*shifts,const double*row_weights,
 const double*column_weights,int infinity,double scale,double*orders,
 std::int64_t*prefix_values,std::uint32_t*chosen_rows,std::uint32_t*chosen_columns){
 try{
  if(n<1||n>1024||rank<1||rank>n||words_per_value<1||words_per_value>words||scale<=0)return 1;
  f4::Matrix matrix(n,std::vector<f4::Series>(n));
  std::vector<std::array<int,4>>rows(n),cols(n);
  std::vector<double>rw(row_weights,row_weights+n),cw(column_weights,column_weights+n);
  for(int i=0;i<n;++i){rows[i]={i,0,0,0};cols[i]={0,0,i,0};}
  for(int i=0;i<n;++i)for(int j=0;j<n;++j){
   int index=i*n+j;f4::Poly p;
   for(int k=0;k<words_per_value;++k){
    p.lo.w[k]=coefficients[(2*index)*words_per_value+k];
    p.hi.w[k]=coefficients[(2*index+1)*words_per_value+k];
   }
   if(shifts[index]<-1000000||shifts[index]>1000000)return 2;
   matrix[i][j]={p,int(shifts[index])};
  }
  auto result=f4::profile_block(nullptr,matrix,rows,cols,bool(infinity),scale,rank,&rw,&cw);
  if(int(result.size())!=rank)return 3;
  std::int64_t cumulative=0;
  for(int i=0;i<rank;++i){
   auto p=result[i];orders[i]=p.value;cumulative+=p.base;prefix_values[i]=cumulative;
   chosen_rows[i]=p.labels[0];chosen_columns[i]=p.labels[2];
  }
  return 0;
 }catch(...){return 4;}
}
