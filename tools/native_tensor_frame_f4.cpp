// Exact polynomial matrix sandwiches modulo z^precision.
#include "native_rectangular_smith_f4.cpp"

extern "C" int tensor_frame_f4(int n, const std::uint64_t* coefficients,
 int input_words, std::uint64_t* result) {
 try {
  if(n<1||n>32||input_words<1||input_words>words)return 1;
  std::vector<f4::Poly> data(3*n*n), intermediate(n*n), output(n*n);
  for(int i=0;i<3*n*n;++i)for(int k=0;k<input_words;++k){
   data[i].lo.w[k]=coefficients[(2*i)*input_words+k];
   data[i].hi.w[k]=coefficients[(2*i+1)*input_words+k];
  }
  for(int i=0;i<n;++i)for(int k=0;k<n;++k)if(data[i*n+k])
   for(int j=0;j<n;++j)if(data[n*n+k*n+j])
    intermediate[i*n+j]^=f4::multiply(data[i*n+k],data[n*n+k*n+j]);
  for(int i=0;i<n;++i)for(int k=0;k<n;++k)if(intermediate[i*n+k])
   for(int j=0;j<n;++j)if(data[2*n*n+k*n+j])
    output[i*n+j]^=f4::multiply(intermediate[i*n+k],data[2*n*n+k*n+j]);
  for(int i=0;i<n*n;++i)for(int k=0;k<words;++k){
   result[(2*i)*words+k]=output[i].lo.w[k];
   result[(2*i+1)*words+k]=output[i].hi.w[k];
  }
  return 0;
 }catch(...){return 2;}
}
