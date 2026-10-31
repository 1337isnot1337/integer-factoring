// Exact selected exterior coefficients over F2[w].
// The determinant degree bound is checked before using any truncated series.
#include "native_rectangular_smith_f4.cpp"

Poly exact_binary_determinant(std::vector<std::vector<Poly>> a){
 const int n=int(a.size());
 if(n==0)return Poly(1);
 int bound=0;
 for(const auto&row:a){int d=-1;for(const auto&p:row)d=std::max(d,degree(p));if(d<0)return Poly();bound+=d;}
 if(bound>=precision)throw std::runtime_error("determinant degree exceeds precision");
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
  if(unit!=Poly(1))for(int width=1;width<precision;width*=2)inverse=multiply(unit,multiply(inverse,inverse));
  for(int i=k+1;i<n;++i)if(a[i][k]){
   Poly factor=multiply(a[i][k]>>v,inverse);
   for(int j=k+1;j<n;++j)if(a[k][j])a[i][j]^=multiply(factor,a[k][j]);
  }
 }
 if(degree(determinant)>bound)throw std::runtime_error("determinant exceeds its proved degree bound");
 return determinant;
}

extern "C" int selected_binary_exterior(int modes,int A,int B,const std::uint64_t*input,int input_words,
 const std::uint32_t*xletters,const std::uint32_t*yletters,std::uint64_t*output,int check_endpoint){
 try{
  if(modes<1||modes>16||A<1||B<1||A>64||B>64||input_words<1||input_words>words)return 1;
  int n=2*modes;std::vector<std::vector<Poly>> one(n,std::vector<Poly>(n));
  for(int i=0;i<n;++i)for(int j=0;j<n;++j)for(int k=0;k<input_words;++k)one[i][j].w[k]=input[(i*n+j)*input_words+k];
  for(int v=0;v<A;++v)for(int u=0;u<B;++u)for(int y=0;y<B;++y)for(int x=0;x<A;++x){
   std::uint32_t rr=xletters[v]|(yletters[y]<<modes),cc=yletters[u]|(xletters[x]<<modes);Poly value;
   if(__builtin_popcount(rr)==__builtin_popcount(cc)){
    std::vector<int>rows,columns;
    for(int i=0;i<n;++i){if((rr>>i)&1)rows.push_back(i);if((cc>>i)&1)columns.push_back(i);}
    int count=int(rows.size());std::vector<std::vector<Poly>> a(count,std::vector<Poly>(count));
    for(int i=0;i<count;++i)for(int j=0;j<count;++j)a[i][j]=one[rows[i]][columns[j]];
    value=exact_binary_determinant(std::move(a));
   }
   std::size_t index=x+std::size_t(A)*(y+std::size_t(B)*(u+std::size_t(B)*v));
   int endpoint=0;for(int k=0;k<words;++k){output[index*words+k]=value.w[k];endpoint^=__builtin_parityll(value.w[k]);}
   if(check_endpoint&&endpoint!=int(x==v&&y==u))return 3;
  }
  return 0;
 }catch(...){return 2;}
}
