// Ordinary F8 polynomial determinant; the degree guard prevents truncation.
#include "coherent_hot_compounds_f8.cpp"
f8::Poly scalar_poly(const f8::Poly&a,unsigned c){
 f8::Poly out,term=a;
 while(c){if(c&1)out^=term;c>>=1;f8::Poly next;next.lo=term.top;next.hi=term.lo^term.top;next.top=term.hi;term=next;}
 return out;
}
unsigned coefficient(const f8::Poly&a,int d){return unsigned((a.lo.w[d/64]>>(d%64))&1)|unsigned(((a.hi.w[d/64]>>(d%64))&1)<<1)|unsigned(((a.top.w[d/64]>>(d%64))&1)<<2);}
f8::Poly exact_quotient(f8::Poly a,const f8::Poly&b){
 int bd=f8::degree(b);if(bd<0)throw std::runtime_error("zero denominator");unsigned inverse[8]={0,1,5,6,7,2,3,4};unsigned inv=inverse[coefficient(b,bd)];auto unit=scalar_poly(b,inv);f8::Poly out;
 for(int ad=f8::degree(a);ad>=bd;ad=f8::degree(a)){
  int shift=ad-bd;unsigned c=coefficient(a,ad);auto term=scalar_poly(unit,c);a^=term<<shift;
  unsigned q=scalar_poly(f8::Poly(inv),c).constant();
  if(q&1)out.lo.w[shift/64]^=std::uint64_t(1)<<(shift%64);
  if(q&2)out.hi.w[shift/64]^=std::uint64_t(1)<<(shift%64);
  if(q&4)out.top.w[shift/64]^=std::uint64_t(1)<<(shift%64);
 }
 if(a)throw std::runtime_error("nonexact Bareiss division");return out;
}
extern "C" int determinant_f8(int n,int D,const std::uint64_t*data,int nw,std::uint64_t*out){
 try{
  if(n<1||n>128||D<0||D>=precision||2*n*D>=precision||nw<1||nw>words)return 1;
  std::vector<std::vector<f8::Poly>>a(n,std::vector<f8::Poly>(n));
  for(int i=0;i<n;++i)for(int j=0;j<n;++j)for(int w=0;w<nw;++w){auto k=std::size_t(i)*n+j;a[i][j].lo.w[w]=data[(3*k)*nw+w];a[i][j].hi.w[w]=data[(3*k+1)*nw+w];a[i][j].top.w[w]=data[(3*k+2)*nw+w];}
  for(const auto&row:a)for(const auto&v:row)if(f8::degree(v)>D)return 1;
  f8::Poly previous(1),pivot;
  for(int k=0;k<n;++k){
   int pi=-1,pj=-1,best=precision;
   for(int i=k;i<n;++i)for(int j=k;j<n;++j)if(a[i][j]){int d=f8::degree(a[i][j]);if(d<best){best=d;pi=i;pj=j;}}
   if(pi<0)return 3;std::swap(a[k],a[pi]);for(auto&row:a)std::swap(row[k],row[pj]);pivot=a[k][k];
   for(int i=k+1;i<n;++i)for(int j=k+1;j<n;++j){auto numerator=f8::multiply(pivot,a[i][j])^f8::multiply(a[i][k],a[k][j]);a[i][j]=exact_quotient(numerator,previous);if(f8::degree(a[i][j])>(k+2)*D)throw std::runtime_error("minor degree");}
   previous=pivot;
  }
  for(int w=0;w<words;++w){out[w]=pivot.lo.w[w];out[words+w]=pivot.hi.w[w];out[2*words+w]=pivot.top.w[w];}return 0;
 }catch(...){return 2;}
}
