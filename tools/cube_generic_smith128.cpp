// Exact truncated-series Smith orders over F_2[z]/(z^128+z^7+z^2+z+1).
// The Python wrapper proves that this field specialization preserves every
// generic determinantal valuation before dispatching a matrix here.
#include <cstdint>
#include <vector>
#include <algorithm>
#if defined(__PCLMUL__)
#include <immintrin.h>
#endif
using U=__uint128_t;
inline U gf(U a,U b){
#if defined(__PCLMUL__)
 __m128i av=_mm_set_epi64x(uint64_t(a>>64),uint64_t(a)),bv=_mm_set_epi64x(uint64_t(b>>64),uint64_t(b));
 __m128i p=_mm_clmulepi64_si128(av,bv,0),q=_mm_clmulepi64_si128(av,bv,17);
 __m128i r=_mm_xor_si128(_mm_clmulepi64_si128(av,bv,1),_mm_clmulepi64_si128(av,bv,16));
 uint64_t pw[2],qw[2],rw[2];_mm_storeu_si128((__m128i*)pw,p);_mm_storeu_si128((__m128i*)qw,q);_mm_storeu_si128((__m128i*)rw,r);
 U lo=U(pw[0])|(U(pw[1]^rw[0])<<64),hi=U(qw[0]^rw[1])|(U(qw[1])<<64);
 U spill=(hi>>127)^(hi>>126)^(hi>>121);
 return lo^hi^(hi<<1)^(hi<<2)^(hi<<7)^spill^(spill<<1)^(spill<<2)^(spill<<7);
#else
 U result=0;
 for(int i=0;i<128;++i){if(b&1)result^=a;b>>=1;bool high=bool(a>>127);a<<=1;if(high)a^=U(0x87);}
 return result;
#endif
}
inline U inverse(U a){U b=a;for(int i=1;i<127;++i)b=gf(gf(b,b),a);return gf(b,b);}
using Poly=std::vector<U>;
int val(const Poly&a){for(int i=0;i<int(a.size());++i)if(a[i])return i;return int(a.size());}
Poly product(const Poly&a,const Poly&b,int precision){
 Poly out(precision);std::vector<int>ai,bi;for(int i=0;i<int(a.size());++i)if(a[i])ai.push_back(i);for(int j=0;j<int(b.size());++j)if(b[j])bi.push_back(j);
 for(int i:ai)for(int j:bi){if(i+j>=precision)break;out[i+j]^=gf(a[i],b[j]);}return out;
}
Poly unit_inverse(const Poly&a,int precision){
 Poly out(1,inverse(a[0]));
 for(int w=1;w<precision;w*=2){int n=std::min(precision,2*w);Poly square(n);for(int i=0;i<int(out.size())&&2*i<n;++i)square[2*i]=gf(out[i],out[i]);out=product(a,square,n);}return out;
}
extern "C" int cube_gf128_smith(int n,int precision,const uint64_t*input,int64_t*out){
 try{
  if(n<1||n>64||precision<1||precision>128)return 3;
  std::vector<std::vector<Poly>>A(n,std::vector<Poly>(n,Poly(precision)));
  for(int i=0;i<n;++i)for(int j=0;j<n;++j)for(int k=0;k<precision;++k){size_t at=((size_t(i)*n+j)*precision+k)*2;A[i][j][k]=U(input[at])|(U(input[at+1])<<64);}
  for(int k=0;k<n;++k){
   int e=precision,pi=-1,pj=-1;for(int i=k;i<n;++i)for(int j=k;j<n;++j){int v=val(A[i][j]);if(v<e){e=v;pi=i;pj=j;}}
   if(pi<0)return 1;out[k]=e;std::swap(A[k],A[pi]);for(auto&r:A)std::swap(r[k],r[pj]);
   Poly unit(A[k][k].begin()+e,A[k][k].end());auto inv=unit_inverse(unit,precision-e);
   for(int j=k+1;j<n;++j)A[k][j]=product(A[k][j],inv,precision);
   for(int i=k+1;i<n;++i){if(val(A[i][k])==precision)continue;Poly factor(A[i][k].begin()+e,A[i][k].end());for(int j=k+1;j<n;++j){auto p=product(factor,A[k][j],precision);for(int z=0;z<precision;++z)A[i][j][z]^=p[z];}}
  }return 0;
 }catch(...){return 2;}
}
extern "C" void cube_gf128_product(const uint64_t*a,const uint64_t*b,uint64_t*out){U p=gf(U(a[0])|(U(a[1])<<64),U(b[0])|(U(b[1])<<64));out[0]=uint64_t(p);out[1]=uint64_t(p>>64);}
