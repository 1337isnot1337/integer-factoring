#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#if defined(__PCLMUL__)
#include <wmmintrin.h>
#endif
using U=std::uint64_t;
using P=std::vector<U>;
void trim(P &a){while(!a.empty()&&!a.back())a.pop_back();}
int deg(const P&a){return a.empty()?-1:int(64*(a.size()-1)+63-__builtin_clzll(a.back()));}
P parse(std::string s){P a;for(int e=int(s.size());e>0;e-=16){int b=std::max(0,e-16);a.push_back(std::stoull(s.substr(b,e-b),nullptr,16));}trim(a);return a;}
std::string hex(const P&a){if(a.empty())return "0";std::ostringstream s;s<<std::hex<<a.back();for(int i=int(a.size())-2;i>=0;--i)s<<std::setfill('0')<<std::setw(16)<<a[i];return s.str();}
void xorshift(P&a,const P&b,int shift){if(b.empty())return;int word=shift/64,bit=shift%64;a.resize(std::max(a.size(),b.size()+word+(bit!=0)),0);for(std::size_t i=0;i<b.size();++i){a[i+word]^=b[i]<<bit;if(bit)a[i+word+1]^=b[i]>>(64-bit);}trim(a);}
P mul(const P&a,const P&b){if(a.empty()||b.empty())return {};P c(a.size()+b.size(),0);for(std::size_t i=0;i<a.size();++i)if(a[i])for(std::size_t j=0;j<b.size();++j)if(b[j]){
#if defined(__PCLMUL__)
__m128i u=_mm_set_epi64x(0,(long long)a[i]),v=_mm_set_epi64x(0,(long long)b[j]);__m128i w=_mm_clmulepi64_si128(u,v,0);alignas(16) U z[2];_mm_store_si128((__m128i*)z,w);c[i+j]^=z[0];c[i+j+1]^=z[1];
#else
U value=a[i];while(value){int bit=__builtin_ctzll(value);value&=value-1;c[i+j]^=b[j]<<bit;if(bit)c[i+j+1]^=b[j]>>(64-bit);}
#endif
}trim(c);return c;}
P exactdiv(P a,const P&b){if(b.empty())throw std::runtime_error("zero divisor");if(b.size()==1&&b[0]==1)return a;P q;int db=deg(b);while(deg(a)>=db){int s=deg(a)-db;if(q.size()<=std::size_t(s/64))q.resize(s/64+1,0);q[s/64]^=U(1)<<(s%64);xorshift(a,b,s);}if(!a.empty())throw std::runtime_error("nonexact Bareiss division");trim(q);return q;}
int vx(const P&a){int v=256;for(std::size_t i=0;i<a.size();++i)if(a[i]){int at=int(i%4)*64+__builtin_ctzll(a[i]);v=std::min(v,at);if(!v)break;}return v;}
P shiftx(const P&a,int v){if(v==0)return a;if(vx(a)<v)throw std::runtime_error("nonregular normalized frame");int word=v/64,bit=v%64;P b(a.size()>std::size_t(word)?a.size()-word:0,0);for(std::size_t j=0;j<b.size();++j){b[j]=a[j+word]>>bit;if(bit&&j+word+1<a.size())b[j]^=a[j+word+1]<<(64-bit);}trim(b);return b;}
int main(){try{int n;std::cin>>n;if(n<1||n>8)throw std::runtime_error("order outside 1..8");std::vector<std::vector<P>> a(n,std::vector<P>(n)),f(n,std::vector<P>(n));std::string s;for(int i=0;i<n;++i){f[i][i]={1};for(int j=0;j<n;++j){std::cin>>s;a[i][j]=parse(s);}}std::vector<bool> ir(n),jc(n);P previous={1};std::vector<int> orders,rows,cols,prefixes;std::vector<std::vector<P>> frames;
for(int k=0;k<n;++k){int bi=-1,bj=-1,bv=100000,pv=vx(previous);for(int i=0;i<n;++i)if(!ir[i])for(int j=0;j<n;++j)if(!jc[j]&&!a[i][j].empty()){int v=vx(a[i][j]);if(v<bv){bv=v;bi=i;bj=j;}}if(bi<0)throw std::runtime_error("generic rank deficit");P pivot=a[bi][bj];orders.push_back(bv-pv);prefixes.push_back(bv);rows.push_back(bi);cols.push_back(bj);std::vector<P> row;for(auto &value:f[bi])row.push_back(shiftx(value,pv));frames.push_back(row);ir[bi]=true;jc[bj]=true;
if(k+1<n){auto next=a;auto nf=f;for(int i=0;i<n;++i)if(!ir[i]){for(int j=0;j<n;++j)if(!jc[j]){P value=mul(pivot,a[i][j]);P other=mul(a[i][bj],a[bi][j]);xorshift(value,other,0);next[i][j]=exactdiv(std::move(value),previous);}for(int j=0;j<n;++j){P value=mul(pivot,f[i][j]);P other=mul(a[i][bj],f[bi][j]);xorshift(value,other,0);nf[i][j]=exactdiv(std::move(value),previous);}}a.swap(next);f.swap(nf);}previous=std::move(pivot);}
std::cout<<"{\"orders\":[";for(int i=0;i<n;++i)std::cout<<(i?",":"")<<orders[i];std::cout<<"],\"rows\":[";for(int i=0;i<n;++i)std::cout<<(i?",":"")<<rows[i];std::cout<<"],\"columns\":[";for(int i=0;i<n;++i)std::cout<<(i?",":"")<<cols[i];std::cout<<"],\"prefixes\":[";for(int i=0;i<n;++i)std::cout<<(i?",":"")<<prefixes[i];std::cout<<"],\"frames\":[";for(int i=0;i<n;++i){if(i)std::cout<<',';std::cout<<'[';for(int j=0;j<n;++j)std::cout<<(j?",":"")<<'\"'<<hex(frames[i][j])<<'\"';std::cout<<']';}std::cout<<"]}\n";return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 2;}}
