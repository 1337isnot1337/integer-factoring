// Exhaust all normalized product densities over F2 for one explicit change of basis.
#include <array>
#include <bit>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <cassert>
#include <random>
#include <vector>
using namespace std;
using V=array<uint64_t,4>;
constexpr int n=16;
array<unsigned,n> W={16669,37713,29213,55025,12094,39882,4857,53562,35306,23446,20077,45287,25564,14040,53253,54155};
array<unsigned,n> inverse(array<unsigned,n>M){array<unsigned,n>I;for(int i=0;i<n;i++)I[i]=1u<<i;for(int j=0;j<n;j++){int p=j;while(p<n&&!((M[p]>>j)&1))p++;if(p==n)abort();swap(M[p],M[j]);swap(I[p],I[j]);for(int i=0;i<n;i++)if(i!=j&&((M[i]>>j)&1)){M[i]^=M[j];I[i]^=I[j];}}return I;}
void xx(V&A,const V&B){for(int k=0;k<4;k++)A[k]^=B[k];}
bool product_rank_one(const V&A){
 unsigned v=0;for(int k=0;k<4;k++){uint64_t z=A[k];for(int j=0;j<4;j++){unsigned r=z&65535;z>>=16;if(r){if(v&&v!=r)return false;v=r;}}}return v!=0;
}
array<unsigned,16> dense_product(const array<unsigned,16>&A,const array<unsigned,16>&B){array<unsigned,16>C{};for(int i=0;i<16;i++)for(int j=0;j<16;j++)if((A[i]>>j)&1)C[i]^=B[j];return C;}
V ordinary_realignment(unsigned a,unsigned b,const array<unsigned,16>&I){
 array<unsigned,16>K{};for(int i=0;i<4;i++)for(int j=0;j<4;j++)for(int k=0;k<4;k++)for(int l=0;l<4;l++)if(((a>>(4*i+j))&1)&&((b>>(4*k+l))&1))K[4*i+k]|=1u<<(4*j+l);
 auto C=dense_product(dense_product(I,K),W);V R{};for(int i=0;i<4;i++)for(int j=0;j<4;j++)for(int k=0;k<4;k++)for(int l=0;l<4;l++)if((C[4*i+k]>>(4*j+l))&1){int bit=(4*i+j)*16+4*k+l;R[bit/64]|=1ull<<(bit%64);}return R;
}
unsigned trace_one(unsigned a){unsigned t=0;for(int i=0;i<4;i++)t^=(a>>(5*i))&1;return a^(t?0u:1u);}
int main(int argc,char**argv){
 int start=argc>1?stoi(argv[1]):0,end=argc>2?stoi(argv[2]):32768;string path=argc>3?argv[3]:"density16.json";
 assert(0<=start && start<=end && end<=32768);
 auto I=inverse(W);array<array<V,16>,16>T{};
 for(int a=0;a<16;a++)for(int b=0;b<16;b++){
  int out=(a/4)*4+(b/4),in=(a%4)*4+(b%4);
  for(int i=0;i<4;i++)for(int j=0;j<4;j++)for(int k=0;k<4;k++)for(int l=0;l<4;l++)
   if(((I[i*4+k]>>out)&1)&&((W[in]>>(j*4+l))&1)){int bit=(i*4+j)*16+k*4+l;T[a][b][bit/64]^=1ull<<(bit%64);}
 }
 // Independent ordinary matrix multiplication and realignment controls.
 mt19937 rng(391574);uint64_t direct_checks=0;
 for(int trial=0;trial<64;trial++){
  unsigned a=trace_one(rng()&65535),b=trace_one(rng()&65535);V packed{};
  for(int i=0;i<16;i++)for(int j=0;j<16;j++)if(((a>>i)&1)&&((b>>j)&1))xx(packed,T[i][j]);
  assert(packed==ordinary_realignment(a,b,I));direct_checks++;
 }
 V zero{};assert(!product_rank_one(zero));
 V one{};one[0]=1;assert(product_rank_one(one));
 V two{};two[0]=1ull|(1ull<<17);assert(!product_rank_one(two));
 for(int trial=0;trial<64;trial++){
  unsigned a=1+(rng()%65535),b=1+(rng()%65535);V R{};
  for(int i=0;i<16;i++)if((a>>i)&1)R[i/4]|=(uint64_t)b<<(16*(i%4));
  assert(product_rank_one(R));
 }
 array<array<V,16>,15>changes{};
 for(int a=1;a<16;a++)for(int b=0;b<16;b++){changes[a-1][b]=T[a][b];if(a%5==0)xx(changes[a-1][b],T[0][b]);}
 array<V,16>X=T[0];unsigned amask=1,g=start^(start>>1);
 for(int a=0;a<15;a++)if((g>>a)&1){amask^=1u<<(a+1);if((a+1)%5==0)amask^=1;for(int b=0;b<16;b++)xx(X[b],changes[a][b]);}
 uint64_t checked=0,found=0;auto began=chrono::steady_clock::now();vector<array<unsigned,2>>witnesses;
 for(int ai=start;ai<end;ai++){
  array<V,15>Y;for(int b=1;b<16;b++){Y[b-1]=X[b];if(b%5==0)xx(Y[b-1],X[0]);}
  V P=X[0];unsigned bmask=1;
  for(int bi=0;bi<32768;bi++){
   checked++;
   if(ai%4096==0 && bi%8192==0){assert(P==ordinary_realignment(amask,bmask,I));direct_checks++;}
   if(product_rank_one(P)){found++;if(witnesses.size()<32)witnesses.push_back({amask,bmask});}
   if(bi+1<32768){int j=countr_zero((unsigned)(bi+1));xx(P,Y[j]);bmask^=1u<<(j+1);if((j+1)%5==0)bmask^=1;}
  }
  if(ai+1<end){int j=countr_zero((unsigned)(ai+1));for(int b=0;b<16;b++)xx(X[b],changes[j][b]);amask^=1u<<(j+1);if((j+1)%5==0)amask^=1;}
  if((ai+1)%1024==0)cout<<"rows="<<ai+1<<" checked="<<checked<<" found="<<found<<" seconds="<<chrono::duration<double>(chrono::steady_clock::now()-began).count()<<'\n'<<flush;
 }
 ofstream f(path);f<<"{\"field\":2,\"start\":"<<start<<",\"end\":"<<end<<",\"checked\":"<<checked<<",\"product_densities\":"<<found<<",\"ordinary_conjugation_checks\":"<<direct_checks<<",\"realignment_rank_controls\":67,\"seconds\":"<<chrono::duration<double>(chrono::steady_clock::now()-began).count()<<",\"witnesses\":[";for(size_t i=0;i<witnesses.size();i++){if(i)f<<',';f<<'['<<witnesses[i][0]<<','<<witnesses[i][1]<<']';}f<<"],\"matrix_rows\":[";for(int i=0;i<16;i++){if(i)f<<',';f<<W[i];}f<<"]}\n";
}
