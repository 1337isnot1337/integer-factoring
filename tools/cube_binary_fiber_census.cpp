#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>
#if defined(__PCLMUL__)
#include <wmmintrin.h>
#endif
using U=std::uint64_t;using V=__uint128_t;
V multiply(V a,V b){
#if defined(__PCLMUL__)
__m128i x=_mm_set_epi64x(U(a>>64),U(a)),y=_mm_set_epi64x(U(b>>64),U(b));auto p=_mm_clmulepi64_si128(x,y,0);if(U(a>>64)||U(b>>64)){auto q=_mm_xor_si128(_mm_clmulepi64_si128(x,y,1),_mm_clmulepi64_si128(x,y,16));p=_mm_xor_si128(p,_mm_slli_si128(q,8));}alignas(16) U out[2];_mm_store_si128((__m128i*)out,p);return V(out[0])|(V(out[1])<<64);
#else
V out=0;while(b){if(b&1)out^=a;b>>=1;a<<=1;}return out;
#endif
}
int valuation(V a){return U(a)?__builtin_ctzll(U(a)):64+__builtin_ctzll(U(a>>64));}
std::array<int,8> profile(const std::uint32_t*input){V a[8][8];for(int i=0;i<8;++i)for(int j=0;j<8;++j)a[i][j]=input[8*i+j];std::array<int,8> out{};for(int k=0;k<8;++k){int best=128,pi=-1,pj=-1;for(int i=k;i<8;++i)for(int j=k;j<8;++j)if(a[i][j]){int v=valuation(a[i][j]);if(v<best){best=v;pi=i;pj=j;}}if(pi<0)throw std::runtime_error("missing promised pivot");for(int j=0;j<8;++j)std::swap(a[k][j],a[pi][j]);for(int i=0;i<8;++i)std::swap(a[i][k],a[i][pj]);out[k]=best;V unit=a[k][k]>>best,inv=1;if(unit!=1)for(int s=1;s<128;s*=2)inv=multiply(unit,multiply(inv,inv));if(multiply(unit,inv)!=1)throw std::runtime_error("unit inverse failed");for(int i=k+1;i<8;++i)if(a[i][k]){V q=multiply(a[i][k]>>best,inv);for(int j=k+1;j<8;++j)if(a[k][j])a[i][j]^=multiply(q,a[k][j]);a[i][k]=0;}}return out;}
struct Count{U count=0,representative=0;};
int main(int argc,char**argv){try{if(argc<3)return 2;std::ifstream in(argv[1]);int d;in>>d;std::vector<std::uint32_t> current(384);std::vector<std::vector<std::uint32_t>> dirs(d,std::vector<std::uint32_t>(384));for(auto&v:current)in>>v;for(auto&row:dirs)for(auto&v:row)in>>v;U limit=U(1)<<d;if(argc>3)limit=std::min(limit,U(std::stoull(argv[3])));double seconds=argc>4?std::stod(argv[4]):1800;auto start=std::chrono::steady_clock::now();std::array<int,16> base={0,2,4,6,8,10,10,12,12,14,14,16,18,20,22,24};std::unordered_map<std::string,Count> classes;std::vector<std::pair<U,std::string>> controls;U done=0,dominated=0;bool complete=true;
for(U step=0;step<limit;++step){if(step){int bit=__builtin_ctzll(step);for(int i=0;i<384;++i)current[i]^=dirs[bit][i];}std::string key(48,'\0');bool dom=true;for(int axis=0;axis<3;++axis){std::array<int,16> orders{};for(int parity=0;parity<2;++parity){auto p=profile(current.data()+64*(2*axis+parity));std::copy(p.begin(),p.end(),orders.begin()+8*parity);}std::sort(orders.begin(),orders.end());for(int i=0;i<16;++i){key[16*axis+i]=char(orders[i]);if(orders[i]>base[i])dom=false;}}if(step<32||step%65537==0)controls.push_back({step^(step>>1),key});bool major=true;for(int axis=0;axis<3;++axis){int sum=0,refsum=0;for(int i=0;i<16;++i){sum+=int((unsigned char)key[16*axis+i]);refsum+=base[i];if(sum>refsum)major=false;}}if(!major){auto&entry=classes[key];if(!entry.count)entry.representative=step^(step>>1);++entry.count;}dominated+=major;done=step+1;if(done%65536==0){double elapsed=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();std::cerr<<"cases "<<done<<" classes "<<classes.size()<<" dominated "<<dominated<<" seconds "<<elapsed<<'\n';if(elapsed>seconds){complete=false;break;}}}
std::ofstream out(argv[2]);out<<"{\"dimension\":"<<d<<",\"cases\":"<<done<<",\"cutwise_weakly_majorized\":"<<dominated<<",\"completed_requested_cases\":"<<(complete?"true":"false")<<",\"complete_binary_space\":"<<(done==(U(1)<<d)?"true":"false")<<",\"exception_classes\":[";bool first=true;for(const auto&[key,val]:classes){if(!first)out<<',';first=false;out<<"{\"count\":"<<val.count<<",\"representative\":"<<val.representative<<",\"profiles\":[";for(int axis=0;axis<3;++axis){if(axis)out<<',';out<<'[';for(int i=0;i<16;++i)out<<(i?",":"")<<int((unsigned char)key[16*axis+i]);out<<']';}out<<"]}";}out<<"],\"profile_controls\":[";for(std::size_t k=0;k<controls.size();++k){if(k)out<<',';out<<"{\"assignment\":"<<controls[k].first<<",\"orders\":[";for(int i=0;i<48;++i)out<<(i?",":"")<<int((unsigned char)controls[k].second[i]);out<<"]}";}out<<"],\"elapsed_seconds\":"<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<"}\n";std::cerr<<"finished "<<done<<" classes "<<classes.size()<<'\n';return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 2;}}
