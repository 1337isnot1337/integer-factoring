#include <array>
#include <bit>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <vector>
using W=uint64_t;
using Q=unsigned __int128;
using Perm=std::array<unsigned,64>;
using Mat=std::array<W,64>;
int rank64(Mat a,int cap=64){
    Mat piv{};int r=0;
    for(W v:a)while(v){int k=63-std::countl_zero(v);if(piv[k])v^=piv[k];else{piv[k]=v;if(++r>cap)return r;break;}}
    return r;
}
int rank128(std::array<Q,128>a){
    std::array<Q,128>piv{};int r=0;
    for(Q v:a)while(v){W hi=v>>64,lo=v;int k=hi?127-std::countl_zero(hi):63-std::countl_zero(lo);if(piv[k])v^=piv[k];else{piv[k]=v;++r;break;}}
    return r;
}
struct Part{unsigned mask;Perm embed,extract;};
Mat phi(const Perm&u,const Part&p){
    Mat a{};
    for(unsigned v=0;v<8;++v)for(unsigned y=0;y<8;++y){unsigned z=p.extract[u[p.embed[8*v+y]]];unsigned x=z/8,w=z%8;a[8*x+y]|=W(1)<<(8*w+v);}
    return a;
}
int main(int argc,char**argv){
    unsigned degree=argc>1?std::atoi(argv[1]):2;
    unsigned maxvars=argc>2?std::atoi(argv[2]):0;
    std::vector<uint32_t>mon;
    for(unsigned mask=0;mask<32;++mask)if(std::popcount(mask)<=degree){uint32_t v=0;for(unsigned x=0;x<32;++x)if((x&mask)==mask)v|=uint32_t(1)<<x;mon.push_back(v);}
    if(maxvars && maxvars<mon.size())mon.resize(maxvars);
    std::vector<Part>parts;
    for(unsigned mask=1;mask<64;mask+=2)if(std::popcount(mask)==3){
        Part p;p.mask=mask;
        for(unsigned z=0;z<64;++z){unsigned x=0,y=0,kx=0,ky=0;for(unsigned k=0;k<6;++k)if(mask>>k&1)x|=((z>>k)&1)<<kx++;else y|=((z>>k)&1)<<ky++;p.extract[z]=8*x+y;p.embed[8*x+y]=z;}
        parts.push_back(p);
    }
    uint64_t trials=uint64_t(1)<<mon.size(),full=0,jettests=0,survivors=0;
    std::array<uint64_t,7>hist{};uint32_t feedback=1;int best=100000;
    for(uint64_t t=0;t<trials;++t){
        if(t)feedback^=mon[std::countr_zero(t)];
        Perm u,u2;
        for(unsigned x=0;x<64;++x)u[x]=(x>>1)|(((x&1)^((feedback>>(x>>1))&1))<<5);
        unsigned x=0,len=0;do{x=u[x];++len;}while(x && len<64);
        if(len!=64 || x)continue;
        ++full;for(unsigned j=0;j<64;++j)u2[j]=u[u[j]];
        for(const Part&p:parts){
            Mat a=phi(u,p);int r=rank64(a,5);++hist[r];if(r>5)continue;
            Mat c=phi(u2,p);int r2=rank64(c);
            int cheap=64*r+63*std::max(0,r2-2*r);
            if(cheap>=best && cheap>358)continue;
            std::array<Q,128>jet{};for(int i=0;i<64;++i){jet[i]=a[i];jet[64+i]=Q(a[i]^c[i])|(Q(a[i])<<64);}
            int ones=rank128(jet)-2*r,cost=64*r+63*ones;++jettests;
            bool record=cost<best||cost<=358;
            if(cost<best){best=cost;std::cerr<<"best "<<best<<" r "<<r<<" r2 "<<r2<<" ones "<<ones<<" full "<<full<<" trial "<<t<<"\n";}
            if(cost<=358)++survivors;
            if(record)std::cout<<"{\"degree\":"<<degree<<",\"trial\":"<<t<<",\"feedback_mask\":"<<feedback<<",\"partition_mask\":"<<p.mask<<",\"rank\":"<<r<<",\"square_rank\":"<<r2<<",\"ones\":"<<ones<<",\"first_jet_floor\":"<<cost<<"}\n";
        }
    }
    std::cerr<<"DONE trials "<<trials<<" full "<<full<<" partitions "<<parts.size()<<" jets "<<jettests<<" survivors "<<survivors<<" best "<<best<<" histogram ";for(int i=0;i<7;++i)if(hist[i])std::cerr<<i<<":"<<hist[i]<<" ";std::cerr<<"\n";
    return 0;
}
