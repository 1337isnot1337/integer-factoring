// Stochastic exact-rank search in the F4 scalar extension. Not an exhaustive audit.
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <random>
#include <string>
#include <vector>

using Rows = std::array<uint16_t, 162>;
struct Point { uint64_t low=0, high=0; Rows a{}, b{}; };
struct Ranks { int row=0, column=0; int product() const {return row*column;} };

int rank_block(const uint16_t* a, const uint16_t* b, int n) {
    uint16_t pa[9]={}, pb[9]={}; int rank=0;
    for (int i=0; i<n; ++i) {
        uint16_t x=a[i], y=b[i];
        while (x|y) {
            int p=__builtin_ctz(x|y);
            unsigned v=((x>>p)&1)|(((y>>p)&1)<<1);
            if (pa[p]|pb[p]) {
                if (v==1) {x^=pa[p]; y^=pb[p];}
                else if (v==2) {x^=pb[p]; y^=pa[p]^pb[p];}
                else {x^=pa[p]^pb[p]; y^=pa[p];}
            } else {
                if (v==1) {pa[p]=x; pb[p]=y;}
                else if (v==2) {pa[p]=x^y; pb[p]=x;}
                else {pa[p]=y; pb[p]=x^y;}
                ++rank; break;
            }
        }
    }
    return rank;
}

int main(int argc, char** argv) {
    if (argc!=5) {std::cerr<<"usage: search blocks.bin seeds.txt iterations seed\n"; return 1;}
    std::ifstream input(argv[1],std::ios::binary); uint32_t dim, blocks; uint64_t target, dual;
    input.read(reinterpret_cast<char*>(&dim),4); input.read(reinterpret_cast<char*>(&blocks),4);
    input.read(reinterpret_cast<char*>(&target),8); input.read(reinterpret_cast<char*>(&dual),8);
    if (!input || dim!=33 || blocks!=25 || !dual) return 2;
    std::vector<uint32_t> sizes(blocks); input.read(reinterpret_cast<char*>(sizes.data()),4*blocks);
    int length=0; for (auto n:sizes) {if (!n || n>9) return 3; length+=n;}
    if (length!=81) return 3;
    std::vector<uint16_t> raw(2*dim*length); input.read(reinterpret_cast<char*>(raw.data()),2*raw.size());
    if (!input) return 4;
    std::vector<Rows> basis(dim);
    for (unsigned j=0;j<dim;++j) for (unsigned side=0;side<2;++side)
        for (int r=0;r<81;++r) basis[j][81*side+r]=raw[81*(side*dim+j)+r];
    auto toggle=[&](Point& p,int j,unsigned v) {
        if (v&1) {p.low^=1ULL<<j;for (int r=0;r<162;++r)p.a[r]^=basis[j][r];}
        if (v&2) {p.high^=1ULL<<j;for (int r=0;r<162;++r)p.b[r]^=basis[j][r];}
    };
    auto ranks=[&](const Point& p) {
        Ranks r;int offset=0;
        for (auto n:sizes) {r.row+=rank_block(p.a.data()+offset,p.b.data()+offset,n);
            r.column+=rank_block(p.a.data()+81+offset,p.b.data()+81+offset,n);offset+=n;}
        return r;
    };
    std::ifstream seed_input(argv[2]); std::vector<Point> seeds; uint64_t low,high;
    while (seed_input>>low>>high) {
        if ((low|high)>>dim || __builtin_parityll(low&dual)!=1 || __builtin_parityll(high&dual)!=0) return 5;
        Point p; for (unsigned j=0;j<dim;++j) toggle(p,j,((low>>j)&1)|(((high>>j)&1)<<1));
        seeds.push_back(p);
    }
    if (seeds.empty()) return 5;
    const uint64_t iterations=std::stoull(argv[3]), seed=std::stoull(argv[4]);
    std::mt19937_64 rng(seed);std::uniform_real_distribution<double> uniform(0,1);
    int pivot=__builtin_ctzll(dual); std::vector<int> free;
    for (unsigned j=0;j<dim;++j)if (int(j)!=pivot)free.push_back(j);
    Point current=seeds.front(), best=current; Ranks cr=ranks(current), br=cr;
    uint64_t queries=1,accepted=0,restarts=0; auto start=std::chrono::steady_clock::now();
    for (const auto& point:seeds) {
        Ranks r=ranks(point);++queries;
        std::cout<<"{\"status\":\"seed_control\",\"coordinates\":["<<point.low<<","<<point.high
                 <<"],\"ranks\":["<<r.row<<","<<r.column<<"]}\n";
    }
    auto output=[&](const char* status,uint64_t steps) {
        std::cout<<"{\"status\":\""<<status<<"\",\"iterations\":"<<steps<<",\"queries\":"<<queries
        <<",\"accepted\":"<<accepted<<",\"restarts\":"<<restarts<<",\"coordinates\":["<<best.low<<","<<best.high
        <<"],\"ranks\":["<<br.row<<","<<br.column<<"],\"product\":"<<br.product()
        <<",\"seconds\":"<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<"}\n"<<std::flush;
    };
    output("initial",0);
    for (uint64_t step=0;step<iterations;++step) {
        const int epoch=2000;
        if (step%epoch==0) {current=(rng()%4?seeds[rng()%seeds.size()]:best);cr=ranks(current);++queries;++restarts;}
        Point proposed=current; unsigned edits=(rng()%10==0?2+rng()%4:1);
        for (unsigned k=0;k<edits;++k) {int j=free[rng()%free.size()]; unsigned v=1+rng()%3;
            toggle(proposed,j,v);if ((dual>>j)&1)toggle(proposed,pivot,v);}
        Ranks nr=ranks(proposed);++queries;
        double delta=(nr.product()-cr.product())/25.0+0.01*(nr.row+nr.column-cr.row-cr.column);
        double temperature=0.05+2.0*(1.0-double(step%epoch)/epoch);
        if (delta<=0 || uniform(rng)<std::exp(-delta/temperature)) {current=proposed;cr=nr;++accepted;}
        if (cr.product()<br.product() || (cr.product()==br.product() && cr.row+cr.column<br.row+br.column)) {
            best=current;br=cr;output("improved",step+1);
        }
        if ((step+1)%1000000==0)output("progress",step+1);
    }
    output("complete",iterations);
    return 0;
}
