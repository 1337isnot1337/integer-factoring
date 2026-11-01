// Complete Gray-code coverage of a declared affine F4 subspace.
// Reuse the independently controlled exact packed F4 rank implementation.
#define main fused_f4_stochastic_main
#include "search_fused_f4_envelope.cpp"
#undef main

int main(int argc,char** argv) {
    if (argc!=9) {std::cerr<<"usage: scan blocks.bin subspace.txt retained.jsonl start stop mask threshold affine_dual_value\n";return 1;}
    std::ifstream input(argv[1],std::ios::binary); uint32_t dim,blocks;uint64_t target,dual;
    input.read(reinterpret_cast<char*>(&dim),4);input.read(reinterpret_cast<char*>(&blocks),4);
    input.read(reinterpret_cast<char*>(&target),8);input.read(reinterpret_cast<char*>(&dual),8);
    if (!input || dim!=33 || blocks!=25 || !dual)return 2;
    std::vector<uint32_t> sizes(blocks);input.read(reinterpret_cast<char*>(sizes.data()),4*blocks);
    int length=0;for(auto n:sizes){if(!n||n>9)return 3;length+=n;}if(length!=81)return 3;
    std::vector<uint16_t> raw(2*dim*81);input.read(reinterpret_cast<char*>(raw.data()),2*raw.size());if(!input)return 4;
    std::vector<Rows> basis(dim);
    for(unsigned j=0;j<dim;++j)for(unsigned s=0;s<2;++s)for(int r=0;r<81;++r)basis[j][81*s+r]=raw[81*(s*dim+j)+r];
    std::ifstream declaration(argv[2]);uint64_t offset,coordinate;unsigned nfree;
    declaration>>offset>>nfree;if(!declaration || !nfree || nfree>30 || offset>>dim)return 5;
    std::vector<uint64_t> directions(nfree);for(auto& x:directions){declaration>>x;if(!declaration||x>>dim)return 5;}
    auto compile=[&](uint64_t code){Rows out{};for(unsigned j=0;j<dim;++j)if((code>>j)&1)for(int r=0;r<162;++r)out[r]^=basis[j][r];return out;};
    std::vector<Rows> generators;for(auto x:directions)generators.push_back(compile(x));
    const uint64_t start=std::stoull(argv[4]),stop=std::stoull(argv[5]),mask=std::stoull(argv[6]);
    const int threshold=std::stoi(argv[7]); // argv[8] fixes an independently checked affine dual value.
    const unsigned affine_value=std::stoul(argv[8]);
    if(start>stop || stop>(1ULL<<(2*nfree)) || mask>>(2*nfree) || threshold<0 || affine_value>1)return 6;
    if(__builtin_parityll(offset&dual)!=int(affine_value))return 7;
    for(auto x:directions)if(__builtin_parityll(x&dual))return 7;
    Point point;point.low=offset;point.a=compile(offset);
    auto toggle=[&](unsigned bit){unsigned j=bit%nfree;bool high=bit>=nfree;
        auto& rows=high?point.b:point.a;auto& code=high?point.high:point.low;code^=directions[j];
        for(int r=0;r<162;++r)rows[r]^=generators[j][r];};
    uint64_t initial=(start^(start>>1))^mask;for(unsigned bit=0;bit<2*nfree;++bit)if((initial>>bit)&1)toggle(bit);
    std::ofstream output(argv[3]);if(!output)return 8;
    uint64_t retained=0,rejected=0,zeros=0,scalar_skips=0;int minimum=100000;auto clock=std::chrono::steady_clock::now();
    for(uint64_t position=start;position<stop;++position){
        if(!point.low && !point.high){++zeros;}
        else if (!affine_value && [&](){unsigned p=63-__builtin_clzll(point.low|point.high);
            return (((point.low>>p)&1)|(((point.high>>p)&1)<<1))!=1;}()) {++scalar_skips;}
        else {
            int rr=0,rc=0,offset=0;bool reject=false;
            for(auto n:sizes){rr+=rank_block(point.a.data()+offset,point.b.data()+offset,n);
                rc+=rank_block(point.a.data()+81+offset,point.b.data()+81+offset,n);offset+=n;
                if(rr*rc>threshold){reject=true;break;}}
            if(reject)++rejected;
            else {++retained;minimum=std::min(minimum,rr*rc);
                output<<"{\"position\":"<<position<<",\"coordinates\":["<<point.low<<","<<point.high
                      <<"],\"ranks\":["<<rr<<","<<rc<<"]}\n";}
        }
        if(position+1<stop)toggle(__builtin_ctzll(position+1));
    }
    output.close();if(!output)return 8;
    std::cout<<"{\"status\":\"complete\",\"start\":"<<start<<",\"stop\":"<<stop<<",\"visited\":"<<stop-start
        <<",\"retained\":"<<retained<<",\"rejected\":"<<rejected<<",\"zeros\":"<<zeros
        <<",\"scalar_orbit_skips\":"<<scalar_skips
        <<",\"minimum_retained_product\":"<<minimum<<",\"seconds\":"
        <<std::chrono::duration<double>(std::chrono::steady_clock::now()-clock).count()<<"}\n";
}
