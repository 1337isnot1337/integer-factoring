#define main debruijn_unused_main
#include "search_debruijn_q8.cpp"
#undef main
#include <set>

int smallrank(std::vector<unsigned>v){
    std::array<unsigned,6>p{};int r=0;
    for(unsigned x:v)while(x){int k=31-std::countl_zero(x);if(p[k])x^=p[k];else{p[k]=x;++r;break;}}
    return r;
}
struct Space{std::array<unsigned,3>a,b;};
int main(){
    std::vector<Space>spaces;std::set<W>seen;
    for(unsigned a=1;a<64;++a)for(unsigned b=a+1;b<64;++b)for(unsigned c=b+1;c<64;++c){
        if(c==(a^b))continue;
        W key=1;
        for(unsigned t=1;t<8;++t){unsigned x=((t&1)?a:0)^((t&2)?b:0)^((t&4)?c:0);key|=W(1)<<x;}
        if(!seen.insert(key).second)continue;
        Space s;s.a={a,b,c};std::vector<unsigned>basis={a,b,c};int j=0;
        for(unsigned x=1;x<64;x<<=1){auto v=basis;v.push_back(x);if(smallrank(v)>int(basis.size())){s.b[j++]=x;basis.push_back(x);}}
        if(j!=3){std::cerr<<"bad complement\n";return 1;}spaces.push_back(s);
    }
    if(spaces.size()!=1395){std::cerr<<"bad space count "<<spaces.size()<<"\n";return 1;}
    std::array<unsigned,6>taps={3,27,33,39,45,51};
    std::array<Perm,6>linear,u,u2;
    for(unsigned t=0;t<6;++t){
        for(unsigned x=0;x<64;++x)linear[t][x]=(x>>1)|((std::popcount(x&taps[t])&1)<<5);
        u[t]=linear[t];u[t][0]=32;u[t][1]=0;
        for(unsigned x=0;x<64;++x)u2[t][x]=u[t][u[t][x]];
        unsigned x=0,l=0;do{x=u[t][x];++l;}while(x&&l<64);if(l!=64||x)return 2;
    }
    uint64_t gauges=0,cases=0,linearpass=0,leadingpass=0,jets=0,survivors=0;
    std::array<uint64_t,6>rankhist{};std::array<int,6>best_by_rank;best_by_rank.fill(100000);
    int best=100000;
    for(const Space&s:spaces){
        std::array<unsigned,8>aemb{};for(unsigned t=1;t<8;++t){unsigned bit=std::countr_zero(t);aemb[t]=aemb[t^(1U<<bit)]^s.a[bit];}
        for(unsigned graph=0;graph<512;++graph){
            std::array<unsigned,3>b;for(unsigned j=0;j<3;++j)b[j]=s.b[j]^aemb[(graph>>(3*j))&7];
            std::array<unsigned,8>bemb{};for(unsigned t=1;t<8;++t){unsigned bit=std::countr_zero(t);bemb[t]=bemb[t^(1U<<bit)]^b[bit];}
            Part p;p.mask=graph;for(unsigned i=0;i<8;++i)for(unsigned j=0;j<8;++j){unsigned x=aemb[i]^bemb[j];p.embed[8*i+j]=x;p.extract[x]=8*i+j;}
            ++gauges;
            for(unsigned t=0;t<6;++t){
                ++cases;std::vector<unsigned>aa,dd;
                for(unsigned j=0;j<3;++j){aa.push_back(p.extract[linear[t][s.a[j]]]>>3);dd.push_back(p.extract[linear[t][b[j]]]&7);}
                unsigned e=p.extract[1],f=p.extract[32];
                int correction_rank=((e/8&&e%8)?2:1)*((f/8&&f%8)?2:1);
                if((1<<(smallrank(aa)+smallrank(dd)))-correction_rank>5)continue;
                ++linearpass;Mat A=phi(u[t],p);int r=rank64(A,5);if(r>5)continue;
                ++leadingpass;++rankhist[r];Mat C=phi(u2[t],p);int r2=rank64(C);
                int cheap=64*r+63*std::max(0,r2-2*r);
                if(cheap>=best_by_rank[r]&&cheap>358)continue;
                std::array<Q,128>jet{};for(unsigned i=0;i<64;++i){jet[i]=A[i];jet[64+i]=Q(A[i]^C[i])|(Q(A[i])<<64);}
                int ones=rank128(jet)-2*r,cost=64*r+63*ones;++jets;
                if(cost<best_by_rank[r])best_by_rank[r]=cost;
                if(cost<best||cost<=358){
                    if(cost<best){best=cost;std::cerr<<"best "<<best<<" r "<<r<<" square "<<r2<<" ones "<<ones<<" cases "<<cases<<"\n";}
                    if(cost<=358)++survivors;
                    std::cout<<"{\"taps\":"<<taps[t]<<",\"row_basis\":["<<s.a[0]<<","<<s.a[1]<<","<<s.a[2]<<"],\"column_basis\":["<<b[0]<<","<<b[1]<<","<<b[2]<<"],\"rank\":"<<r<<",\"square_rank\":"<<r2<<",\"ones\":"<<ones<<",\"first_jet_floor\":"<<cost<<"}\n";
                }
            }
        }
    }
    std::cerr<<"DONE spaces "<<spaces.size()<<" gauges "<<gauges<<" cases "<<cases<<" linearpass "<<linearpass<<" leadingpass "<<leadingpass<<" jets "<<jets<<" survivors "<<survivors<<" best "<<best<<" histogram ";for(int i=0;i<6;++i)if(rankhist[i])std::cerr<<i<<":"<<rankhist[i]<<":"<<best_by_rank[i]<<" ";std::cerr<<"\n";
}
