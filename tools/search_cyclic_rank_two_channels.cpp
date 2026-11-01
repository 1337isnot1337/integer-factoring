// Exact F2 four-label search of normalized two-channel maps.
// The optional full mode includes every commuting module over F2.
#include <array>
#include <cstdint>
#include <iostream>
#include <vector>
#include <chrono>
#include <set>
#include <string>

using Mat = uint16_t;
using Op = std::array<uint16_t,16>;
constexpr Mat ID = 0x8421;

Mat product(Mat a, Mat b) {
    Mat c=0;
    for(int i=0;i<4;i++) for(int k=0;k<4;k++)
        if((a>>(4*i+k))&1) c ^= ((b>>(4*k))&15)<<(4*i);
    return c;
}
Mat transpose(Mat a) {
    Mat b=0;
    for(int i=0;i<4;i++)for(int j=0;j<4;j++) b|=((a>>(4*i+j))&1)<<(4*j+i);
    return b;
}
bool inverse(Mat a, Mat &inv) {
    uint8_t rows[4];
    for(int i=0;i<4;i++)rows[i]=((a>>(4*i))&15)|(1<<(4+i));
    for(int j=0;j<4;j++) {
        int p=j;while(p<4&&!((rows[p]>>j)&1))p++;
        if(p==4)return false;
        std::swap(rows[p],rows[j]);
        for(int i=0;i<4;i++)if(i!=j&&((rows[i]>>j)&1))rows[i]^=rows[j];
    }
    inv=0;for(int i=0;i<4;i++)inv|=(rows[i]>>4)<<(4*i);
    return true;
}
uint16_t apply(const Op &a, uint16_t v) {
    uint16_t out=0;
    while(v){int k=__builtin_ctz(v);out^=a[k];v&=v-1;}
    return out;
}
Op compose(const Op&a,const Op&b) {
    Op out;for(int i=0;i<16;i++)out[i]=apply(a,b[i]);return out;
}
Op identity(){Op out;for(int i=0;i<16;i++)out[i]=1<<i;return out;}
bool nonzero(const Op&a){for(auto x:a)if(x)return true;return false;}
std::array<uint16_t,16> phi(const Op&a){
    std::array<uint16_t,16> out{};
    for(int x=0;x<4;x++)for(int y=0;y<4;y++)for(int u=0;u<4;u++)for(int v=0;v<4;v++)
        out[4*x+y]|=((a[4*v+y]>>(4*x+u))&1)<<(4*u+v);
    return out;
}
template<class T,int N> int rank(std::array<T,N> rows){
    T pivots[N]{};int count=0;
    for(T x:rows)while(x){int k=63-__builtin_clzll(uint64_t(x));if(pivots[k])x^=pivots[k];else{pivots[k]=x;count++;break;}}
    return count;
}
Op channel(Mat a,Mat b){
    Mat c=ID^product(a,b);Op out{};
    for(int i=0;i<4;i++)for(int j=0;j<4;j++) {
        out[4*i+j]=((c>>(4*i))&15)<<(4*j);
        for(int x=0;x<4;x++)if((a>>(4*x+j))&1)out[4*i+j]^=((b>>(4*i))&15)<<(4*x);
    }
    return out;
}
struct Conjugator {Mat p,inv;Op action;};
int main(int argc,char**argv){
    bool all_commuting=argc>1&&std::string(argv[1])=="--all-commuting";
    auto started=std::chrono::steady_clock::now();
    std::vector<Conjugator> conjugators;
    for(unsigned p=0;p<65536;p++){Mat inv;if(!inverse(p,inv))continue;Conjugator g{Mat(p),inv,{}};
        for(int i=0;i<16;i++)g.action[i]=product(product(p,Mat(1<<i)),inv);
        conjugators.push_back(g);
    }
    std::vector<std::pair<Mat,Mat>> factors;
    int conjugacy_classes=0;
    uint64_t commuting_pairs=0;
    if(all_commuting){
        std::array<bool,65536> seen{};
        for(unsigned a=0;a<65536;a++)if(!seen[a]){
            conjugacy_classes++;
            uint64_t orbit_size=0,centralizer_size=0;
            for(const auto&g:conjugators){Mat x=apply(g.action,Mat(a));if(!seen[x]){seen[x]=true;orbit_size++;}}
            if(a==0||a==ID){commuting_pairs+=65536;
                std::cout<<"{\"class\":true,\"A\":"<<a<<",\"orbit_size\":"<<orbit_size<<",\"commuting_matrices\":65536}"<<std::endl;
                continue;}
            for(unsigned b=0;b<65536;b++)if(product(a,b)==product(b,a)){
                commuting_pairs++;centralizer_size++;factors.emplace_back(Mat(a),Mat(b));
            }
            std::cout<<"{\"class\":true,\"A\":"<<a<<",\"orbit_size\":"<<orbit_size
                     <<",\"commuting_matrices\":"<<centralizer_size<<"}"<<std::endl;
        }
    }else{
        for(int polynomial=0;polynomial<16;polynomial++){
            Mat a=0;for(int i=0;i<3;i++)a|=1<<(4*(i+1)+i);
            for(int i=0;i<4;i++)if((polynomial>>i)&1)a|=1<<(4*i+3);
            Mat powers[4]={ID,a,product(a,a),0};powers[3]=product(powers[2],a);
            for(int bmask=0;bmask<16;bmask++){
                Mat b=0;for(int i=0;i<4;i++)if((bmask>>i)&1)b^=powers[i];
                factors.emplace_back(a,b);
            }
        }
    }
    uint64_t tested=0,regular=0,trace_rejected=0,histogram[17]{},duplicates=0,rank_one_bases=0;int best=17;
    std::set<Op> distinct_bases;
    Op one=identity();
    for(size_t factor_index=0;factor_index<factors.size();factor_index++){
            auto [a,b]=factors[factor_index];
            Op base=channel(a,b);
            if(all_commuting){
                if(rank<uint16_t,16>(phi(base))<2){rank_one_bases++;continue;}
                if(!distinct_bases.insert(base).second){duplicates++;continue;}
            }
            for(const auto &g:conjugators){
                tested++;Op u,s;unsigned trace=0;
                for(int i=0;i<16;i++){u[i]=apply(g.action,base[i]);s[i]=u[i]^(1<<i);trace^=(s[i]>>i)&1;}
                if(trace){trace_rejected++;continue;}
                Op s2=compose(s,s),s4=compose(s2,s2),s8=compose(s4,s4),s16=compose(s8,s8);
                if(nonzero(s16))continue;
                Op s15=compose(s8,compose(s4,compose(s2,s)));
                if(!nonzero(s15))continue;
                regular++;
                Op u2;for(int i=0;i<16;i++)u2[i]=s2[i]^(1<<i);
                auto h0=phi(u),h1=phi(u2);
                for(int i=0;i<16;i++)h1[i]^=h0[i];
                int r0=rank<uint16_t,16>(h0);
                std::array<uint32_t,32> toeplitz{};
                for(int i=0;i<16;i++){toeplitz[i]=h0[i];toeplitz[i+16]=uint32_t(h1[i])|(uint32_t(h0[i])<<16);}
                int ones=rank<uint32_t,32>(toeplitz)-2*r0;
                histogram[ones]++;
                if(all_commuting)std::cout<<"{\"regular\":true,\"A\":"<<a<<",\"B\":"<<b<<",\"P\":"<<g.p
                    <<",\"zero_orders\":"<<r0<<",\"order_one_count\":"<<ones<<"}\n";
                if(ones<best || ones<=1){
                    best=std::min(best,ones);
                    std::cout<<"{\"candidate\":true,\"factor_index\":"<<factor_index
                             <<",\"A\":"<<a<<",\"B\":"<<b<<",\"P\":"<<g.p<<",\"P_inverse\":"<<g.inv
                             <<",\"zero_orders\":"<<r0<<",\"order_one_count\":"<<ones<<",\"U_columns\":[";
                    for(int i=0;i<16;i++)std::cout<<(i?",":"")<<u[i];
                    std::cout<<"]}"<<std::endl;
                }
            }
        if(factor_index%128==127||factor_index+1==factors.size())
            std::cerr<<"completed factor "<<factor_index+1<<"/"<<factors.size()<<", tested "<<tested<<", regular "<<regular<<", best ones "<<best<<std::endl;
    }
    std::cout<<"{\"summary\":true,\"all_commuting\":"<<(all_commuting?"true":"false")
             <<",\"conjugacy_classes\":"<<conjugacy_classes<<",\"canonical_commuting_pairs\":"<<commuting_pairs
             <<",\"factor_pairs\":"<<factors.size()<<",\"distinct_rank_two_bases\":"<<distinct_bases.size()
             <<",\"duplicate_bases\":"<<duplicates<<",\"rank_one_bases_skipped\":"<<rank_one_bases
             <<",\"conjugators\":"<<conjugators.size()<<",\"tested\":"<<tested
             <<",\"regular\":"<<regular<<",\"trace_rejected\":"<<trace_rejected<<",\"order_one_histogram\":[";
    for(int i=0;i<17;i++)std::cout<<(i?",":"")<<histogram[i];
    std::cout<<"],\"seconds\":"<<std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count()<<"}"<<std::endl;
}
