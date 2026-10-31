#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>
#include <limits>

using U=std::uint64_t;
using Clock=std::chrono::steady_clock;

U multiply(U a,U b) {
    U r=0;
    while(b) { unsigned k=__builtin_ctzll(b);r^=a<<k;b&=b-1; }
    return r;
}

unsigned binary_rank(std::vector<U>& a,unsigned rows,unsigned words) {
    std::vector<int> piv(rows,-1);unsigned rank=0;
    for(unsigned i=0;i<rows;++i) {
        U* row=a.data()+std::size_t(i)*words;
        unsigned w=0;
        while(w<words) {
            while(w<words&&!row[w])++w;
            if(w==words)break;
            unsigned col=64*w+__builtin_ctzll(row[w]);
            if(piv[col]<0) {piv[col]=i;++rank;break;}
            const U* other=a.data()+std::size_t(piv[col])*words;
            for(unsigned j=w;j<words;++j)row[j]^=other[j];
        }
    }
    return rank;
}

#include <functional>
#include <map>

struct WordPair {std::vector<unsigned> a,b;};
const int charges[4]={0,1,-1,0};

std::uint64_t factorial(unsigned n){std::uint64_t a=1;for(unsigned i=2;i<=n;++i)a*=i;return a;}

int main(int argc,char**argv) {
    if(argc<3){std::cerr<<"usage: probe_typed_gaussian_charge_blocks seed.txt c0,c1,c2,c3 [maxblock] [all]\n";return 2;}
    std::array<U,256> seed{};std::ifstream in(argv[1]);for(U&v:seed)if(!(in>>v))return 2;
    // These local checks justify every omitted off-charge block and the
    // type-dependent polynomial degree bound used by the packed arithmetic.
    const unsigned occupation[4]={0,1,1,2};
    for(unsigned x=0;x<4;++x)for(unsigned y=0;y<4;++y)
      for(unsigned u=0;u<4;++u)for(unsigned v=0;v<4;++v) {
        U p=seed[((x*4+y)*4+u)*4+v];
        if(p&&charges[x]+charges[u]!=charges[y]+charges[v]) {
            std::cerr<<"seed violates integer charge conservation\n";return 2;
        }
        // In this gauge deg K(x,y,u,v) <= occ(x)+3 occ(v).
        // A common word type therefore has total degree <= 4 total_occ.
        if(p&&63-__builtin_clzll(p)>occupation[x]+3*occupation[v]) {
            std::cerr<<"seed violates the occupied-degree bound\n";return 2;
        }
        if(__builtin_parityll(p)!=int(x==v&&y==u)) {
            std::cerr<<"seed endpoint is not the crossing\n";return 2;
        }
    }
    std::array<unsigned,4> counts{};std::stringstream cs(argv[2]);std::string item;
    for(unsigned i=0;i<4;++i){if(!std::getline(cs,item,','))return 2;counts[i]=std::stoul(item);}
    unsigned m=std::accumulate(counts.begin(),counts.end(),0U),degree=4*(counts[1]+counts[2]+2*counts[3]);
    if(m==0||m>20||degree>62){std::cerr<<"arithmetic guard\n";return 2;}
    unsigned maxblock=argc>3?std::stoul(argv[3]):20000;
    bool allcodes=argc>4&&std::string(argv[4])=="all";
    if(argc>4&&!allcodes){std::cerr<<"the optional final argument must be all\n";return 2;}
    if(allcodes&&degree>16){std::cerr<<"exhaustive-code degree guard\n";return 2;}
    unsigned codes=allcodes?(1U<<(degree/2+1)):0;
    std::array<std::vector<std::uint64_t>,2> code_ranks;for(auto&r:code_ranks)r.assign(codes,0);
    std::uint64_t dimension=factorial(m);for(unsigned a:counts)dimension/=factorial(a);
    if(dimension<2||dimension>std::numeric_limits<std::uint64_t>::max()/dimension) {
        std::cerr<<"dimension or squared-dimension guard\n";return 2;
    }
    std::array<std::vector<std::uint64_t>,2> ranks;for(auto&r:ranks)r.assign(degree+1,0);
    auto start=Clock::now();unsigned block_count=0;unsigned biggest=0;std::uint64_t pair_check[2]={0,0};
    for(unsigned cut=0;cut<2;++cut) {
        int required=cut?2*(int(counts[1])-int(counts[2])):0;
        for(unsigned n0=0;n0<=m;++n0)for(unsigned n1=0;n1+n0<=m;++n1)
          for(unsigned n2=0;n2+n1+n0<=m;++n2)for(unsigned n3=0;n3+n2+n1+n0<=m;++n3) {
            unsigned n4=m-n0-n1-n2-n3;
            if(-2*int(n0)-int(n1)+int(n3)+2*int(n4)!=required)continue;
            std::array<unsigned,5> signature={n0,n1,n2,n3,n4};std::vector<int> target;
            std::uint64_t multiplicity=factorial(m);
            for(unsigned k=0;k<5;++k){multiplicity/=factorial(signature[k]);for(unsigned j=0;j<signature[k];++j)target.push_back(int(k)-2);}
            std::array<unsigned,4> left=counts,right=counts;WordPair current;current.a.resize(m);current.b.resize(m);
            std::vector<WordPair> pairs;bool overflow=false;
            std::function<void(unsigned)> visit=[&](unsigned pos){
                if(overflow)return;
                if(pos==m){pairs.push_back(current);if(pairs.size()>maxblock)overflow=true;return;}
                for(unsigned a=0;a<4;++a)if(left[a])for(unsigned b=0;b<4;++b)if(right[b]) {
                    int q=cut?charges[a]+charges[b]:charges[a]-charges[b];
                    if(q!=target[pos])continue;
                    --left[a];--right[b];current.a[pos]=a;current.b[pos]=b;visit(pos+1);++left[a];++right[b];
                }
            };
            visit(0);
            if(overflow){std::cout<<"{\"blocked\":true,\"cut\":"<<cut<<",\"maxblock\":"<<maxblock<<"}"<<std::endl;return 3;}
            unsigned n=pairs.size();if(!n)continue;
            biggest=std::max(biggest,n);pair_check[cut]+=multiplicity*n;++block_count;
            unsigned words=(n+63)/64;std::size_t cells=std::size_t(n)*words;
            if(double(cells)*(degree+1)*8>1.5e9){std::cerr<<"allocation guard\n";return 3;}
            std::vector<std::vector<U>> matrices(degree+1,std::vector<U>(cells));
            for(unsigned i=0;i<n;++i)for(unsigned j=0;j<n;++j) {
                U p=1;
                for(unsigned pos=0;pos<m&&p;++pos) {
                    unsigned x=pairs[i].a[pos],y,u,v;
                    if(!cut){y=pairs[i].b[pos];u=pairs[j].b[pos];v=pairs[j].a[pos];}
                    else {u=pairs[i].b[pos];y=pairs[j].a[pos];v=pairs[j].b[pos];}
                    p=multiply(p,seed[((x*4+y)*4+u)*4+v]);
                }
                bool endpoint=!cut?i==j:(pairs[i].a==pairs[j].b&&pairs[i].b==pairs[j].a);
                if(__builtin_parityll(p)!=endpoint){std::cerr<<"endpoint failed\n";return 4;}
                while(p){unsigned k=__builtin_ctzll(p);p&=p-1;matrices[k][std::size_t(i)*words+j/64]^=U(1)<<(j%64);}
            }
            std::cout<<"{\"block\":"<<block_count<<",\"cut\":"<<cut<<",\"signature\":[";
            for(unsigned k=0;k<5;++k){if(k)std::cout<<',';std::cout<<signature[k];}
            std::cout<<"],\"size\":"<<n<<",\"multiplicity\":"<<multiplicity<<",\"ranks\":[";
            std::vector<unsigned> local_ranks(degree+1);
            if(allcodes) {
                for(unsigned k=1;k<=degree;k+=2)if(std::any_of(matrices[k].begin(),matrices[k].end(),[](U x){return x;}))return 5;
                std::vector<U> combined(cells),work(cells);unsigned previous=0;
                for(unsigned step=0;step<codes;++step) {
                    unsigned code=step^(step>>1);
                    if(step){unsigned k=2*__builtin_ctz(code^previous);for(std::size_t i=0;i<cells;++i)combined[i]^=matrices[k][i];}
                    work=combined;unsigned r=binary_rank(work,n,words);code_ranks[cut][code]+=multiplicity*r;
                    if(code&&!(code&(code-1)))local_ranks[2*__builtin_ctz(code)]=r;
                    previous=code;
                }
            } else for(unsigned k=0;k<=degree;++k)local_ranks[k]=binary_rank(matrices[k],n,words);
            for(unsigned k=0;k<=degree;++k){unsigned r=local_ranks[k];ranks[cut][k]+=multiplicity*r;if(k)std::cout<<',';std::cout<<r;}
            std::cout<<"],\"seconds\":"<<std::chrono::duration<double>(Clock::now()-start).count()<<"}"<<std::endl;
          }
        if(pair_check[cut]!=dimension*dimension){std::cerr<<"pair count failed "<<pair_check[cut]<<" "<<dimension*dimension<<"\n";return 4;}
    }
    long double cost=0;
    for(unsigned k=0;k<=degree;++k){cost+=std::sqrt((long double)ranks[0][k]*ranks[1][k]);std::cout<<"{\"degree\":"<<k<<",\"R\":"<<ranks[0][k]<<",\"C\":"<<ranks[1][k]<<"}"<<std::endl;}
    if(allcodes) {
        if(code_ranks[0].back()!=dimension*dimension||code_ranks[1].back()!=dimension*dimension){std::cerr<<"endpoint code rank failed\n";return 5;}
        for(unsigned code=0;code<codes;++code)std::cout<<"{\"code\":"<<code<<",\"R\":"<<code_ranks[0][code]<<",\"C\":"<<code_ranks[1][code]<<"}"<<std::endl;
    }
    std::cout<<"{\"done\":true,\"D\":"<<dimension<<",\"cost\":"<<(double)cost<<",\"rho\":"<<(double)(std::log(cost)/(2*std::log((long double)dimension)))
             <<",\"blocks\":"<<block_count<<",\"largest_block\":"<<biggest<<",\"seconds\":"<<std::chrono::duration<double>(Clock::now()-start).count()<<"}"<<std::endl;
    return 0;
}
