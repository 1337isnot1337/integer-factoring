// Complete F3 binary cubic crossings, pure at 0,2,infinity and Swap at 1.
// Local Smith orders use minimum-valuation pivots, with precision 4*3+1.
#include <algorithm>
#include <array>
#include <cassert>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <vector>

using Tensor = std::array<int,16>;
using Series = std::array<int,13>;
using Matrix = std::array<std::array<Series,4>,4>;
using Signature = std::array<int,24>;
int mod(int x) { x %= 3; return x < 0 ? x+3 : x; }
int index(int x,int y,int u,int v) { return 8*x+4*y+2*u+v; }

Series multiply(const Series &a,const Series &b) {
    Series c{};
    for(int i=0;i<13;++i) if(a[i])
        for(int j=0;i+j<13;++j) if(b[j]) c[i+j]=(c[i+j]+a[i]*b[j])%3;
    return c;
}

std::array<int,4> smith(Matrix a) {
    std::array<int,4> orders{};
    for(int k=0;k<4;++k) {
        int order=13,pi=-1,pj=-1;
        for(int i=k;i<4;++i) for(int j=k;j<4;++j) {
            int v=0;
            while(v<13 && !a[i][j][v]) ++v;
            if(v<order) { order=v;pi=i;pj=j; }
        }
        assert(order<13);
        orders[k]=order;
        std::swap(a[k],a[pi]);
        for(auto &row:a) std::swap(row[k],row[pj]);
        Series unit{},inverse{};
        for(int i=0;i+order<13;++i) unit[i]=a[k][k][i+order];
        inverse[0]=unit[0]; // Every nonzero element of F3 is its own inverse.
        for(int i=1;i<13;++i) {
            int sum=0;
            for(int j=1;j<=i;++j) sum+=unit[j]*inverse[i-j];
            inverse[i]=mod(-inverse[0]*sum);
        }
        for(int j=k+1;j<4;++j) a[k][j]=multiply(a[k][j],inverse);
        for(int i=k+1;i<4;++i) {
            Series factor{};
            for(int j=0;j+order<13;++j) factor[j]=a[i][k][j+order];
            for(int j=k+1;j<4;++j) {
                auto correction=multiply(factor,a[k][j]);
                for(int v=0;v<13;++v) a[i][j][v]=mod(a[i][j][v]-correction[v]);
            }
        }
    }
    assert(std::is_sorted(orders.begin(),orders.end()));
    assert(std::accumulate(orders.begin(),orders.end(),0)<=12);
    return orders;
}

Signature profiles(const Tensor &a,const Tensor &b,const Tensor &c) {
    std::array<Tensor,4> coefficient{};
    for(int x=0;x<2;++x) for(int y=0;y<2;++y)
    for(int u=0;u<2;++u) for(int v=0;v<2;++v) {
        int k=index(x,y,u,v),e=(x==v && y==u);
        coefficient[0][k]=a[k];
        coefficient[1][k]=mod(b[k]-e-c[k]);
        coefficient[2][k]=mod(-a[k]-b[k]-e);
        coefficient[3][k]=c[k];
        assert(mod(coefficient[0][k]+coefficient[1][k]+coefficient[2][k]+coefficient[3][k])==e);
    }
    const int binomial[4][4]={{1,0,0,0},{1,1,0,0},{1,2,1,0},{1,0,0,1}};
    Signature signature{};
    int output=0;
    for(int place=0;place<3;++place) {
        std::array<Tensor,4> local{};
        for(int k=0;k<16;++k) for(int j=0;j<4;++j) {
            if(place==2) local[j][k]=coefficient[3-j][k];
            else for(int l=j;l<4;++l) {
                int power=1;
                for(int t=0;t<l-j;++t) power*=place==0?0:2;
                local[j][k]=mod(local[j][k]+coefficient[l][k]*binomial[l][j]*power);
            }
        }
        for(int cut=0;cut<2;++cut) {
            Matrix matrix{};
            for(int i=0;i<4;++i) for(int j=0;j<4;++j) {
                int k=cut==0?index(i/2,i%2,j/2,j%2):index(i/2,j/2,i%2,j%2);
                for(int t=0;t<4;++t) matrix[i][j][t]=local[t][k];
            }
            for(int order:smith(matrix)) signature[output++]=order;
        }
    }
    return signature;
}

struct Certificate { int n,d,k,l; };
Certificate certificate(const Signature &s,int place) {
    Certificate best{100,1,0,0};
    for(int k=1;k<=4;++k) for(int l=1;l<=4;++l) if(k*l>=8) {
        int a=0,b=0;
        for(int i=0;i<k;++i) a+=s[8*place+i];
        for(int i=0;i<l;++i) b+=s[8*place+4+i];
        int n=a,d=k;
        if(b*k>a*l) { n=b;d=l; }
        int g=std::gcd(n,d);n/=g;d/=g;
        if(n*best.d<best.n*d) best={n,d,k,l};
    }
    return best;
}

struct Record { long long count=0;std::array<int,3> witness{}; };
int main(int argc,char **argv) {
    assert(argc==2 || argc==3);
    long long limit=argc==3?std::stoll(argv[2]):0;
    std::vector<Tensor> values(1),representatives(1);
    const int directions[4][2]={{1,0},{0,1},{1,1},{1,2}};
    for(int i=0;i<4;++i) for(int j=0;j<4;++j)
    for(int k=0;k<4;++k) for(int l=0;l<4;++l) for(int scale=1;scale<=2;++scale) {
        Tensor value{};
        for(int x=0;x<2;++x) for(int y=0;y<2;++y)
        for(int u=0;u<2;++u) for(int v=0;v<2;++v)
            value[index(x,y,u,v)]=scale*directions[i][x]*directions[j][y]*directions[k][u]*directions[l][v]%3;
        values.push_back(value);
    }
    std::sort(values.begin(),values.end());
    assert(values.size()==513 && std::adjacent_find(values.begin(),values.end())==values.end());
    for(int u=0;u<2;++u) for(int v=0;v<2;++v)
        for(int scale=1;scale<=(u==0 && v==0?2:1);++scale) {
            Tensor value{};value[index(0,0,u,v)]=scale;representatives.push_back(value);
        }
    std::map<Signature,Record> catalogue;
    long long processed=0,total=6LL*513*513;
    for(int a=0;a<6;++a) for(int b=0;b<513;++b) for(int c=0;c<513;++c) {
        if(limit && processed>=limit) goto done;
        auto signature=profiles(representatives[a],values[b],values[c]);
        auto &record=catalogue[signature];
        if(!record.count) record.witness={a,b,c};
        ++record.count;++processed;
        if(processed%100000==0) std::cerr<<processed<<" / "<<total<<" seeds\n";
    }
done:
    std::ofstream out(argv[1]);
    assert(out);
    out<<"{\"degree\":3,\"processed\":"<<processed<<",\"total\":"<<total
       <<",\"complete\":"<<(processed==total?"true":"false")<<",\"catalogue\":[\n";
    bool first=true;int failures=0;
    for(const auto &[s,record]:catalogue) {
        if(!first) out<<",\n";
        first=false;
        out<<"{\"profiles\":[";
        for(int cut=0;cut<6;++cut) {
            if(cut) out<<",";
            out<<"[";
            for(int i=0;i<4;++i) { if(i) out<<",";out<<s[4*cut+i]; } out<<"]";
        }
        out<<"],\"count\":"<<record.count<<",\"indices\":["<<record.witness[0]<<","<<record.witness[1]<<","<<record.witness[2]
           <<"],\"prefix_certificates\":[";
        int numerator=0,denominator=1;
        for(int place=0;place<3;++place) {
            auto c=certificate(s,place);
            numerator=numerator*c.d+c.n*denominator;denominator*=c.d;
            if(place) out<<",";
            out<<"["<<c.n<<","<<c.d<<","<<c.k<<","<<c.l<<"]";
        }
        if(numerator>3*denominator) ++failures;
        out<<"]}";
    }
    out<<"],\"signatures\":"<<catalogue.size()<<",\"prefix_certificate_failures\":"<<failures<<"}\n";
    std::cerr<<"processed "<<processed<<", signatures "<<catalogue.size()<<", prefix failures "<<failures<<"\n";
}
