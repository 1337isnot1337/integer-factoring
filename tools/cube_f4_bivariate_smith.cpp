// Exact x-adic Smith row certificates over F4(y), with degree-six inputs.
// F4 = F2[alpha]/(alpha^2 + alpha + 1).  No specialization is used.
#include <algorithm>
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

using Word = std::uint64_t;
using Binary = std::vector<Word>;
constexpr int STRIDE = 256;

void trim(Binary &a) { while (!a.empty() && !a.back()) a.pop_back(); }
int degree(const Binary &a) {
    return a.empty() ? -1 : int(64 * (a.size()-1) + 63 - __builtin_clzll(a.back()));
}
Binary parse(const std::string &s) {
    if (s.empty() || s.find_first_not_of("0123456789abcdefABCDEF") != std::string::npos)
        throw std::runtime_error("invalid hexadecimal polynomial");
    Binary a;
    for (int end=int(s.size()); end>0; end-=16) {
        int begin=std::max(0,end-16);
        a.push_back(std::stoull(s.substr(begin,end-begin),nullptr,16));
    }
    trim(a); return a;
}
std::string hexadecimal(const Binary &a) {
    if (a.empty()) return "0";
    std::ostringstream out; out << std::hex << a.back();
    for (int i=int(a.size())-2;i>=0;--i)
        out << std::setfill('0') << std::setw(16) << a[i];
    return out.str();
}
void xor_shift(Binary &a, const Binary &b, int shift) {
    if (b.empty()) return;
    int word=shift/64, bit=shift%64;
    a.resize(std::max(a.size(),b.size()+word+(bit!=0)),0);
    for (std::size_t i=0;i<b.size();++i) {
        a[i+word]^=b[i]<<bit;
        if (bit) a[i+word+1]^=b[i]>>(64-bit);
    }
    trim(a);
}
Binary add(Binary a,const Binary &b) { xor_shift(a,b,0);return a; }
Binary multiply(const Binary &a,const Binary &b) {
    if (a.empty() || b.empty()) return {};
    Binary c(a.size()+b.size(),0);
    for (std::size_t i=0;i<a.size();++i) if (a[i])
        for (std::size_t j=0;j<b.size();++j) if (b[j]) {
#if defined(__PCLMUL__)
            __m128i u=_mm_set_epi64x(0,(long long)a[i]);
            __m128i v=_mm_set_epi64x(0,(long long)b[j]);
            __m128i w=_mm_clmulepi64_si128(u,v,0);
            alignas(16) Word z[2]; _mm_store_si128((__m128i*)z,w);
            c[i+j]^=z[0];c[i+j+1]^=z[1];
#else
            Word value=a[i];
            while (value) {
                int bit=__builtin_ctzll(value);value&=value-1;
                c[i+j]^=b[j]<<bit;
                if (bit) c[i+j+1]^=b[j]>>(64-bit);
            }
#endif
        }
    trim(c);return c;
}
int x_order(const Binary &a) {
    int order=STRIDE;
    for (std::size_t i=0;i<a.size();++i) if (a[i]) {
        order=std::min(order,int(i%4)*64+__builtin_ctzll(a[i]));
        if (!order) break;
    }
    return order;
}
Binary divide_x(const Binary &a,int power) {
    if (!power) return a;
    if (x_order(a)<power) throw std::runtime_error("nonregular normalized row frame");
    int word=power/64,bit=power%64;
    Binary result(a.size()>std::size_t(word)?a.size()-word:0,0);
    for (std::size_t j=0;j<result.size();++j) {
        result[j]=a[j+word]>>bit;
        if (bit && j+word+1<a.size()) result[j]^=a[j+word+1]<<(64-bit);
    }
    trim(result);return result;
}
void validate_input(const Binary &a) {
    for (std::size_t i=0;i<a.size();++i) {
        Word value=a[i];
        while (value) {
            int at=int(i*64)+__builtin_ctzll(value);value&=value-1;
            if (at%STRIDE+at/STRIDE>6)
                throw std::runtime_error("input total degree exceeds six");
        }
    }
}

struct Poly {
    Binary real,alpha;
    bool zero() const { return real.empty() && alpha.empty(); }
};
Poly add(Poly a,const Poly &b) {
    xor_shift(a.real,b.real,0);xor_shift(a.alpha,b.alpha,0);return a;
}
Poly multiply(const Poly &a,const Poly &b) {
    if (a.zero() || b.zero()) return {};
    if (a.alpha.empty()) return {multiply(a.real,b.real),multiply(a.real,b.alpha)};
    if (b.alpha.empty()) return {multiply(a.real,b.real),multiply(a.alpha,b.real)};
    Binary p=multiply(a.real,b.real),q=multiply(a.alpha,b.alpha);
    Binary cross=multiply(add(a.real,a.alpha),add(b.real,b.alpha));
    return {add(p,q),add(cross,p)};
}
int degree(const Poly &a) { return std::max(degree(a.real),degree(a.alpha)); }
int leading(const Poly &a,int at) {
    auto coefficient=[&](const Binary &b) {
        return at>=0 && std::size_t(at/64)<b.size() ? int((b[at/64]>>(at%64))&1) : 0;
    };
    return coefficient(a.real)+2*coefficient(a.alpha);
}
int scalar_product(int a,int b) {
    int p=(a&1)*(b&1),q=(a>>1)*(b>>1);
    return (p^q)+2*(((a&1)^(a>>1))*((b&1)^(b>>1))^p);
}
Poly scale(const Poly &a,int coefficient) {
    if (coefficient==0) return {};
    if (coefficient==1) return a;
    if (coefficient==2) return {a.alpha,add(a.real,a.alpha)};
    return {add(a.real,a.alpha),a.real};
}
Poly divide_exact(Poly a,const Poly &b) {
    if (b.zero()) throw std::runtime_error("zero polynomial divisor");
    if (b.alpha.empty() && b.real.size()==1 && b.real[0]==1) return a;
    int db=degree(b),lb=leading(b,db),inverse=scalar_product(lb,lb);
    Poly multiples[4];for(int c=1;c<4;++c)multiples[c]=scale(b,c);
    Poly quotient;
    while (degree(a)>=db) {
        int da=degree(a),shift=da-db,c=scalar_product(leading(a,da),inverse);
        if (c&1) { if(quotient.real.size()<=std::size_t(shift/64))quotient.real.resize(shift/64+1,0);
                   quotient.real[shift/64]^=Word(1)<<(shift%64); }
        if (c&2) { if(quotient.alpha.size()<=std::size_t(shift/64))quotient.alpha.resize(shift/64+1,0);
                   quotient.alpha[shift/64]^=Word(1)<<(shift%64); }
        xor_shift(a.real,multiples[c].real,shift);
        xor_shift(a.alpha,multiples[c].alpha,shift);
    }
    if (!a.zero()) throw std::runtime_error("nonexact F4 Bareiss division");
    trim(quotient.real);trim(quotient.alpha);return quotient;
}
int x_order(const Poly &a) { return std::min(x_order(a.real),x_order(a.alpha)); }
Poly divide_x(const Poly &a,int power) { return {divide_x(a.real,power),divide_x(a.alpha,power)}; }

int main() {
    try {
        int n;if (!(std::cin>>n) || n<1 || n>16) throw std::runtime_error("matrix order outside 1..16");
        using Matrix=std::vector<std::vector<Poly>>;
        Matrix a(n,std::vector<Poly>(n)),frame(n,std::vector<Poly>(n));
        for (int i=0;i<n;++i) {
            frame[i][i].real={1};
            for (int j=0;j<n;++j) {
                std::string u,v;if (!(std::cin>>u>>v))throw std::runtime_error("missing matrix entry");
                a[i][j]={parse(u),parse(v)};
                validate_input(a[i][j].real);validate_input(a[i][j].alpha);
            }
        }
        std::vector<bool> used_row(n),used_column(n);
        Poly previous{{1},{}};
        std::vector<int> orders,rows,columns,prefixes;
        Matrix certificate;
        for (int k=0;k<n;++k) {
            int bi=-1,bj=-1,bv=100000,pv=x_order(previous);
            for (int i=0;i<n;++i) if (!used_row[i])
                for (int j=0;j<n;++j) if (!used_column[j] && !a[i][j].zero()) {
                    int value=x_order(a[i][j]);
                    if (value<bv) { bv=value;bi=i;bj=j; }
                }
            if (bi<0) throw std::runtime_error("generic rank deficit");
            Poly pivot=a[bi][bj];
            orders.push_back(bv-pv);prefixes.push_back(bv);rows.push_back(bi);columns.push_back(bj);
            std::vector<Poly> row;
            for (const auto &value:frame[bi]) row.push_back(divide_x(value,pv));
            certificate.push_back(std::move(row));used_row[bi]=true;used_column[bj]=true;
            if (k+1<n) {
                Matrix next=a,nf=frame;
                for (int i=0;i<n;++i) if (!used_row[i]) {
                    for (int j=0;j<n;++j) if (!used_column[j])
                        next[i][j]=divide_exact(add(multiply(pivot,a[i][j]),multiply(a[i][bj],a[bi][j])),previous);
                    for (int j=0;j<n;++j)
                        nf[i][j]=divide_exact(add(multiply(pivot,frame[i][j]),multiply(a[i][bj],frame[bi][j])),previous);
                }
                a.swap(next);frame.swap(nf);
            }
            previous=std::move(pivot);
        }
        std::cout<<"{\"field\":\"F4\",\"stride\":256,\"orders\":[";
        for(int i=0;i<n;++i)std::cout<<(i?",":"")<<orders[i];
        std::cout<<"],\"rows\":[";for(int i=0;i<n;++i)std::cout<<(i?",":"")<<rows[i];
        std::cout<<"],\"columns\":[";for(int i=0;i<n;++i)std::cout<<(i?",":"")<<columns[i];
        std::cout<<"],\"prefixes\":[";for(int i=0;i<n;++i)std::cout<<(i?",":"")<<prefixes[i];
        std::cout<<"],\"frames\":[";
        for(int i=0;i<n;++i) {
            if(i)std::cout<<',';std::cout<<'[';
            for(int j=0;j<n;++j) {
                if(j)std::cout<<',';
                std::cout<<"[\""<<hexadecimal(certificate[i][j].real)<<"\",\""
                         <<hexadecimal(certificate[i][j].alpha)<<"\"]";
            }
            std::cout<<']';
        }
        std::cout<<"]}\n";return 0;
    } catch (const std::exception &error) {
        std::cerr<<error.what()<<'\n';return 2;
    }
}
