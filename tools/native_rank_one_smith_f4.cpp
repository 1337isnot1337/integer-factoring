// Structured finite-series search engine for D+u*v^T and arrowhead matrices.
// Every pivot is a global minimum of the represented matrix's valuations.
// This is a search accelerator: retained profiles require full-matrix replay.
#include "native_rectangular_smith_f4.cpp"

namespace rank_one_detail {
using S=f4::Series;
S norm(S a){
 if(!a.bits)return {};
 int v=f4::valuation(a.bits);a.base+=v;a.bits=a.bits>>v;return a;
}
S add(S a,S b){
 if(!a.bits)return b;if(!b.bits)return a;
 if(a.base>b.base)std::swap(a,b);
 a.bits^=b.bits<<(b.base-a.base);return norm(a);
}
S mul(S a,S b){
 if(!a.bits||!b.bits)return {};
 return norm({f4::multiply(a.bits,b.bits),a.base+b.base});
}
S div(S a,S b){
 if(!b.bits)throw std::runtime_error("zero structured divisor");
 if(!a.bits)return {};
 unsigned invs[4]={0,1,3,2};f4::Poly inverse(invs[b.bits.constant()]);
 if(b.bits!=f4::Poly(1))for(int k=1;k<precision;k*=2)
  inverse=f4::multiply(b.bits,f4::multiply(inverse,inverse));
 return norm({f4::multiply(a.bits,inverse),a.base-b.base});
}
void erase(std::vector<S>&a,int i){a.erase(a.begin()+i);}
struct Result{std::vector<int>orders;std::array<int,6>counts{};};
Result profile(std::vector<S>d,std::vector<S>u,std::vector<S>v,int cutoff){
 Result out;bool arrow=false;S center;int original=d.size();
 auto record=[&](S pivot,int kind){
  if(!pivot.bits||pivot.base<0||pivot.base>=cutoff)throw std::runtime_error("unresolved structured pivot");
  if(!out.orders.empty()&&pivot.base<out.orders.back())throw std::runtime_error("nonmonotone structured pivot");
  out.orders.push_back(pivot.base);++out.counts[kind];
 };
 while(!d.empty()||arrow){
  int n=d.size(),best=1000000000,at=-1,other=-1,kind=-1;S pivot;
  auto consider=[&](S x,int k,int i,int j=-1){if(x.bits&&x.base<best){best=x.base;kind=k;at=i;other=j;pivot=x;}};
  if(!arrow){
   for(int i=0;i<n;++i)consider(add(d[i],mul(u[i],v[i])),0,i);
   int first=-1,second=-1;
   for(int j=0;j<n;++j)if(v[j].bits){
    if(first<0||v[j].base<v[first].base){second=first;first=j;}
    else if(second<0||v[j].base<v[second].base)second=j;
   }
   for(int i=0;i<n;++i)if(u[i].bits){int j=first==i?second:first;if(j>=0)consider(mul(u[i],v[j]),1,i,j);}
   if(kind<0)throw std::runtime_error("singular structured matrix");
   record(pivot,kind);
   if(kind==0){
    S factor=div(d[at],pivot);
    for(int k=0;k<n;++k)if(k!=at)u[k]=mul(u[k],factor);
    erase(d,at);erase(u,at);erase(v,at);
   }else{
    int i=at,j=other;
    center=add(add(div(mul(d[i],d[j]),pivot),div(mul(d[j],v[i]),v[j])),div(mul(d[i],u[j]),u[i]));
    S uf=div(d[i],u[i]),vf=div(d[j],v[j]);
    std::vector<S>nd,nu,nv;
    for(int k=0;k<n;++k)if(k!=i&&k!=j){nd.push_back(d[k]);nu.push_back(mul(u[k],uf));nv.push_back(mul(v[k],vf));}
    d=std::move(nd);u=std::move(nu);v=std::move(nv);arrow=true;
   }
  }else{
   consider(center,2,-1);
   for(int i=0;i<n;++i)consider(d[i],3,i);
   for(int i=0;i<n;++i){consider(u[i],4,i);consider(v[i],5,i);}
   if(kind<0)throw std::runtime_error("singular structured arrowhead");
   record(pivot,kind);
   if(kind==2){
    for(auto&entry:v)entry=div(entry,center);
    arrow=false;center={};
   }else if(kind==3){
    center=add(center,div(mul(u[at],v[at]),d[at]));
    erase(d,at);erase(u,at);erase(v,at);
   }else if(kind==4){
    S factor=div(d[at],u[at]);center=add(v[at],mul(center,factor));
    for(int k=0;k<n;++k)if(k!=at)u[k]=mul(u[k],factor);
    erase(d,at);erase(u,at);erase(v,at);
   }else{
    S factor=div(d[at],v[at]);center=add(u[at],mul(center,factor));
    for(int k=0;k<n;++k)if(k!=at)v[k]=mul(v[k],factor);
    erase(d,at);erase(u,at);erase(v,at);
   }
  }
 }
 if(int(out.orders.size())!=original)throw std::runtime_error("structured profile length");
 return out;
}
}

extern "C" int rank_one_f4_profile(int n,const std::uint64_t*coefficients,int input_words,
 int cutoff,std::int64_t*orders,std::int64_t*counts){
 try{
  if(n<1||n>4096||input_words<1||input_words>words||cutoff<1||cutoff>precision)return 1;
  std::vector<f4::Series>all(3*n);
  for(int i=0;i<3*n;++i){
   f4::Poly p;
   for(int k=0;k<input_words;++k){p.lo.w[k]=coefficients[(2*i)*input_words+k];p.hi.w[k]=coefficients[(2*i+1)*input_words+k];}
   all[i]=rank_one_detail::norm({p,0});
  }
  std::vector<f4::Series>d(all.begin(),all.begin()+n),u(all.begin()+n,all.begin()+2*n),v(all.begin()+2*n,all.end());
  auto result=rank_one_detail::profile(d,u,v,cutoff);
  for(int i=0;i<n;++i)orders[i]=result.orders[i];
  for(int i=0;i<6;++i)counts[i]=result.counts[i];
  return 0;
 }catch(...){return 4;}
}
