// Structured Smith elimination with certified absolute series error bounds.
#include "native_rectangular_smith_f4.cpp"

namespace certified_rank_one {
constexpr int infinity=1000000000;
struct S {f4::Poly bits;int base=0,error=infinity;};
int plus(int e,std::int64_t shift){
 if(e==infinity)return infinity;
 auto result=std::int64_t(e)+shift;
 if(result<=-infinity||result>=infinity)throw std::runtime_error("series exponent outside certified integer range");
 return int(result);
}
bool exact_zero(const S&a){return !a.bits&&a.error==infinity;}
int lower(const S&a){return a.bits?a.base:a.error;}
S norm(S a){
 if(a.error<infinity&&a.bits){
  int known=a.error-a.base;
  if(known<=0)a.bits={};else if(known<precision)a.bits&=::low_mask(known);
 }
 if(!a.bits){a.base=a.error<infinity?a.error:0;return a;}
 int v=f4::valuation(a.bits);a.base+=v;a.bits=a.bits>>v;
 if(a.base>=a.error)throw std::runtime_error("uncertified stored coefficient");
 return a;
}
S add(S a,S b){
 if(exact_zero(a))return b;if(exact_zero(b))return a;
 if(a.base>b.base)std::swap(a,b);
 int error=std::min(a.error,b.error),shift=b.base-a.base;
 if(b.bits&&std::int64_t(f4::degree(b.bits))+shift>=precision)
  error=std::min(error,plus(a.base,precision));
 return norm({a.bits^(b.bits<<shift),a.base,error});
}
S mul(S a,S b){
 if(exact_zero(a)||exact_zero(b))return {};
 int va=lower(a),vb=lower(b),base=plus(va,vb);
 int error=std::min(plus(a.error,vb),plus(b.error,va));
 if(a.bits&&b.bits&&f4::degree(a.bits)+f4::degree(b.bits)>=precision)
  error=std::min(error,plus(base,precision));
 return norm({f4::multiply(a.bits,b.bits),base,error});
}
S div(S a,S b){
 if(!b.bits||b.base>=b.error)throw std::runtime_error("divisor leading term unknown");
 if(exact_zero(a))return {};
 if(!a.bits)return norm({{},plus(a.error,-b.base),plus(a.error,-b.base)});
 int base=plus(a.base,-std::int64_t(b.base));
 int error=std::min(plus(a.error,-std::int64_t(b.base)),plus(b.error,std::int64_t(a.base)-2*std::int64_t(b.base)));
 unsigned inverse_scalar[4]={0,1,3,2};f4::Poly inverse(inverse_scalar[b.bits.constant()]);
 if(f4::degree(b.bits)>0){
  for(int k=1;k<precision;k*=2)inverse=f4::multiply(b.bits,f4::multiply(inverse,inverse));
  error=std::min(error,plus(base,precision));
 }
 return norm({f4::multiply(a.bits,inverse),base,error});
}
void erase(std::vector<S>&a,int i){a.erase(a.begin()+i);}
struct Result{std::vector<int>orders;std::array<int,6>counts{};int minimum_error_margin=infinity;};
Result profile(std::vector<S>d,std::vector<S>u,std::vector<S>v,int cutoff){
 Result out;bool arrow=false;S center;int original=d.size();
 while(!d.empty()||arrow){
  int n=d.size(),best=infinity,certified_lower=infinity,at=-1,other=-1,kind=-1;S pivot;
  auto consider=[&](S x,int k,int i,int j=-1){
   certified_lower=std::min(certified_lower,lower(x));
   if(x.bits&&x.base<best){best=x.base;kind=k;at=i;other=j;pivot=x;}
  };
  if(!arrow){
   for(int i=0;i<n;++i)consider(add(d[i],mul(u[i],v[i])),0,i);
   int first=-1,second=-1;
   for(int j=0;j<n;++j)if(!exact_zero(v[j])){
    if(first<0||lower(v[j])<lower(v[first])){second=first;first=j;}
    else if(second<0||lower(v[j])<lower(v[second]))second=j;
   }
   // A minimum lower bound covers every off-diagonal entry. Every known
   // off-diagonal minimum is also seen unless an unresolved lower bound
   // is strictly smaller, in which case the pivot gate correctly fails.
   for(int i=0;i<n;++i)if(!exact_zero(u[i])){
    int j=first==i?second:first;if(j>=0)consider(mul(u[i],v[j]),1,i,j);
   }
  }else{
   consider(center,2,-1);
   for(int i=0;i<n;++i)consider(d[i],3,i);
   for(int i=0;i<n;++i){consider(u[i],4,i);consider(v[i],5,i);}
  }
  if(kind<0||best>certified_lower||best<0||best>=cutoff)
   throw std::runtime_error("pivot not certified inside output window");
  if(pivot.error<=best)throw std::runtime_error("pivot leading term not certified");
  if(!out.orders.empty()&&best<out.orders.back())throw std::runtime_error("nonmonotone pivot");
  if(pivot.error<infinity)out.minimum_error_margin=std::min(out.minimum_error_margin,pivot.error-best);
  out.orders.push_back(best);++out.counts[kind];
  if(!arrow){
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
  }else if(kind==2){
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
 if(int(out.orders.size())!=original)throw std::runtime_error("incomplete profile");
 return out;
}
}

extern "C" int certified_rank_one_f4_profile(int n,const std::uint64_t*coefficients,int input_words,
 int cutoff,std::int64_t*orders,std::int64_t*counts,std::int64_t*margin){
 try{
  using namespace certified_rank_one;
  if(n<1||n>4096||input_words<1||input_words>words||cutoff<1||cutoff>precision)return 1;
  std::vector<S>all(3*n);
  for(int i=0;i<3*n;++i){
   f4::Poly p;
   for(int k=0;k<input_words;++k){p.lo.w[k]=coefficients[(2*i)*input_words+k];p.hi.w[k]=coefficients[(2*i+1)*input_words+k];}
   all[i]=norm({p,0,infinity});
  }
  std::vector<S>d(all.begin(),all.begin()+n),u(all.begin()+n,all.begin()+2*n),v(all.begin()+2*n,all.end());
  auto result=profile(d,u,v,cutoff);
  for(int i=0;i<n;++i)orders[i]=result.orders[i];
  for(int i=0;i<6;++i)counts[i]=result.counts[i];
  *margin=result.minimum_error_margin;
  return 0;
 }catch(...){return 4;}
}
