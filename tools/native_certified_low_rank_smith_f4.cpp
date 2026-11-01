// Consecutive rank-one diagonalizations, propagating all remaining probes.
#include "native_certified_rank_one_smith_f4.cpp"

namespace certified_low_rank {
using namespace certified_rank_one;
struct Probe {
 std::vector<S>left,right,outleft,outright;
 S left_center,right_center;
};
struct PhaseResult {
 std::vector<int>orders;
 std::vector<S>diagonal;
 std::array<int,6>counts{};
 int minimum_error_margin=infinity;
};

S divide(S a,S b){
 if(!b.bits||b.base>=b.error)throw std::runtime_error("divisor leading term unknown");
 if(exact_zero(a))return {};
 if(!a.bits)return norm({{},plus(a.error,-b.base),plus(a.error,-b.base)});
 int base=plus(a.base,-std::int64_t(b.base));
 int error=std::min(plus(a.error,-std::int64_t(b.base)),plus(b.error,std::int64_t(a.base)-2*std::int64_t(b.base)));
 static thread_local f4::Poly cached_unit,cached_inverse;
 if(cached_unit!=b.bits){
  unsigned inverse_scalar[4]={0,1,3,2};cached_inverse=f4::Poly(inverse_scalar[b.bits.constant()]);
  if(f4::degree(b.bits)>0)for(int k=1;k<precision;k*=2)
   cached_inverse=f4::multiply(b.bits,f4::multiply(cached_inverse,cached_inverse));
  cached_unit=b.bits;
 }
 if(f4::degree(b.bits)>0)error=std::min(error,plus(base,precision));
 return norm({f4::multiply(a.bits,cached_inverse),base,error});
}
struct PivotMaps {std::vector<S>column,row;int pivot_row,pivot_column;};
PivotMaps pivot_maps(int kind,int at,int other,S pivot,const std::vector<S>&d,
 const std::vector<S>&u,const std::vector<S>&v,S center){
 int n=d.size();bool arrow=kind>=2;int size=n+int(arrow);
 int i=kind==2||kind==5?n:at,j=kind==2||kind==4?n:(kind==1?other:at);
 PivotMaps maps;maps.column.resize(size);maps.row.resize(size);maps.pivot_row=i;maps.pivot_column=j;
 auto entry=[&](int row,int col){
  if(!arrow){auto value=mul(u[row],v[col]);return row==col?add(d[row],value):value;}
  if(row==n&&col==n)return center;
  if(row==n)return v[col];if(col==n)return u[row];
  return row==col?d[row]:S{};
 };
 for(int k=0;k<size;++k){
  if(k!=i)maps.column[k]=divide(entry(k,j),pivot);
  if(k!=j)maps.row[k]=divide(entry(i,k),pivot);
 }
 return maps;
}
void propagate(Probe&probe,int kind,int at,int other,const PivotMaps&maps){
 int n=probe.left.size();auto left=std::move(probe.left),right=std::move(probe.right);
 probe.left.clear();probe.right.clear();
 if(kind>=2){left.push_back(probe.left_center);right.push_back(probe.right_center);}
 S lp=left[maps.pivot_row],rp=right[maps.pivot_column];
 probe.outleft.push_back(lp);probe.outright.push_back(rp);
 for(int k=0;k<int(left.size());++k){
  if(k!=maps.pivot_row)left[k]=add(left[k],mul(maps.column[k],lp));
  if(k!=maps.pivot_column)right[k]=add(right[k],mul(maps.row[k],rp));
 }
 if(kind==0){erase(left,at);erase(right,at);probe.left=std::move(left);probe.right=std::move(right);}
 else if(kind==1){
  probe.left_center=left[other];probe.right_center=right[at];
  for(int k=0;k<n;++k)if(k!=at&&k!=other){probe.left.push_back(left[k]);probe.right.push_back(right[k]);}
 }else if(kind==2){
  left.pop_back();right.pop_back();probe.left=std::move(left);probe.right=std::move(right);
  probe.left_center={};probe.right_center={};
 }else{
  probe.left_center=left[kind==5?at:n];probe.right_center=right[kind==4?at:n];
  for(int k=0;k<n;++k)if(k!=at){probe.left.push_back(left[k]);probe.right.push_back(right[k]);}
 }
}

PhaseResult phase(std::vector<S>d,std::vector<S>u,std::vector<S>v,int cutoff,std::vector<Probe>&probes){
 PhaseResult out;bool arrow=false;S center;int original=d.size();
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
  out.orders.push_back(best);out.diagonal.push_back(pivot);++out.counts[kind];
  if(!probes.empty()){auto maps=pivot_maps(kind,at,other,pivot,d,u,v,center);
   for(auto&probe:probes)propagate(probe,kind,at,other,maps);}
  if(!arrow){
   if(kind==0){
    S factor=divide(d[at],pivot);
    for(int k=0;k<n;++k)if(k!=at)u[k]=mul(u[k],factor);
    erase(d,at);erase(u,at);erase(v,at);
   }else{
    int i=at,j=other;
    center=add(add(divide(mul(d[i],d[j]),pivot),divide(mul(d[j],v[i]),v[j])),divide(mul(d[i],u[j]),u[i]));
    S uf=divide(d[i],u[i]),vf=divide(d[j],v[j]);
    std::vector<S>nd,nu,nv;
    for(int k=0;k<n;++k)if(k!=i&&k!=j){nd.push_back(d[k]);nu.push_back(mul(u[k],uf));nv.push_back(mul(v[k],vf));}
    d=std::move(nd);u=std::move(nu);v=std::move(nv);arrow=true;
   }
  }else if(kind==2){
   for(auto&entry:v)entry=divide(entry,center);
   arrow=false;center={};
  }else if(kind==3){
   center=add(center,divide(mul(u[at],v[at]),d[at]));
   erase(d,at);erase(u,at);erase(v,at);
  }else if(kind==4){
   S factor=divide(d[at],u[at]);center=add(v[at],mul(center,factor));
   for(int k=0;k<n;++k)if(k!=at)u[k]=mul(u[k],factor);
   erase(d,at);erase(u,at);erase(v,at);
  }else{
   S factor=divide(d[at],v[at]);center=add(u[at],mul(center,factor));
   for(int k=0;k<n;++k)if(k!=at)v[k]=mul(v[k],factor);
   erase(d,at);erase(u,at);erase(v,at);
  }
 }
 if(int(out.orders.size())!=original)throw std::runtime_error("incomplete profile");
 return out;
}
}

extern "C" int certified_low_rank_f4_profile(int n,int rank,const std::uint64_t*coefficients,int input_words,
 int cutoff,std::int64_t*orders,std::int64_t*counts,std::int64_t*margin){
 try{
  using namespace certified_low_rank;
  if(n<1||n>4096||rank<1||rank>32||input_words<1||input_words>words||cutoff<1||cutoff>precision)return 1;
  std::vector<S>all((1+2*rank)*n);
  for(int i=0;i<(1+2*rank)*n;++i){
   f4::Poly p;
   for(int k=0;k<input_words;++k){p.lo.w[k]=coefficients[(2*i)*input_words+k];p.hi.w[k]=coefficients[(2*i+1)*input_words+k];}
   all[i]=norm({p,0,infinity});
  }
  std::vector<S>d(all.begin(),all.begin()+n);
  std::vector<std::vector<S>>left(rank),right(rank);
  for(int k=0;k<rank;++k){
   left[k]=std::vector<S>(all.begin()+(1+2*k)*n,all.begin()+(2+2*k)*n);
   right[k]=std::vector<S>(all.begin()+(2+2*k)*n,all.begin()+(3+2*k)*n);
  }
  std::array<int,6>total{},probe_total{};int minimum=infinity;
  for(int k=0;k<rank;++k){
   std::vector<Probe>probes;
   for(int j=k+1;j<rank;++j){Probe p;p.left=left[j];p.right=right[j];probes.push_back(std::move(p));}
   auto result=phase(d,left[k],right[k],k+1==rank?cutoff:precision,probes);
   d=result.diagonal;minimum=std::min(minimum,result.minimum_error_margin);
   for(int i=0;i<6;++i){total[i]+=result.counts[i];if(k+1<rank)probe_total[i]+=result.counts[i];}
   for(int j=k+1;j<rank;++j){
    left[j]=std::move(probes[j-k-1].outleft);right[j]=std::move(probes[j-k-1].outright);
    if(int(left[j].size())!=n||int(right[j].size())!=n)throw std::runtime_error("incomplete transformed probe");
   }
   if(k+1==rank)for(int i=0;i<n;++i)orders[i]=result.orders[i];
  }
  for(int i=0;i<6;++i){counts[i]=total[i];counts[6+i]=probe_total[i];}*margin=minimum;
  return 0;
 }catch(...){return 4;}
}
