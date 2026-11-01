// Adaptive regular LDU cones for fixed Laurent matrices over F4.
// The legacy source supplies field arithmetic; its profile implementation
// is unchanged. Explicit coefficient knowledge guards every cached frame.
#include "native_rectangular_smith_f4.cpp"
#include <limits>
#include <memory>
#include <mutex>
#include <unordered_map>

namespace chamber {
using I=std::int64_t;
constexpr I INF=I(1)<<50;
I plus(I a,I b){return a>=INF||b>=INF?INF:std::min(INF,a+b);}
struct Series { f4::Poly bits; I base=0,known=INF; };
bool exact_zero(const Series&a){return !a.bits&&a.known==INF;}
I lower(const Series&a){return a.bits?a.base:a.known;}
Series normalize(Series a){
 if(a.bits){
  I width=a.known==INF?precision:std::min<I>(precision,a.known-a.base);
  if(width<=0)a.bits={};else if(width<precision)a.bits&=::low_mask(int(width));
  if(a.bits){int v=f4::valuation(a.bits);a.base+=v;a.bits=a.bits>>v;}
 }
 if(!a.bits)a.base=0;
 return a;
}
Series trim(Series a,I bound){
 if(a.bits&&a.base+f4::degree(a.bits)>=bound){a.known=std::min(a.known,bound);return normalize(a);}
 return a;
}
Series add(Series a,Series b){
 if(exact_zero(a))return b;if(exact_zero(b))return a;
 Series out;out.known=std::min(a.known,b.known);
 if(!a.bits&&!b.bits)return out;
 out.base=a.bits?(b.bits?std::min(a.base,b.base):a.base):b.base;
 I top=out.base;
 if(a.bits)top=std::max(top,a.base+f4::degree(a.bits));
 if(b.bits)top=std::max(top,b.base+f4::degree(b.bits));
 if(top>=out.base+precision)out.known=std::min(out.known,out.base+precision);
 if(a.bits&&a.base-out.base<precision)out.bits^=a.bits<<int(a.base-out.base);
 if(b.bits&&b.base-out.base<precision)out.bits^=b.bits<<int(b.base-out.base);
 return normalize(out);
}
Series multiply(const Series&a,const Series&b){
 if(exact_zero(a)||exact_zero(b))return {};
 Series out;out.known=std::min(plus(a.known,lower(b)),plus(b.known,lower(a)));
 if(a.bits&&b.bits){
  out.base=a.base+b.base;
  if(f4::degree(a.bits)+f4::degree(b.bits)>=precision)
   out.known=std::min(out.known,out.base+precision);
  out.bits=f4::multiply(a.bits,b.bits);
 }
 return normalize(out);
}
Series inverse_unit(const Series&a,int maximum=precision){
 if(!a.bits)throw std::runtime_error("unknown pivot");
 unsigned inverses[4]={0,1,3,2};
 f4::Poly inv(inverses[a.bits.constant()]);
 if(!inv)throw std::runtime_error("zero pivot unit");
 if(a.known==INF&&f4::degree(a.bits)==0)return {inv,0,INF};
 int width=int(std::min<I>(maximum,a.known==INF?precision:a.known-a.base));
 if(width<1)throw std::runtime_error("unresolved pivot unit");
 for(int p=1;p<width;p*=2)inv=f4::multiply(a.bits,f4::multiply(inv,inv));
 return normalize({inv,0,width});
}
Series shifted(Series a,I amount){
 if(a.bits)a.base+=amount;
 if(a.known!=INF)a.known+=amount;
 return a;
}

// Sign of sums of input doubles and scale times an exact small integer.
// A long-double filter handles separated cases; TwoSum/FMA expansions
// resolve boundaries without treating roundoff as a valid cone inequality.
int exact_sign(I base,double ri,double cj,double rp,double cp,double scale){
 if(!base&&ri==rp&&cj==cp)return 0;
 long double value=(long double)scale*base+(long double)ri-rp+(long double)cj-cp;
 long double magnitude=std::abs((long double)scale*base)+std::abs((long double)ri)
   +std::abs((long double)rp)+std::abs((long double)cj)+std::abs((long double)cp);
 long double error=32*std::numeric_limits<long double>::epsilon()*magnitude;
 if(value>error)return 1;if(value<-error)return -1;
 if(std::abs(base)>I(1)<<52)throw std::runtime_error("exponent comparison range");
 double product=scale*double(base);
 if(!std::isfinite(product))throw std::runtime_error("weight overflow");
 double terms[6]={product,std::fma(scale,double(base),-product),ri,-rp,cj,-cp};
 double expansion[16]={},next[16]={};int count=0;
 for(double term:terms){
  double q=term;int used=0;
  for(int i=0;i<count;++i){
   double x=q+expansion[i];
   double virtual_b=x-q;
   double y=(q-(x-virtual_b))+(expansion[i]-virtual_b);
   if(y)next[used++]=y;q=x;
  }
  if(q||!used)next[used++]=q;
  count=used;for(int i=0;i<count;++i)expansion[i]=next[i];
 }
 for(int i=count-1;i>=0;--i)if(expansion[i])return expansion[i]>0?1:-1;
 return 0;
}
struct Condition{bool row;int index,pivot;I difference;};
struct Pivot{I base;int row,column;};
struct Checkpoint{
 int done;std::size_t condition_count;
 std::vector<std::vector<Series>>residual;std::vector<int>rows,columns;
};
struct Chamber{
 std::vector<Condition>conditions;std::vector<Pivot>pivots;
 std::vector<std::size_t>step_offsets;std::shared_ptr<Checkpoint>checkpoint;
};
struct Oracle{
 int n,rank;std::vector<std::vector<Series>>source;
 std::array<std::vector<Chamber>,2>cones;
 std::uint64_t hits=0,builds=0,unresolved=0,reused_pivots=0;
};
int valid_prefix(const Chamber&cone,const double*rw,const double*cw,double scale){
 for(std::size_t index=0;index<cone.conditions.size();++index){
  const auto&test=cone.conditions[index];
  const double*weights=test.row?rw:cw;
  if(exact_sign(test.difference,weights[test.index],0,weights[test.pivot],0,scale)<0)
   return int(std::upper_bound(cone.step_offsets.begin(),cone.step_offsets.end(),index)
              -cone.step_offsets.begin())-1;
 }
 return int(cone.pivots.size());
}
void condition(Chamber&cone,bool row,int index,int pivot,const Series&a,I pivotbase){
 if(exact_zero(a))return;
 I bound=lower(a);
 if(bound==INF)return;
 cone.conditions.push_back({row,index,pivot,bound-pivotbase});
}
Chamber build(const Oracle&o,const double*rw,const double*cw,double scale,bool infinity,
              const Chamber*predecessor,int&reused){
 std::vector<std::vector<Series>>matrix;std::vector<int>rows,columns;Chamber cone;
 int done=0;
 if(predecessor){
  const auto&saved=*predecessor->checkpoint;done=saved.done;
  matrix=saved.residual;rows=saved.rows;columns=saved.columns;
  cone.conditions.assign(predecessor->conditions.begin(),predecessor->conditions.begin()+saved.condition_count);
  cone.pivots.assign(predecessor->pivots.begin(),predecessor->pivots.begin()+done);
  cone.step_offsets.assign(predecessor->step_offsets.begin(),predecessor->step_offsets.begin()+done);
  cone.checkpoint=predecessor->checkpoint;
 }else{
  matrix=o.source;rows.resize(o.n);columns.resize(o.n);
  for(int i=0;i<o.n;++i)rows[i]=columns[i]=i;
 }
 reused=done;int n=int(rows.size());
 if(!predecessor)for(auto&line:matrix)for(auto&entry:line){
  if(infinity&&entry.bits){
   int d=f4::degree(entry.bits);f4::Poly reversed;
   for(int k=0;k<=d;++k){
    if(entry.bits.lo.w[k/64]>>(k%64)&1)reversed.lo.w[(d-k)/64]|=std::uint64_t(1)<<((d-k)%64);
    if(entry.bits.hi.w[k/64]>>(k%64)&1)reversed.hi.w[(d-k)/64]|=std::uint64_t(1)<<((d-k)%64);
   }
   entry.bits=reversed;entry.base=-d-entry.base;
  }
  entry=normalize(entry);
 }
 double minimum=std::numeric_limits<double>::infinity();
 for(int i=0;i<n;++i)for(int j=0;j<n;++j)if(matrix[i][j].bits)
  minimum=std::min(minimum,scale*matrix[i][j].base+rw[rows[i]]+cw[columns[j]]);
 double cap=minimum+precision*scale;
 if(!std::isfinite(cap))throw std::runtime_error("unresolved finite coefficient window");
 std::vector<std::vector<I>>bounds(n,std::vector<I>(n));
 for(int i=0;i<n;++i)for(int j=0;j<n;++j){
  bounds[i][j]=I(std::floor((cap-rw[rows[i]]-cw[columns[j]])/scale))-1;
  matrix[i][j]=trim(matrix[i][j],bounds[i][j]);
 }
 int checkpoint_step=-1;
 for(int numerator:{1,2,3}){
  int candidate=o.rank*numerator/4;
  std::uint64_t dimension=o.n-candidate;
  if(candidate>0&&candidate<o.rank&&4*dimension*dimension*sizeof(Series)<=64*1024*1024){
   checkpoint_step=candidate;break;
  }
 }
 for(int k=0;k<o.rank-done;++k){
  int global_step=k+done;
  if(global_step==checkpoint_step&&!cone.checkpoint){
   auto saved=std::make_shared<Checkpoint>();saved->done=global_step;
   saved->condition_count=cone.conditions.size();
   saved->rows.assign(rows.begin()+k,rows.end());saved->columns.assign(columns.begin()+k,columns.end());
   int dimension=n-k;saved->residual.resize(dimension,std::vector<Series>(dimension));
   for(int i=0;i<dimension;++i)for(int j=0;j<dimension;++j)saved->residual[i][j]=matrix[k+i][k+j];
   cone.checkpoint=std::move(saved);
  }
  cone.step_offsets.push_back(cone.conditions.size());
  int pi=-1,pj=-1;
  if(predecessor&&global_step<int(predecessor->pivots.size())){
   const auto&hint=predecessor->pivots[global_step];
   for(int i=k;i<n;++i)if(rows[i]==hint.row)pi=i;
   for(int j=k;j<n;++j)if(columns[j]==hint.column)pj=j;
   if(pi<0||pj<0||!matrix[pi][pj].bits)pi=pj=-1;
  }
  if(pi<0)for(int i=k;i<n&&pi<0;++i)for(int j=k;j<n;++j)if(matrix[i][j].bits){pi=i;pj=j;break;}
  if(pi<0)throw std::runtime_error("finite knowledge did not resolve every pivot");
  // A local row/column minimum suffices: integral row and column clearing
  // diagonalizes the matrix, and the final diagonal is sorted afterwards.
  // This avoids an unnecessary cubic number of exact pivot comparisons.
  bool changed=true;
  while(changed){
   changed=false;
   for(int i=k;i<n;++i)if(matrix[i][pj].bits&&exact_sign(
      matrix[i][pj].base-matrix[pi][pj].base,rw[rows[i]],0,rw[rows[pi]],0,scale)<0){pi=i;changed=true;}
   for(int j=k;j<n;++j)if(matrix[pi][j].bits&&exact_sign(
      matrix[pi][j].base-matrix[pi][pj].base,0,cw[columns[j]],0,cw[columns[pj]],scale)<0){pj=j;changed=true;}
  }
  std::swap(matrix[k],matrix[pi]);std::swap(bounds[k],bounds[pi]);std::swap(rows[k],rows[pi]);
  for(auto&line:matrix)std::swap(line[k],line[pj]);
  for(auto&line:bounds)std::swap(line[k],line[pj]);std::swap(columns[k],columns[pj]);
  I pivotbase=matrix[k][k].base;
  cone.pivots.push_back({pivotbase,rows[k],columns[k]});
  int width=int(std::clamp<I>(bounds[k][k]-pivotbase,1,precision));
  auto inverse=inverse_unit(matrix[k][k],width);
  for(int j=k+1;j<n;++j)matrix[k][j]=trim(multiply(matrix[k][j],inverse),bounds[k][j]);
  // Division by the actual pivot unit makes the pivot exactly monomial.
  matrix[k][k]={f4::Poly(1),pivotbase,INF};
  for(int i=k+1;i<n;++i)condition(cone,true,rows[i],rows[k],matrix[i][k],pivotbase);
  for(int j=k+1;j<n;++j)condition(cone,false,columns[j],columns[k],matrix[k][j],pivotbase);
  for(int i=k+1;i<n;++i){
   if(exact_zero(matrix[i][k]))continue;
   auto factor=shifted(matrix[i][k],-pivotbase);
   for(int j=k+1;j<n;++j)if(!exact_zero(matrix[k][j]))
    matrix[i][j]=trim(add(matrix[i][j],multiply(factor,matrix[k][j])),bounds[i][j]);
   matrix[i][k]={};
  }
 }
 return cone;
}
void output(const Chamber&cone,const double*rw,const double*cw,double scale,
            double*orders,I*prefix,std::uint32_t*rows,std::uint32_t*columns){
 auto pivots=cone.pivots;
 std::stable_sort(pivots.begin(),pivots.end(),[&](const Pivot&a,const Pivot&b){
  return exact_sign(a.base-b.base,rw[a.row],cw[a.column],rw[b.row],cw[b.column],scale)<0;
 });
 I cumulative=0;
 for(std::size_t i=0;i<pivots.size();++i){
  const auto&p=pivots[i];orders[i]=scale*p.base+rw[p.row]+cw[p.column];
  if(!std::isfinite(orders[i]))throw std::runtime_error("profile overflow");
  cumulative+=p.base;prefix[i]=cumulative;rows[i]=p.row;columns[i]=p.column;
 }
}
std::mutex mutex;
std::uint64_t next_handle=1;
std::unordered_map<std::uint64_t,std::unique_ptr<Oracle>>objects;
}

extern "C" std::uint64_t smith_chamber_create(int n,int rank,const std::uint64_t*coefficients,
 int words_per_value,const std::int64_t*shifts){
 std::lock_guard<std::mutex>lock(chamber::mutex);
 try{
  if(n<1||n>1024||rank<1||rank>n||words_per_value<1||words_per_value>words)return 0;
  auto object=std::make_unique<chamber::Oracle>();object->n=n;object->rank=rank;
  object->source.resize(n,std::vector<chamber::Series>(n));
  for(int i=0;i<n;++i)for(int j=0;j<n;++j){
   int index=i*n+j;auto&entry=object->source[i][j];entry.base=shifts[index];
   if(std::abs(entry.base)>1000000)return 0;
   for(int k=0;k<words_per_value;++k){
    entry.bits.lo.w[k]=coefficients[(2*index)*words_per_value+k];
    entry.bits.hi.w[k]=coefficients[(2*index+1)*words_per_value+k];
   }
  }
  auto handle=chamber::next_handle++;chamber::objects.emplace(handle,std::move(object));return handle;
 }catch(...){return 0;}
}
extern "C" void smith_chamber_destroy(std::uint64_t handle){
 std::lock_guard<std::mutex>lock(chamber::mutex);chamber::objects.erase(handle);
}
extern "C" int smith_chamber_compare(std::int64_t base,double ri,double cj,
 double rp,double cp,double scale){
 try{return chamber::exact_sign(base,ri,cj,rp,cp,scale);}catch(...){return -2;}
}
extern "C" int smith_chamber_profile(std::uint64_t handle,const double*rw,const double*cw,
 int infinity,double scale,double*orders,std::int64_t*prefix,std::uint32_t*rows,
 std::uint32_t*columns,std::uint64_t*stats){
 std::lock_guard<std::mutex>lock(chamber::mutex);
 try{
  auto found=chamber::objects.find(handle);if(found==chamber::objects.end()||scale<=0)return 1;
  auto&o=*found->second;auto&cones=o.cones[bool(infinity)];
  const chamber::Chamber*predecessor=nullptr;int best_prefix=-1;
  for(auto current=cones.rbegin();current!=cones.rend();++current){
   int allowed=chamber::valid_prefix(*current,rw,cw,scale);
   if(allowed==o.rank){
    ++o.hits;chamber::output(*current,rw,cw,scale,orders,prefix,rows,columns);
    stats[0]=1;stats[1]=o.hits;stats[2]=o.builds;stats[3]=o.unresolved;stats[4]=current->conditions.size();
    stats[5]=0;stats[6]=o.reused_pivots;stats[7]=o.rank;return 0;
   }
   if(current->checkpoint&&allowed>=current->checkpoint->done&&allowed>best_prefix){
    predecessor=&*current;best_prefix=allowed;
   }
  }
  ++o.builds;
  try{
   int reused=0;auto cone=chamber::build(o,rw,cw,scale,bool(infinity),predecessor,reused);
   bool allowed=chamber::valid_prefix(cone,rw,cw,scale)==o.rank;o.reused_pivots+=reused;
   stats[5]=reused;stats[6]=o.reused_pivots;stats[7]=best_prefix<0?0:best_prefix;
   if(cones.size()==4)cones.erase(cones.begin());cones.push_back(std::move(cone));
   if(!allowed){++o.unresolved;stats[1]=o.hits;stats[2]=o.builds;stats[3]=o.unresolved;return 3;}
   chamber::output(cones.back(),rw,cw,scale,orders,prefix,rows,columns);
   stats[0]=0;stats[1]=o.hits;stats[2]=o.builds;stats[3]=o.unresolved;stats[4]=cones.back().conditions.size();return 0;
  }catch(...){++o.unresolved;stats[1]=o.hits;stats[2]=o.builds;stats[3]=o.unresolved;return 4;}
 }catch(...){return 5;}
}
