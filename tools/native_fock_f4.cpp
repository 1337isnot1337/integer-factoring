// Exact truncated polynomial evaluation of a supplied fermionic program.
// The Python frontend clears only polynomial units at zero, tracks all
// Laurent shifts, and checks the common precision before Smith reduction.
#ifndef F4_SERIES_WORDS
#define F4_SERIES_WORDS 16
#endif
#include "native_rectangular_smith_f4.cpp"

namespace fock {
using P=f4::Poly;
P mul(const P&a,const P&b){
 if(!a||!b)return {};
 if(a==P(1))return b;
 if(a==P(2)){P r;r.lo=b.hi;r.hi=b.lo^b.hi;return r;}
 if(a==P(3)){P r;r.lo=b.lo^b.hi;r.hi=b.lo;return r;}
 return f4::multiply(a,b);
}
}

// Instruction records: kind, source1, source2, coefficient1, coefficient2.
// Kind 0 is a basis vector, kind 1 a CAR sum (coefficient1=CAR id), and
// kind 2 a polynomial linear combination of two earlier vectors.
// CAR records: mode, creation flag, polynomial coefficient id.
extern "C" int fock_f4_evaluate_v2(int k,int count,const int*program,
 int coefficient_count,const std::uint64_t*coefficients,int car_count,
 const int*car_offsets,const int*car_terms,int final_register,
 const int*diagonal_orders,int common_base,int output_words,std::uint64_t*output,
 std::int64_t*shifts,std::int64_t*metadata){
 try{
  using namespace fock;
  if(k<1||k>10||count<1||count>10000||final_register<0||final_register>=count||output_words<2||output_words>words)return 1;
  int n=1<<k;
  std::vector<P>cs(coefficient_count);
  for(int j=0;j<coefficient_count;++j)for(int w=0;w<words;++w){
   cs[j].lo.w[w]=coefficients[(2*j)*words+w];
   cs[j].hi.w[w]=coefficients[(2*j+1)*words+w];
  }
  std::vector<int>last(count,-1);
  for(int j=0;j<count;++j){auto p=program+5*j;
   if(p[0]==0){if(j!=0)return 2;continue;}
   if(p[1]<0||p[1]>=j)return 2;last[p[1]]=j;
   if(p[0]==2){if(p[2]<0||p[2]>=j)return 2;last[p[2]]=j;}
   if(p[0]<1||p[0]>2)return 2;
  }
  last[final_register]=count;
  std::int64_t smallest=1000000000,known=1000000000,maximum_degree=-1;
  for(int column=0;column<n;++column){
   std::vector<std::vector<P>>reg(count);
   for(int j=0;j<count;++j){
    auto p=program+5*j;auto&dst=reg[j];dst.resize(n);
    if(p[0]==0)dst[column]=P(1);
    else if(p[0]==1){
     int op=p[3];if(op<0||op>=car_count)return 3;
     const auto&src=reg[p[1]];
     for(int t=car_offsets[op];t<car_offsets[op+1];++t){
      int mode=car_terms[3*t],creation=car_terms[3*t+1],co=car_terms[3*t+2];
      if(mode<0||mode>=k||co<0||co>=coefficient_count)return 3;
      int bit=1<<mode;
      for(int a=0;a<n;++a)if(bool(a&bit)!=bool(creation)&&src[a])
       dst[a^bit]^=mul(cs[co],src[a]);
     }
    }else{
     if(p[3]<0||p[3]>=coefficient_count||p[4]<0||p[4]>=coefficient_count)return 3;
     const auto&a=reg[p[1]];const auto&b=reg[p[2]];
     for(int i=0;i<n;++i)dst[i]=mul(cs[p[3]],a[i])^mul(cs[p[4]],b[i]);
    }
    if(p[0]&&last[p[1]]==j)std::vector<P>().swap(reg[p[1]]);
    if(p[0]==2&&p[2]!=p[1]&&last[p[2]]==j)std::vector<P>().swap(reg[p[2]]);
   }
   int shift=common_base+diagonal_orders[column];
   known=std::min(known,std::int64_t(shift+precision));
   for(int row=0;row<n;++row){
    auto value=reg[final_register][row];int index=row*n+column;
    if(!value){shifts[index]=0;continue;}
    int v=f4::valuation(value);value=value>>v;shifts[index]=shift+v;
    maximum_degree=std::max(maximum_degree,std::int64_t(f4::degree(value)));
    smallest=std::min(smallest,shifts[index]);
    for(int w=0;w<output_words;++w){output[(2*index)*output_words+w]=value.lo.w[w];output[(2*index+1)*output_words+w]=value.hi.w[w];}
   }
  }
  metadata[0]=smallest;metadata[1]=known;metadata[2]=precision;metadata[3]=maximum_degree;
  if(smallest==1000000000||known<smallest+64*output_words)return 4;
  return 0;
 }catch(...){return 5;}
}

// The coordinate version evaluates B on selected vacuum input vectors,
// inserts Gamma(Sigma) explicitly (instruction kind 3), and reads L on
// selected vacuum output coordinates. Every output register has the same
// tracked Laurent shift and unit denominator, supplied by the frontend.
extern "C" int fock_f4_coordinates_v2(int k,int count,const int*program,
 int coefficient_count,const std::uint64_t*coefficients,int car_count,
 const int*car_offsets,const int*car_terms,const int*diagonal_offsets,
 int left_multiplicity,const int*vacuum_rows,int right_multiplicity,
 const int*vacuum_columns,const int*output_registers,int common_base,int output_words,
 std::uint64_t*output,std::int64_t*shifts,std::int64_t*metadata){
 try{
  using namespace fock;
  if(k<1||k>10||count<1||count>30000||output_words<2||output_words>words)return 1;
  int n=1<<k;
  if(left_multiplicity<1||right_multiplicity<1||n%left_multiplicity||n%right_multiplicity)return 1;
  int left_states=n/left_multiplicity,right_states=n/right_multiplicity;
  std::vector<P>cs(coefficient_count);
  for(int j=0;j<coefficient_count;++j)for(int w=0;w<words;++w){
   cs[j].lo.w[w]=coefficients[(2*j)*words+w];cs[j].hi.w[w]=coefficients[(2*j+1)*words+w];
  }
  std::vector<int>last(count,-1);
  for(int j=0;j<count;++j){auto p=program+5*j;
   if(p[0]==0){if(j!=0)return 2;continue;}
   if(p[1]<0||p[1]>=j||p[0]<1||p[0]>3)return 2;last[p[1]]=j;
   if(p[0]==2){if(p[2]<0||p[2]>=j)return 2;last[p[2]]=j;}
  }
  for(int j=0;j<left_states*right_states;++j){int r=output_registers[j];if(r<0||r>=count)return 2;last[r]=count;}
  std::int64_t smallest=1000000000,maximum_degree=-1;
  for(int column_copy=0;column_copy<right_multiplicity;++column_copy){
   std::vector<std::vector<P>>reg(count);
   for(int j=0;j<count;++j){auto p=program+5*j;auto&dst=reg[j];dst.resize(n);
    if(p[0]==0){if(vacuum_columns[column_copy]<0||vacuum_columns[column_copy]>=n)return 2;dst[vacuum_columns[column_copy]]=P(1);}
    else if(p[0]==1){
     int op=p[3];if(op<0||op>=car_count)return 3;const auto&src=reg[p[1]];
     for(int t=car_offsets[op];t<car_offsets[op+1];++t){
      int mode=car_terms[3*t],creation=car_terms[3*t+1],co=car_terms[3*t+2];
      if(mode<0||mode>=k||co<0||co>=coefficient_count)return 3;int bit=1<<mode;
      for(int a=0;a<n;++a)if(bool(a&bit)!=bool(creation)&&src[a])dst[a^bit]^=mul(cs[co],src[a]);
     }
    }else if(p[0]==2){
     if(p[3]<0||p[3]>=coefficient_count||p[4]<0||p[4]>=coefficient_count)return 3;
     const auto&a=reg[p[1]];const auto&b=reg[p[2]];
     for(int i=0;i<n;++i)dst[i]=mul(cs[p[3]],a[i])^mul(cs[p[4]],b[i]);
    }else{
     const auto&src=reg[p[1]];
     for(int i=0;i<n;++i){if(diagonal_offsets[i]<0)return 3;dst[i]=src[i]<<diagonal_offsets[i];}
    }
    if(p[0]&&last[p[1]]==j)std::vector<P>().swap(reg[p[1]]);
    if(p[0]==2&&p[2]!=p[1]&&last[p[2]]==j)std::vector<P>().swap(reg[p[2]]);
   }
   for(int right_state=0;right_state<right_states;++right_state)for(int left_state=0;left_state<left_states;++left_state){
    int regid=output_registers[left_state*right_states+right_state];
    for(int row_copy=0;row_copy<left_multiplicity;++row_copy){
     int source_row=vacuum_rows[row_copy];if(source_row<0||source_row>=n)return 2;
     int row=left_state*left_multiplicity+row_copy,column=right_state*right_multiplicity+column_copy,index=row*n+column;
     auto value=reg[regid][source_row];if(!value){shifts[index]=0;continue;}
     int v=f4::valuation(value);value=value>>v;shifts[index]=common_base+v;
     smallest=std::min(smallest,shifts[index]);maximum_degree=std::max(maximum_degree,std::int64_t(f4::degree(value)));
     for(int w=0;w<output_words;++w){output[(2*index)*output_words+w]=value.lo.w[w];output[(2*index+1)*output_words+w]=value.hi.w[w];}
    }
   }
  }
  metadata[0]=smallest;metadata[1]=common_base+precision;metadata[2]=precision;metadata[3]=maximum_degree;
  if(smallest==1000000000||metadata[1]<smallest+64*output_words)return 4;
  return 0;
 }catch(...){return 5;}
}
