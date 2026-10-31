// Exact coordinate export into the two conserved Fock parity blocks.
// This separate entry point permits eleven active modes without changing
// the established ten-mode interfaces in native_fock_f4.cpp.
#include "native_fock_f4.cpp"

extern "C" int fock_f4_coordinates_parity(int k,int count,const int*program,
 int coefficient_count,const std::uint64_t*coefficients,int car_count,
 const int*car_offsets,const int*car_terms,const int*diagonal_offsets,
 int left_multiplicity,const int*vacuum_rows,int right_multiplicity,
 const int*vacuum_columns,const int*output_registers,int common_base,int output_words,
 const int*row_parities,const int*row_positions,const int*column_parities,const int*column_positions,
 std::uint64_t*output,std::int64_t*shifts,std::int64_t*metadata){
 try{
  using namespace fock;
  if(k<1||k>11||count<1||count>30000||output_words<2||output_words>words)return 1;
  int n=1<<k,half=n/2;
  if(left_multiplicity<1||right_multiplicity<1||n%left_multiplicity||n%right_multiplicity)return 1;
  int left_states=n/left_multiplicity,right_states=n/right_multiplicity;
  std::vector<bool>seen_rows(n,false),seen_columns(n,false);
  for(int i=0;i<n;++i){
   if(row_parities[i]<0||row_parities[i]>1||column_parities[i]<0||column_parities[i]>1
      ||row_positions[i]<0||row_positions[i]>=half||column_positions[i]<0||column_positions[i]>=half)return 1;
   int a=row_parities[i]*half+row_positions[i],b=column_parities[i]*half+column_positions[i];
   if(seen_rows[a]||seen_columns[b])return 1;seen_rows[a]=seen_columns[b]=true;
  }
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
     int row=left_state*left_multiplicity+row_copy,column=right_state*right_multiplicity+column_copy;
     auto value=reg[regid][source_row];
     if(row_parities[row]!=column_parities[column]){if(value)return 6;continue;}
     int index=row_parities[row]*half*half+row_positions[row]*half+column_positions[column];
     if(!value){shifts[index]=0;continue;}
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
