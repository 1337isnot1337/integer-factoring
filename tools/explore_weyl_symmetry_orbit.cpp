#define main rank_one_orbit_base_main
#include "explore_weyl_rank_one_orbit.cpp"
#undef main

int main(int argc,char** argv){
 if(argc!=2)return 2;
 auto started=std::chrono::steady_clock::now();Data data;
 std::string prefix=argv[1];read_dump(data,prefix+".row.bin",true);read_dump(data,prefix+".col.bin",false);
 for(int i=0;i<2;++i){auto reverse=data.edges[i];for(auto& e:reverse)std::swap(e.row,e.col);data.edges.push_back(std::move(reverse));}
 auto closure=orbit(data,false);int gram=pairing_rank(closure,closure);
 std::cout<<"{\"energy\":"<<data.h<<",\"labels\":"<<data.d<<",\"parity_zero_dimension\":"<<data.n
  <<",\"transpose_closed_orbit\":"<<closure.basis.size()<<",\"pairing_rank\":"<<gram
  <<",\"maximum_word_length\":"<<closure.maxdepth
  <<",\"seconds\":"<<std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count()<<"}\n";
 return 0;
}
