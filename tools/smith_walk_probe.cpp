#define SMITH_WALK_NO_MAIN
#include "smith_walk_search.cpp"

// Input: 64 polynomial bitmasks in label-major, auxiliary-row-major order.
// Args: input degree m a b p hx hy hu hv pivot_rule basis.
int main(int argc,char** argv) {
  assert(argc==13);Local local{};std::ifstream input(argv[1]);assert(input);
  for(auto& row:local)for(auto& value:row){uint64_t raw;input>>raw;assert(input);value=raw;}
  int degree=std::stoi(argv[2]),m=std::stoi(argv[3]);
  std::array<int,3> lattice{std::stoi(argv[4]),std::stoi(argv[5]),std::stoi(argv[6])};
  std::array<int,4> physical{std::stoi(argv[7]),std::stoi(argv[8]),std::stoi(argv[9]),std::stoi(argv[10])};
  int pivot_rule=std::stoi(argv[11]),basis=std::stoi(argv[12]);
  auto factors=local_walks(auxiliary_conjugate(local,basis),degree,lattice,m,physical,pivot_rule);
  std::cout<<"[";
  for(int index=0;index<8;index++) {
    if(index)std::cout<<",";
    WalkChoice c{factors[index],lattice,physical,index,basis,pivot_rule,true,0,1};
    write_choice(std::cout,c);
  }
  std::cout<<"]\n";
}
