#include <algorithm>
#include <bit>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <vector>

using U = std::uint64_t;
using Matrix = std::vector<U>;

int rank(Matrix rows, int dimension) {
  Matrix pivots(dimension);
  int answer = 0;
  for (U row : rows) {
    while (row) {
      int k = 63 - std::countl_zero(row);
      if (!pivots[k]) { pivots[k] = row; ++answer; break; }
      row ^= pivots[k];
    }
  }
  return answer;
}

U floor_sqrt(U value) {
  U r = U(std::sqrt((long double)value));
  while ((r+1)*(r+1) <= value) ++r;
  while (r*r > value) --r;
  return r;
}

int main(int argc, char** argv) {
  assert(argc == 2);
  std::ifstream in(argv[1]);
  int q, m, target;
  in >> q >> m >> target;
  int n = q*q, count = 1 << m;
  assert(in && n <= 64 && m <= 20 && target > 0 && target < count);
  std::vector<Matrix> row_basis(m, Matrix(n)), column_basis(m, Matrix(n));
  for (int i=0; i<m; ++i) {
    for (U &row : row_basis[i]) in >> row;
    for (int x=0; x<q; ++x) for (int y=0; y<q; ++y)
      for (int u=0; u<q; ++u) for (int v=0; v<q; ++v)
        if ((row_basis[i][x*q+y] >> (u*q+v)) & 1)
          column_basis[i][x*q+u] |= U(1) << (y*q+v);
  }
  assert(in);
  Matrix endpoint(n);
  for (int i=0; i<m; ++i) if ((target >> i)&1)
    for (int a=0; a<n; ++a) endpoint[a] ^= row_basis[i][a];
  for (int x=0; x<q; ++x) for (int y=0; y<q; ++y)
    assert(endpoint[x*q+y] == (U(1) << (y*q+x)));

  constexpr U scale = 1000000;
  std::vector<int> rr(count), cc(count);
  std::vector<U> costs(count);
  Matrix row(n), column(n);
  int old = 0;
  for (int index=1; index<count; ++index) {
    int mask = index ^ (index >> 1);
    int changed = std::countr_zero(unsigned(mask ^ old));
    for (int a=0; a<n; ++a) {
      row[a] ^= row_basis[changed][a];
      column[a] ^= column_basis[changed][a];
    }
    rr[mask] = rank(row,n); cc[mask] = rank(column,n);
    assert(rr[mask] && cc[mask]);
    costs[mask] = floor_sqrt(U(rr[mask])*cc[mask]*scale*scale);
    old = mask;
  }
  long double best_two = std::sqrt((long double)rr[target]*cc[target]);
  int split = 0;
  for (int mask=1; mask<count; ++mask) {
    long double value = std::sqrt((long double)rr[mask]*cc[mask])
                        + std::sqrt((long double)rr[mask^target]*cc[mask^target]);
    if (value < best_two) { best_two = value; split = mask; }
  }

  // Complete Dijkstra with integer edge costs rounded down. Its distance
  // is a rigorous lower bound for every real radical-cost decomposition.
  constexpr U infinity = std::numeric_limits<U>::max()/4;
  std::vector<U> distance(count,infinity);
  std::vector<int> previous(count,-1);
  std::vector<bool> used(count,false);
  distance[0] = 0;
  int popped = 0;
  for (int step=0; step<count; ++step) {
    U best = infinity; int u=-1;
    for (int v=0; v<count; ++v)
      if (!used[v] && distance[v]<best) { best=distance[v]; u=v; }
    assert(u>=0);
    used[u]=true; ++popped;
    if (u==target) break;
    for (int v=0; v<count; ++v) if (!used[v]) {
      U candidate=best+costs[u^v];
      if (candidate<distance[v]) { distance[v]=candidate; previous[v]=u; }
    }
  }
  std::vector<int> pieces;
  for (int u=target; u; u=previous[u]) {
    assert(previous[u]>=0);
    pieces.push_back(u^previous[u]);
  }
  long double upper=0;
  U upper_integer=0;
  for (int piece:pieces) {
    U square=U(rr[piece])*cc[piece]*scale*scale;
    U lower=floor_sqrt(square);
    upper_integer += lower + (lower*lower != square);
    upper += std::sqrt((long double)rr[piece]*cc[piece]);
  }
  std::cout<<std::setprecision(18)
    <<"{\"q\":"<<q<<",\"span_dimension\":"<<m
    <<",\"combinations\":"<<count<<",\"target\":"<<target
    <<",\"integer_scale\":"<<scale
    <<",\"minimum_cost_lower_integer\":"<<distance[target]
    <<",\"selected_cost_upper_integer\":"<<upper_integer
    <<",\"selected_cost\":"<<upper
    <<",\"best_two_piece_cost\":"<<best_two
    <<",\"best_two_piece_split\":"<<split
    <<",\"dijkstra_popped\":"<<popped<<",\"pieces\":[";
  for (std::size_t i=0;i<pieces.size();++i) {
    int mask=pieces[i];
    std::cout<<(i?",":"")<<"{\"mask\":"<<mask
      <<",\"row_rank\":"<<rr[mask]<<",\"column_rank\":"<<cc[mask]<<"}";
  }
  std::cout<<"]}\n";
}
