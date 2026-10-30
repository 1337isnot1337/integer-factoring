#pragma once
#include <array>
#include <bit>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
using U=uint64_t;
struct T {U a=0,b=0;};
T neg(T x){return {x.b,x.a};}
T add(T x,T y){return {(x.a&~(y.a|y.b))|(y.a&~(x.a|x.b))|(x.b&y.b),(x.b&~(y.a|y.b))|(y.b&~(x.a|x.b))|(x.a&y.a)};}
T scale(T x,int c){return c==1?x:c==2?neg(x):T{};}
struct Row {std::vector<U>a,b;};
using Matrix=std::vector<Row>;
int rank3(Matrix M,int n){int rank=0,w=(n+63)/64;for(int c=0;c<n;c++){int word=c/64,bit=c%64,p=rank;while(p<n&&!(((M[p].a[word]|M[p].b[word])>>bit)&1))p++;if(p==n)continue;std::swap(M[p],M[rank]);if((M[rank].b[word]>>bit)&1)std::swap(M[rank].a,M[rank].b);for(int i=rank+1;i<n;i++){int co=((M[i].a[word]>>bit)&1)+2*((M[i].b[word]>>bit)&1);if(!co)continue;for(int h=word;h<w;h++){T x{M[i].a[h],M[i].b[h]},y{M[rank].a[h],M[rank].b[h]};if(co==1)y=neg(y);auto z=add(x,y);M[i].a[h]=z.a;M[i].b[h]=z.b;}}rank++;if(rank==n)break;}return rank;}
int mod(int x){x%=3;return x<0?x+3:x;}
#include <algorithm>
#include <cassert>
#include <chrono>
int J=16;U PMASK=(U(1)<<16)-1;
T truncp(T x){return{x.a&PMASK,x.b&PMASK};}
T pmul(T x,T y){T out{};U bits=x.a|x.b;while(bits){int k=std::countr_zero(bits);bits&=bits-1;T z{(y.a<<k)&PMASK,(y.b<<k)&PMASK};if((x.b>>k)&1)z=neg(z);out=add(out,z);}return out;}
T invunit(T x){assert((x.a|x.b)&1);T r=(x.a&1)?T{1,0}:T{0,1};for(int k=1;k<J;k*=2){auto ar=pmul(x,r);r=pmul(r,add(T{0,1},neg(ar)));}auto z=pmul(x,r);assert(z.a==1&&z.b==0);return r;}
using PolyMatrix=std::vector<std::vector<T>>;
std::vector<int>smith(PolyMatrix A,int n){int m=A.size();std::vector<int>ans;for(int k=0;k<n;k++){int best=J,pi=-1,pj=-1;for(int i=k;i<m;i++){for(int j=k;j<n;j++){U bits=A[i][j].a|A[i][j].b;if(bits){int v=std::countr_zero(bits);if(v<best){best=v;pi=i;pj=j;}}}if(best==0)break;}if(pi<0)throw std::runtime_error("precision exhausted");std::swap(A[k],A[pi]);if(pj!=k)for(auto&row:A)std::swap(row[k],row[pj]);T unit{A[k][k].a>>best,A[k][k].b>>best};T inverse=invunit(unit);std::vector<int>nz;for(int j=k;j<n;j++){A[k][j]=pmul(A[k][j],inverse);if(j>k&&(A[k][j].a|A[k][j].b))nz.push_back(j);}assert(A[k][k].a==(U(1)<<best)&&!A[k][k].b);for(int i=k+1;i<m;i++){T v=A[i][k];if(!(v.a|v.b))continue;T factor{v.a>>best,v.b>>best};A[i][k]={};for(int j:nz)A[i][j]=add(A[i][j],neg(pmul(factor,A[k][j])));}ans.push_back(best);}return ans;}
using PM=std::array<T,4>;
PM pmatmul(PM A,PM B){PM C{};for(int i=0;i<2;i++)for(int j=0;j<2;j++)for(int k=0;k<2;k++)C[2*i+j]=add(C[2*i+j],pmul(A[2*i+k],B[2*k+j]));return C;}
T coeffpoly(int c,int power){return c==1?T{U(1)<<power,0}:c==2?T{0,U(1)<<power}:T{};}
T reversep(T x,int L){T r{};for(int j=0;j<=L;j++){if(x.a>>j&1)r.a|=U(1)<<(L-j);if(x.b>>j&1)r.b|=U(1)<<(L-j);}return r;}
std::vector<PM> blocked(int q,int b,int&d){d=1;for(int i=0;i<b;i++)d*=q;int kappa=mod(1+2*q);if(!kappa)throw std::runtime_error("bad q");int alpha=1,beta=mod(2*kappa);std::vector<PM>local(q*q*q*q);for(int x=0;x<q;x++)for(int y=0;y<q;y++)for(int u=0;u<q;u++)for(int v=0;v<q;v++){int R=x==y&&u==v,C=x==u&&y==v,E=x==v&&y==u;local[((x*q+y)*q+u)*q+v]={coeffpoly(R,0),coeffpoly(mod(beta*(R+2*E)),1),coeffpoly(mod(C-2*kappa*E),0),coeffpoly(mod(beta*C),1)};}std::vector<PM>all(d*d*d*d);for(int x=0;x<d;x++)for(int y=0;y<d;y++)for(int u=0;u<d;u++)for(int v=0;v<d;v++){int xx=x,yy=y,uu=u,vv=v;PM M{T{1,0},T{},T{},T{1,0}};for(int j=0;j<b;j++){M=pmatmul(M,local[((xx%q*q+yy%q)*q+uu%q)*q+vv%q]);xx/=q;yy/=q;uu/=q;vv/=q;}all[((x*d+y)*d+u)*d+v]=M;}return all;}
PolyMatrix flatten(std::vector<PM>const&K,int d,int cut,int side,int degree,bool reverse){int n=d*d;PolyMatrix A(4*n,std::vector<T>(n));for(int x=0;x<d;x++)for(int y=0;y<d;y++)for(int u=0;u<d;u++)for(int v=0;v<d;v++){int r=cut?x*d+u:x*d+y,c=cut?y*d+v:u*d+v;if(side)std::swap(r,c);auto M=K[((x*d+y)*d+u)*d+v];for(int a=0;a<2;a++)for(int b=0;b<2;b++)A[4*r+2*a+b][c]=reverse?reversep(M[2*a+b],degree):M[2*a+b];}return A;}
double objective(std::vector<int>const&d,std::vector<int>const&nu,double t){double A=0,B=0;for(int x:d)A+=std::exp(x*t);for(int x:nu)B+=std::exp(-x*t);return A*B;}
