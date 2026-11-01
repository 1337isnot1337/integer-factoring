// Standalone exact F3 polynomial/Smith probe; standard C++20 only.
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
using GM=std::vector<T>;
GM gmul(GM const&A,GM const&B,int D){GM C(D*D);for(int i=0;i<D;i++)for(int j=0;j<D;j++)for(int k=0;k<D;k++)C[D*i+j]=add(C[D*i+j],pmul(A[D*i+k],B[D*k+j]));return C;}
std::vector<GM> resolvent_block(int q,int b,int&d){d=1;for(int j=0;j<b;j++)d*=q;int D=q;std::vector<GM>local(q*q*q*q,GM(D*D));T PA{7,0},PB{4,2},PE{4,0};for(int u=0;u<q;u++)for(int v=0;v<q;v++)for(int c=0;c<D;c++){int x=c,y=v,a=u;auto&i=local[((x*q+y)*q+u)*q+v][D*a+c];i=add(i,PA);x=u;y=c;a=v;auto&j=local[((x*q+y)*q+u)*q+v][D*a+c];j=add(j,PB);x=v;y=u;a=c;auto&k=local[((x*q+y)*q+u)*q+v][D*a+c];k=add(k,PE);}
for(int x=0;x<q;x++)for(int y=0;y<q;y++)for(int u=0;u<q;u++)for(int v=0;v<q;v++)for(int a=0;a<D;a++)for(int c=0;c<D;c++){T f=local[((x*q+y)*q+u)*q+v][D*a+c];if((std::popcount(f.a)+2*std::popcount(f.b))%3!=(x==v&&y==u&&a==c))throw std::runtime_error("endpoint");}
std::vector<GM>all(d*d*d*d,GM(D*D));for(int x=0;x<d;x++)for(int y=0;y<d;y++)for(int u=0;u<d;u++)for(int v=0;v<d;v++){int xx=x,yy=y,uu=u,vv=v;GM M(D*D);for(int a=0;a<D;a++)M[D*a+a]={1,0};for(int j=0;j<b;j++){M=gmul(M,local[((xx%q*q+yy%q)*q+uu%q)*q+vv%q],D);xx/=q;yy/=q;uu/=q;vv/=q;}all[((x*d+y)*d+u)*d+v]=std::move(M);}return all;}
PolyMatrix gflatten(std::vector<GM>const&K,int d,int D,int cut,int side,int degree,bool reverse){int n=d*d;PolyMatrix A(D*D*n,std::vector<T>(n));for(int x=0;x<d;x++)for(int y=0;y<d;y++)for(int u=0;u<d;u++)for(int v=0;v<d;v++){int r=cut?x*d+u:x*d+y,c=cut?y*d+v:u*d+v;if(side)std::swap(r,c);auto&M=K[((x*d+y)*d+u)*d+v];for(int a=0;a<D;a++)for(int b=0;b<D;b++)A[D*D*r+D*a+b][c]=reverse?reversep(M[D*a+b],degree):M[D*a+b];}return A;}
#include <map>
#include <random>
using Small=std::vector<unsigned char>;
using V9=std::array<int,9>;
bool nilpotent_certificate(V9 const&U,V9 const&V){int B[3][3]{},C[3][3]{};for(int i=0;i<3;i++)for(int j=0;j<3;j++){for(int k=0;k<3;k++){B[i][j]+=U[3*i+k]*V[3*k+j];C[i][j]+=V[3*j+k]*U[3*k+i];}B[i][j]%=3;C[i][j]%=3;}if((B[0][0]+B[1][1]+B[2][2])%3)return false;int e2=0;for(int i=0;i<3;i++)for(int j=i+1;j<3;j++)e2+=B[i][i]*B[j][j]-B[i][j]*B[j][i];if(mod(e2))return false;
int M[3][3][4]{};for(int i=0;i<3;i++)for(int j=0;j<3;j++){M[i][j][0]=mod(B[i][j]+C[i][j]);M[i][j][1]=mod(B[i][j]-C[i][j]);M[i][j][2]=B[i][j];M[i][j][3]=i==j;}int det[10]{};int perms[6][4]={{0,1,2,1},{0,2,1,-1},{1,0,2,-1},{1,2,0,1},{2,0,1,1},{2,1,0,-1}};for(auto&p:perms)for(int i=0;i<4;i++)for(int j=0;j<4;j++)for(int k=0;k<4;k++)det[i+j+k]+=p[3]*M[0][p[0]][i]*M[1][p[1]][j]*M[2][p[2]][k];for(int k=0;k<9;k++)if(mod(det[k]))return false;return mod(det[9])==1;}
std::vector<GM> candidate(V9 const&U,V9 const&V,int&degree){constexpr int q=3,D=3,n=27;Small B(n*n),N(n*n);for(int x=0;x<q;x++)for(int y=0;y<q;y++)for(int a=0;a<D;a++)for(int u=0;u<q;u++)for(int v=0;v<q;v++)for(int b=0;b<D;b++){int r=(x*q+y)*D+a,c=(u*q+v)*D+b;if(y==v)B[r*n+c]=mod(int(x==b&&a==u)+U[3*x+a]*V[3*b+u]);}
for(int r=0;r<n;r++)for(int u=0;u<q;u++)for(int v=0;v<q;v++)for(int b=0;b<D;b++){int c=(u*q+v)*D+b,ce=(v*q+u)*D+b;N[r*n+c]=mod(int(r==c)-B[r*n+ce]);}
std::vector<GM>K(q*q*q*q,GM(D*D));Small C=B;degree=0;for(int k=0;k<=27;k++){bool nonzero=false;for(auto x:C)nonzero|=x!=0;if(!nonzero){degree=k-1;break;}if(k==27)throw std::runtime_error("nilpotence certificate failed");for(int x=0;x<q;x++)for(int y=0;y<q;y++)for(int u=0;u<q;u++)for(int v=0;v<q;v++)for(int a=0;a<D;a++)for(int b=0;b<D;b++){int r=(x*q+y)*D+a,c=(u*q+v)*D+b,co=C[r*n+c];auto&t=K[((x*q+y)*q+u)*q+v][D*a+b];if(co==1)t.a|=uint64_t(1)<<k;if(co==2)t.b|=uint64_t(1)<<k;}Small next(n*n);for(int i=0;i<n;i++)for(int j=0;j<n;j++){int value=0;for(int h=0;h<n;h++)value+=N[i*n+h]*C[h*n+j];next[i*n+j]=value%3;}C=std::move(next);}
for(int x=0;x<q;x++)for(int y=0;y<q;y++)for(int u=0;u<q;u++)for(int v=0;v<q;v++)for(int a=0;a<D;a++)for(int b=0;b<D;b++){T f=K[((x*q+y)*q+u)*q+v][D*a+b];if((std::popcount(f.a)+2*std::popcount(f.b))%3!=(x==v&&y==u&&a==b))throw std::runtime_error("endpoint");}return K;}

// Reproduce two fixed exact probes. All field work above is integer mod 3.
int main(int argc,char**argv){
 std::string path=argc>1?argv[1]:"";J=64;PMASK=~uint64_t(0);
 auto started=std::chrono::steady_clock::now();
 std::array<V9,2> Us{V9{},V9{1,2,0,0,2,0,1,2,1}};
 std::array<V9,2> Vs{V9{},V9{1,1,0,2,2,0,2,0,0}};
 std::ofstream file;if(!path.empty())file.open(path);
 auto output=[&](std::ostream&o){
  o<<"{\"field\":3,\"q\":3,\"auxiliary\":3,\"precision\":64,\"records\":[";
  for(int model=0;model<2;model++){
   if(model)o<<",";int L;assert(nilpotent_certificate(Us[model],Vs[model]));
   auto K=candidate(Us[model],Vs[model],L);
   o<<"{\"model\":"<<model<<",\"degree\":"<<L<<",\"U\":[";
   for(int j=0;j<9;j++)o<<(j?",":"")<<Us[model][j];o<<"],\"V\":[";
   for(int j=0;j<9;j++)o<<(j?",":"")<<Vs[model][j];o<<"],\"profiles\":[";
   for(int cut=0;cut<2;cut++)for(int side=0;side<2;side++){
    int i=2*cut+side;if(i)o<<",";
    auto ns=smith(gflatten(K,3,3,cut,side,L,false),9);
    auto ds=smith(gflatten(K,3,3,cut,side,L,true),9);for(auto&v:ds)v=L-v;
    o<<"{\"cut_side\":"<<i<<",\"d\":[";
    for(size_t j=0;j<ds.size();j++)o<<(j?",":"")<<ds[j];o<<"],\"nu\":[";
    for(size_t j=0;j<ns.size();j++)o<<(j?",":"")<<ns[j];o<<"]}";
   }o<<"]}";
  }
  o<<"],\"seconds\":"<<std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count()<<"}\n";
 };
 if(file)output(file);else output(std::cout);return 0;
}
