#pragma once
#include "f3_mpo_core.hpp"
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
