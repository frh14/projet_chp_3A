#include<vector>
#include<cstdio>
#include<cstdlib>
#include<math.h>
#include"solveur.hpp"
#include"charge.hpp"

  //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  //--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  //fonction qui résoud par la méthode du gradient bi-conjugué stabilisé un systéme linéaire Ax=b
  //entrée: A matrice du systéme sous stockage indice ligne/colonne + valeur
  //second membre b
  //point de départ du BCG x0

  void BICG(std::vector<int> row,std::vector<int> col,std::vector<double> val,std::vector<double> &X,std::vector<double> &F,double e,int kmax,int Nx,int Ny){

/* //reconditionnement de la matrice
    for (int i = 0; i < row.size(); i++) {
      if (row[i]==col[i]) {
        F[row[i]]=F[row[i]]/val[i],val[i]=1.0;}}*/

    //initialisation du reste
    std::vector<double> r(Nx*Ny,0),mu(Nx*Ny,0),p(Nx*Ny,0),r0(Nx*Ny,0);
    mulSparseMatrix(row,col,val,r,X);

    for (int i = 0; i < Nx*Ny; i++) r[i]=F[i]-r[i];

    double rho(1.0),alpha(1.0),omega(1.0);

    r0=r;

    int k=0;
    double norm=norm2(r);

    while (k<=kmax && norm>=e) {

      double rhoplus=ps(r0,r);

      double beta=rhoplus/rho*alpha/omega;

      for (int i = 0; i < Nx*Ny; i++) p[i]=r[i]+beta*(p[i]-omega*mu[i]);

      mulSparseMatrix(row,col,val,mu,p);

      alpha=rhoplus/ps(r0,mu);
      std::vector<double> h(Nx*Ny,0);
      std::vector<double> s(Nx*Ny,0);
      for (int i = 0; i < Nx*Ny; i++) {
        h[i]=X[i]+alpha*p[i];
        s[i]=r[i]-alpha*mu[i];}

      std::vector<double> t(Nx*Ny,0);
      mulSparseMatrix(row,col,val,t,s);

      omega=ps(t,s)/ps(t,t);

      for (int i = 0; i < Nx*Ny; i++) {
        X[i]=h[i]+omega*s[i];
        r[i]=s[i]-omega*t[i];}

      k=k+1;
      norm=norm2(r);
      rho=rhoplus;
    }

  }


//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//fonction qui prend en entrée un vecteur x
//et qui calcule la norme quadratique de ce vecteur

double norm2(std::vector<double> x){
  double norm=0.;
  for (int i = 0; i < x.size(); i++) norm+=x[i]*x[i];
  norm=sqrt(norm);
  return norm;
}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------------
//fonction qui prend en entrée deux vecteurs
// et qui en calcule le produit scalaire canonique

double ps(std::vector<double> x,std::vector<double> y){
  double ps=0.;
  for (int i = 0; i < x.size(); i++) ps+=x[i]*y[i];
  return ps;
}

//--------------------------------------------------------------------------------------------------------------------------------------------------------------------
//produit matrice/vecteur classique pour un stockage creux
//fonction qui prend en entrée une matrice A sous la forme:
//vecteur valeur->contient la valeur des termes non nuls de la matrice
//vecteur ligne et colonne->contiennent les indices des termes non nuls
//le vecteur x
//et calcul le vecteur y=AX

void mulSparseMatrix(std::vector<int> row, std::vector<int> col, std::vector<double> val, std::vector<double> &y, std::vector<double> x) {

  for (int i = 0; i < y.size(); i++) y[i]=0.;

  for (int k = 0; k < row.size(); k++) y[row[k]]+=val[k]*x[col[k]];

}
