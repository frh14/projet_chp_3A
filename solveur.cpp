#include<vector>
#include<cstdio>
#include<cstdlib>
#include<math.h>
//#include<mpi.h>
#include"solveur.hpp"
#include"charge.hpp"

//--------------------------------------------------------------------------------------------------------------------------------------------------------------
//fonction qui résoud par la méthode du gradient conjugué un systéme linéaire Ax=b
//entrée: A matrice du systéme sous stockage indice ligne/colonne + valeur
//second membre b
//point de départ du GC x0

void GC(std::vector<int> row,std::vector<int> col,std::vector<double> value,std::vector<double> &U,std::vector<double> F,double e,int kmax,int Nx,int Ny)
{

//stockage du vecteur inconnu dans un vecteur intermédiaire x
  std::vector<double> x(Nx*Ny,0);
  for (int i = 0; i < Nx*Ny; i++) {
    x[i]=U[i];}

//calcul de la première direction de descente p0
  std::vector<double> pprime(Nx*Ny,0);
  mulSparseMatrix(row,col,value,pprime,x);

  std::vector<double> p(Nx*Ny,0);
  for (int i = 0; i < Nx*Ny; i++) {
    p[i]=F[i]-pprime[i];}

//calcul du reste initiale r0
  std::vector<double> r(Nx*Ny,0);
  for (int i = 0; i < Nx*Ny; i++) {
    r[i]=p[i];}

//initialisation
  double beta=norme2(r);

  int k=0;

//début de la méthode itérative
  while (beta>e && k<=kmax) {

//calcul du vecteur  z=Ap
    std::vector<double> z(Nx*Ny,0);
    mulSparseMatrix(row,col,value,z,p);

//calcul du pas de descente alpha
    double w=ps(r,r);
    double alpha=w/ps(z,p);

//mise à jour du vecteur inconnu x
    for (int i = 0; i < Nx*Ny; i++) {
      x[i]+=alpha*p[i];}

//calcul du nouveau reste
    std::vector<double> rplus(Nx*Ny,0);
    for (int i = 0; i < Nx*Ny; i++) {
      rplus[i]=r[i]-alpha*z[i];}

//mise à jour de la nouvelle direction de descente
    double gamma=ps(rplus,rplus)/w;

    for (int i=0; i < Nx*Ny; i++) {
      p[i]=rplus[i]+gamma*p[i];}

//mise à jour du reste
    for (int i = 0; i < Nx*Ny; i++) {
      r[i]=rplus[i];}

//mise à jour de beta norme du reste
    beta=norme2(r);

//incrémentation
    k+=1;}

//fin des  itérations

//mise à jour du vecteur inconnu avec le résultat du GC
  for (int i = 0; i < Nx*Ny; i++) {
    U[i]=x[i];}
  }

//-------------------------------------------------------------------------------------------------------------------------------------------------------------------
//fonction qui prend en entrée un vecteur x
//et qui calcule la norme quadratique de ce vecteur

double norme2(std::vector<double> x){
  double norme=0.;
  for (int i = 0; i < x.size(); i++) {
    int n=x.size();
    norme+=pow(x[i],2);
  }
  norme=sqrt(norme);
  return norme;
}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------------
//fonction qui prend en entrée deux vecteurs
// et qui en calcule le produit scalaire canonique

double ps(std::vector<double> x,std::vector<double> y){
  double ps=0.;
  for (int i = 0; i < x.size(); i++) {
    ps+=x[i]*y[i];}
  return ps;
}

//--------------------------------------------------------------------------------------------------------------------------------------------------------------------
//produit matrice/vecteur classique pour un stockage creux
//fonction qui prend en entrée une matrice A sous la forme:
//vecteur valeur->contient la valeur des termes non nuls de la matrice
//vecteur ligne et colonne->contiennent les indices des termes non nuls
//le vecteur x
//et calcul le vecteur y=AX

void mulSparseMatrix(std::vector<int> row, std::vector<int> col, std::vector<double> value, std::vector<double> &y, std::vector<double> x) {
  for (int k = 0; k < row.size(); k++){
    y[row[k]]+=value[k]*x[col[k]];}
}

//-----------------------------------------------------------------------------------------------------------------
//produit matrice/vecteur pour une application en paralléle
//fonction qui prend en entrée une matrice A sous la forme:
//vecteur valeur->contient la valeur des termes non nuls de la matrice
//vecteur ligne et colonne->contiennent les indices des termes non nuls
//le vecteur x sous la forme
//x_val->valeur des composantes de x
//x_indice -> indices des composantes de x
//et calcul le vecteur y=AX

void mulSparseMatrix_2(std::vector<int> row, std::vector<int> col, std::vector<double> value, std::vector<double> &y,int iBeg,int iEnd,std::vector<double> x_val,std::vector<int> x_indice){
  for (int k = 0; k < row.size(); k++) {
    for (int i = 0; i < x_indice.size(); i++) {
      if (x_indice[i]==col[k] && row[k]<=iEnd && row[k]>=iBeg) {
        y[row[k]-iBeg]+=value[k]*x_val[i];}}}
  }
