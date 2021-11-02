#include <vector>
#include<cstdio>
#include<cstdlib>
#include"math.h"
#include "fonctions.hpp"
#include "matrix.hpp"

//-----------------------------------------------------------------------------------------------------------------------------------------------
//fonction qui remplit la matrice de diffsuion du problème

void sparseMatrix(std::vector<int> &row,std::vector<int> &col,std::vector<double> &value, int Nx, int Ny, double Lx, double Ly, double D, double dt){
  // Fonction qui remplt la matrice associée au problème avec un stockage coordonnées
  double di(0),sdi(0),ssdi(0);
  double dx=Lx/(Nx+1),dy=Ly/(Ny+1);
  di=1+2*dt*(1/pow(dx,2)+1/pow(dy,2));
  sdi=-dt*D/pow(dx,2);
  //ssdi=-dt*D/pow(dy,2);
  for(int i=0; i<(Nx*Ny); i++){
    col.push_back(i);
    row.push_back(i);
    value.push_back(di);
    if((i+1)%Nx==0){}
    else{
      col.push_back(i); //stockage de la sous-diagonale
      row.push_back(i+1);
      value.push_back(sdi);
      col.push_back(i+1); //stockage de la sur-diagonale (symetrie)
      row.push_back(i);
      value.push_back(sdi);
    }
    /*if(i+Nx<Nx*Nx){
    col.push_back(i); //stockage de la sous-sous-diagonale
    row.push_back(i+Nx);
    value.push_back(ssdi);
    col.push_back(i+Nx); //stockage de la sur-sur-diagonale (symetrie)
    row.push_back(i);
    value.push_back(ssdi);
  }*/
  }
}

//-----------------------------------------------------------------------------------------------------------------------
//fonction qui remplit la matrice de diffusion partielle

void sparseMatrix_parallel(std::vector<int> &row,std::vector<int> &col,std::vector<double> &value, int Nx, int Ny, double Lx, double Ly, double D, double dt,int iBeg,int iEnd){
  // Fonction qui remplt la matrice associée au problème avec un stockage coordonnées
  double di(0),sdi(0),ssdi(0),dx=Lx/(Nx+1),dy=Ly/(Ny+1);
  di=1+2*dt*(1/pow(dx,2)+1/pow(dy,2));
  sdi=-dt*D/pow(dx,2);
  ssdi=-dt*D/pow(dy,2);

  for(int i=iBeg; i<iEnd+1; i++){
    col.push_back(i);
    row.push_back(i);
    value.push_back(di);}

  for (int i = 0; i < Nx*Ny ; i++) {
    if((i+1)%Nx==0){}
    else{
      if (i+1<=iEnd && i+1>=iBeg) {
        col.push_back(i); //stockage de la sous-diagonale
        row.push_back(i+1);
        value.push_back(sdi);}

      if(i<=iEnd && i>=iBeg)
      {col.push_back(i+1); //stockage de la sur-diagonale (symetrie)
      row.push_back(i);
      value.push_back(sdi);
    }
  }

      if(i+Nx<Nx*Nx && i+Nx<=iEnd && i+Nx>=iBeg){
      col.push_back(i); //stockage de la sous-sous-diagonale
      row.push_back(i+Nx);
      value.push_back(ssdi);}

      if(i+Nx<Nx*Nx && i<=iEnd && i>=iBeg){
      col.push_back(i+Nx); //stockage de la sur-sur-diagonale (symetrie)
      row.push_back(i);
      value.push_back(ssdi);
      }
  }
}

// ---------------------------------------------------------
//fonction qui remplit le vecteur second mambre du problème

void secondMembre(std::vector<double> &S,std::vector<double> U, int Nx, int Ny,double dt, double t, double Lx, double Ly, double D, int mode){
  // Fonction qui remplit le second membre
  double dx=Lx/(Nx+1),dy=Ly/(Ny+1);
  for(int i=0; i<Nx; i++){
    for(int j=0; j<Ny; j++){
      S[i+j*Nx]=U[i+j*Nx]+dt*f((i+1)*dx,(j+1)*dy,t,Lx,Ly,mode);
      if (j==0) {
        S[i+j*Nx]+=D*dt/pow(dy,2)*g((i+1)*dx,0,mode);
      }
      else if (j==Ny-1) {
        S[i+j*Nx]+=D*dt/pow(dy,2)*g((i+1)*dx,Ny+1,mode);
      }
      else if (i==0) {
        S[i+j*Nx]+=D*dt/pow(dx,2)*h(0,(j+1)*dy,mode);
      }
      else if (i==Nx-1) {
        S[i+j*Nx]+=D*dt/pow(dx,2)*h(Nx+1,(j+1)*dy,mode);
      }
    }
  }
}

//-----------------------------------------------------------------------------------
//fonction qui remplit le vecteur second membre partiel

void secondMembre_parallel(std::vector<double> &S,std::vector<double> U, int Nx, int Ny,double dt, double t, double Lx, double Ly, double D, int mode,int iBeg,int iEnd){
  double dx=Lx/(Nx+1),dy=Ly/(Ny+1);

  for (int k = 0; k < iEnd-iBeg+1; k++) {
    int i=(k+iBeg)%Nx;
    int j=(k+iBeg)/Ny;
    S[k]=U[k]+dt*f((i+1)*dx,(j+1)*dy,t,Lx,Ly,mode);
    if (j==0) {
      S[k]+=D*dt/pow(dy,2)*g((i+1)*dx,0,mode);
    }
    else if (j==Ny-1) {
      S[k]+=D*dt/pow(dy,2)*g((i+1)*dx,Ny+1,mode);
    }
    else if (i==0) {
      S[k]+=D*dt/pow(dx,2)*h(0,(j+1)*dy,mode);
    }
    else if (i==Nx-1) {
      S[k]+=D*dt/pow(dx,2)*h(Nx+1,(j+1)*dy,mode);
    }
  }
}
