#include <vector>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <math.h>
#include <iostream>
#include <sstream>
#include <string>

#include "fonctions.hpp"
#include "matrix.hpp"
#include "charge.hpp"
#include "solveur.hpp"
#include "update.hpp"

void Update(std::vector<double> U,std::vector<double> V,int Nx,int Ny,int Nu,int Nv,double dt,double Lx,double Ly,double D,int mode,int h_part,double alpha,double beta,int Nt,double e,int kmax,double errschwz,int maxschwz){

  //construction des stencils
  std::vector<double> U0(3*Ny); //domaine 1 pour le domaine 2
  std::vector<double> V0(3*Ny); //domaine 2 pour le domaine 1

  //construction des matrices de resolution sur chaque sous domaine
  std::vector<int> rowu, rowv;
  std::vector<int> colu, colv;
  std::vector<double> valu, valv;

  double error(1.); int iteschwz(0);
  double t;
  
  //construction matrices de resolution
  //domaine 1
  Matrix(rowu,colu,valu,Nx,Ny,Nu,Lx,Ly,D,dt,alpha,beta,0);
  //domaine 2
  Matrix(rowv,colv,valv,Nx,Ny,Nv,Lx,Ly,D,dt,alpha,beta,1);

  /*
  for (int k=0;k<rowu.size();k++){
    if(rowu[k]%Nu==Nu-1) printf("col[%d], row[%d]  val[%f]\n",colu[k],rowu[k],valu[k]);
  }
  printf("\n");
  for (int k=0;k<rowv.size();k++){
    if(rowv[k]%Nv==0) printf("col[%d], row[%d]  val[%f]\n",colv[k],rowv[k],valv[k]);
  }
  printf("\n");
  */
  
  //schema en temps
  for (int k = 0; k < Nt; k++){

    t = k*dt; //temps de l'experience

    //Initialisation des stencils
    for (int i = 0; i < U0.size(); i++) U0[i]=1.,V0[i]=1.;

    //Boucle de Schwarz
    while ((iteschwz <= maxschwz) && (error >= errschwz)){
      
      //construction du second membre sur chacun des sous-domaines
      std::vector<double> Su(Nu*Ny,0.);
      std::vector<double> Sv(Nv*Ny,0.);

      secondMembre(Su,U,U0,Nx,Ny,Nu,dt,t,Lx,Ly,D,mode,alpha,beta,0);
      BICGStab(rowu,colu,valu,U,Su,e,kmax,Nu,Ny);

      secondMembre(Sv,V,V0,Nx,Ny,Nv,dt,t,Lx,Ly,D,mode,alpha,beta,1);
      BICGStab(rowv,colv,valv,V,Sv,e,kmax,Nv,Ny);

      //mise a jour des stencils
      for (int j = 0; j < Ny; j++){
	V0[j] = U[Nu*j+Nu-h_part-1];
	V0[j+Ny] = U[Nu*j+Nu-h_part];
	V0[j+2*Ny] = U[Nu*j+Nu-h_part+1];
	U0[j] = V[Nv*j+h_part-2];
	U0[j+Ny] = V[Nv*j+h_part-1];
	U0[j+2*Ny] = V[Nv*j+h_part];
      }
      
      // Evaluation de l'erreur pour la condition d'arret error
      error = maj_error(U,V,h_part,Ny,Nu,Nv);
      //printf("ite %d   erreur de schwarz=%f \n",iteschwz,error);
      iteschwz++;
    }
    if (iteschwz>=maxschwz) printf("Schwarz n'a pas convergé ! Erreur Schwarz : %f\n",error);

    //ecriture dans les  fichiers
    Write(U,V,Nx,Ny,Nu,Nv,Lx,Ly,h_part,k);
    
  }
}

void Write(std::vector<double> U,std::vector<double> V,int Nx,int Ny,int Nu,int Nv,double Lx,double Ly,int h_part,int k){
  
  std::string prefixeu="solutionU";
  std::ostringstream oss;
  oss<<prefixeu<<k<<".txt";
  std::ofstream NewFichieru(oss.str().c_str(), std::ios::out | std::ios::trunc);

  std::string prefixev="solutionV";
  std::ostringstream oss2;
  oss2<<prefixev<<k<<".txt";
  std::ofstream NewFichierv(oss2.str().c_str(), std::ios::out | std::ios::trunc);

  std::vector<double> x_tab(Nx);
  std::vector<double> y_tab(Ny);

  for (int i = 0; i < Nx; i++) x_tab[i]=i*Lx/(Nx+1);
  for (int j = 0; j < Ny; j++) y_tab[j]=j*Ly/(Ny+1);

  for (int i = 0; i < Nu; i++){
    for (int j = 0; j < Ny; j++) NewFichieru<< x_tab[i] <<" "<<  y_tab[j] << " " << U[i+j*Nu] << std::endl;
  }

  for (int i = 0; i < Nv; i++){
    for (int j = 0; j < Ny; j++) NewFichierv<< x_tab[i+Nu-h_part] <<" "<<  y_tab[j] << " " << V[i+j*Nv] << std::endl;
    }
}
