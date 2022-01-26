#include <vector>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <math.h>
#include <iostream>
#include <sstream>
#include <string>
#include <mpi.h>

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

    printf("\n Begin time iteration %d\n",k);

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
    if (iteschwz>=maxschwz) printf("\n Schwarz n'a pas convergé ! Erreur Schwarz : %f\n",error);

    printf("\n End time iteration %d , time = %f\n",k,t);

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

void Update_(int argc, char** argv,int Nx,int Ny,double dt,double Lx,double Ly,double D,int mode,int h_part,double alpha,double beta,int Nt,double e,int kmax,double errschwz,int maxschwz){

  //construction des stencils
  std::vector<double> U0(3*Ny); //domaine 1 pour le domaine 2
  std::vector<double> V0(3*Ny); //domaine 2 pour le domaine 1

  //schema en temps
  for (int k = 0; k < Nt; k++){

    double t = k*dt; //temps de l'experience

    //Initialisation des stencils
    for (int i = 0; i < U0.size(); i++) U0[i]=1.,V0[i]=1.;

    MPI_Init(&argc,&argv);
    int me,Nproc;
    MPI_Comm_rank(MPI_COMM_WORLD,&me);  //processeur local
    MPI_Comm_size(MPI_COMM_WORLD,&Nproc); //nombre de processeur total utilisés

    printf("thread %d on %d \n",me,Nproc);

    int IBeg, IEnd;

    charge(Nx,Nproc,me,IBeg,IEnd);

    if (me==0) {}
    else{IBeg-=h_part;}

    int N=IEnd-IBeg+1;

    printf("thread %d N==%d \n",me,N);

    std::vector<double> U_loc(N*Ny,1);

    std::vector<double> U0_loc(3*Ny,0.);
    std::vector<double> V0_loc(3*Ny,0.);

    std::vector<int> row_loc, col_loc;
    std::vector<double> val_loc;

    Matrix_(row_loc,col_loc,val_loc,Nx,Ny,N,Nproc,Lx,Ly,D,dt,alpha,beta,me);

    MPI_Status Status;

    int iteschwz_loc=0.;
    double error_loc=1.;

    while ((iteschwz_loc <= maxschwz) && (error_loc >= errschwz)){

      std::vector<double> S_loc(N*Ny,0.);

      secondMembre_(S_loc,U_loc,U0_loc,V0_loc,Nx,Ny,N,Nproc,dt,t,Lx,Ly,D,mode,alpha,beta,me);

      BICGStab(row_loc,col_loc,val_loc,U_loc,S_loc,e,kmax,N,Ny);

      if (me==0) { // Proc 0
        for (int j = 0; j < Ny; j++){
          // Proc 0 envoit U a proc 1
          MPI_Send(&U_loc[N*j+N-h_part-1],1,MPI_DOUBLE,me+1,0,MPI_COMM_WORLD);
          MPI_Send(&U_loc[N*j+N-h_part],1,MPI_DOUBLE,me+1,1,MPI_COMM_WORLD);
          MPI_Send(&U_loc[N*j+N-h_part+1],1,MPI_DOUBLE,me+1,2,MPI_COMM_WORLD);

          // Proc 0 recoit stencil droite U0 de proc 1
          MPI_Recv(&U0_loc[j],1,MPI_DOUBLE,me+1,0,MPI_COMM_WORLD,&Status);
          MPI_Recv(&U0_loc[j+Ny],1,MPI_DOUBLE,me+1,1,MPI_COMM_WORLD,&Status);
          MPI_Recv(&U0_loc[j+2*Ny],1,MPI_DOUBLE,me+1,2,MPI_COMM_WORLD,&Status);
        }
      }

      else if (me==Nproc-1) { // Proc Nproc-1
        for (int j = 0; j < Ny; j++){
          // Proc Nproc-1 envoit U a proc Nproc-2
          MPI_Send(&U_loc[N*j+h_part-2],1,MPI_DOUBLE,me-1,0,MPI_COMM_WORLD);
          MPI_Send(&U_loc[N*j+h_part-1],1,MPI_DOUBLE,me-1,1,MPI_COMM_WORLD);
          MPI_Send(&U_loc[N*j+h_part],1,MPI_DOUBLE,me-1,2,MPI_COMM_WORLD);

          // Proc Nproc-1 recoit stencil gauche V0 de proc Nproc-2
          MPI_Recv(&V0_loc[j],1,MPI_DOUBLE,me-1,0,MPI_COMM_WORLD,&Status);
          MPI_Recv(&V0_loc[j+Ny],1,MPI_DOUBLE,me-1,1,MPI_COMM_WORLD,&Status);
          MPI_Recv(&V0_loc[j+2*Ny],1,MPI_DOUBLE,me-1,2,MPI_COMM_WORLD,&Status);
        }
      }

      else {
        for (int j = 0; j < Ny; j++){
          // Proc me envoit U a proc me-1
          MPI_Send(&U_loc[N*j+h_part-2],1,MPI_DOUBLE,me-1,0,MPI_COMM_WORLD);
          MPI_Send(&U_loc[N*j+h_part-1],1,MPI_DOUBLE,me-1,1,MPI_COMM_WORLD);
          MPI_Send(&U_loc[N*j+h_part],1,MPI_DOUBLE,me-1,2,MPI_COMM_WORLD);

          // Proc me recoit stencil gauche V0 de proc me-1
          MPI_Recv(&V0_loc[j],1,MPI_DOUBLE,me-1,0,MPI_COMM_WORLD,&Status);
          MPI_Recv(&V0_loc[j+Ny],1,MPI_DOUBLE,me-1,1,MPI_COMM_WORLD,&Status);
          MPI_Recv(&V0_loc[j+2*Ny],1,MPI_DOUBLE,me-1,2,MPI_COMM_WORLD,&Status);

          // Proc me envoit U a proc me+1
          MPI_Send(&U_loc[N*j+N-h_part-1],1,MPI_DOUBLE,me+1,0,MPI_COMM_WORLD);
          MPI_Send(&U_loc[N*j+N-h_part],1,MPI_DOUBLE,me+1,1,MPI_COMM_WORLD);
          MPI_Send(&U_loc[N*j+N-h_part+1],1,MPI_DOUBLE,me+1,2,MPI_COMM_WORLD);

          // Proc me recoit stencil droite U0 de proc me+1
          MPI_Recv(&U0_loc[j],1,MPI_DOUBLE,me+1,0,MPI_COMM_WORLD,&Status);
          MPI_Recv(&U0_loc[j+Ny],1,MPI_DOUBLE,me+1,1,MPI_COMM_WORLD,&Status);
          MPI_Recv(&U0_loc[j+2*Ny],1,MPI_DOUBLE,me+1,2,MPI_COMM_WORLD,&Status);
        }
      }

      std::vector<double> V_loc(h_part*Ny,0.);

      if (me!=0) {
        for (int j = 0; j < Ny; j++) {
          for (int i = 0; i < h_part; i++) {
            MPI_Send(&U_loc[j*N+i],1,MPI_DOUBLE,me-1,10,MPI_COMM_WORLD);
          }
        }
      }

      if (me!=Nproc-1) {
        double error_max(0.);
        for (int j=0; j<Ny; j++){
          for (int i=0; i<h_part; i++ ){
            MPI_Recv(&V_loc[i+j*h_part],1,MPI_DOUBLE,me+1,10,MPI_COMM_WORLD,&Status);
          }
        }

        for (int j=0; j<Ny; j++){
          for (int i=0; i<h_part; i++ ){
            error_loc = abs(U_loc[j*N+N-h_part+i]-V_loc[i+j*h_part]);
            if (error_loc>error_max) error_max=error_loc;
          }
        }
        error_loc=error_max;

        MPI_Send(&error_loc,1,MPI_DOUBLE,me+1,10,MPI_COMM_WORLD);
      }

      if (me!=0) {
        MPI_Recv(&error_loc,1,MPI_DOUBLE,me-1,10,MPI_COMM_WORLD,&Status);
      }

      iteschwz_loc++;
    }

    Write_(U_loc,IBeg,IEnd,Nx,Ny,N,Lx,Ly,h_part,k,me);

  }
  MPI_Finalize();

}

void Write_(std::vector<double> U_loc,int IBeg,int IEnd,int Nx,int Ny,int N,double Lx,double Ly,int h_part,int k,int me){

  std::string prefix="solution";
  std::ostringstream oss;
  oss<<prefix<<me<<".txt";
  std::ofstream NewFichier(oss.str().c_str(), std::ios::out | std::ios::trunc);

  std::vector<double> x_tab_loc,y_tab_loc;
  for (int i = IBeg; i < IEnd; i++) x_tab_loc.push_back(i*Lx/(Nx+1));
  for (int j = 0; j < Ny; j++) y_tab_loc.push_back(j*Ly/(Ny+1));

  for (int i = 0; i < N-1; i++){
    for (int j = 0; j < Ny; j++){
      NewFichier<< x_tab_loc[i] <<" "<<  y_tab_loc[j] << " " << U_loc[i+j*N] << std::endl;
    }
  }

}
