#include <vector>
#include<cstdio>
#include<cstdlib>
#include"math.h"
#include "fonctions.hpp"
#include "matrix.hpp"

//------------------------------------------------------------------------------------------------------------------------
//fonction qui remplit la matrice de diffusion sur un sous domaine

void Matrix(std::vector<int> &row,std::vector<int> &col,std::vector<double> &val,int Nx, int Ny, int Nu, int Nv, double Lx, double Ly, double D, double dt, double alpha, double beta, int me){


  double dx=Lx/(Nx+1),dy=Ly/(Ny+1); //pas d'espace

  double di(0),sdi(0),ssdi(0);
  sdi=-dt*D/(dx*dx);
  ssdi=-dt*D/(dy*dy);
  di=1-2*dt*(sdi+ssdi);

  int N;
  if (me==0) N=Nu; //domaine 1
  else if (me==1 ) N=Nv;//domaine 2

  //construction de la matrice sur le domaine
  for(int i=0; i< N*Ny; i++){

    //remplissage des termes diagonaux du bloc diagonal
    if ((i+1)%N==0 && me==0) {//domaine 1 //modification de la derniere ligne
      col.push_back(i);
      row.push_back(i);

      val.push_back(di+D*dt*beta/(dx*alpha));
    }

    else if((i+1)%N==1 && me==1){//domaine 2 //modification de la premiere ligne
      col.push_back(i);
      row.push_back(i);
      val.push_back(di+D*dt*beta/(dx*alpha));
    }


    else{
      col.push_back(i);
      row.push_back(i);
      val.push_back(di);
    }

    if((i+1)>=N*Ny){} //test de la derniere ligne/colonne

    //remplissage des termes de la sur-diagonale et sous-diagonale du bloc diagonal
    else{

      col.push_back(i+1);//sur-diag inchangee quelque soit le domaine
      row.push_back(i);
      val.push_back(sdi);

      if ((i+1)%N==N-1 && me==0) {//domaine 1
        col.push_back(i);//sous-diag
        row.push_back(i+1);
        val.push_back(2*sdi);
      }

      else if (i%N==N-1 && me==1) {// domaine 2
        col.push_back(i);//sous-diag
        row.push_back(i+1);
        val.push_back(2*sdi);
      }

      else{
        col.push_back(i); //sous-diagonale
        row.push_back(i+1);
        val.push_back(sdi);
      }
    }

    //remplissage des blocs sur la sur-diagonale et la sous-diagonale par blocs
    if(i+Ny<N*Ny){
      col.push_back(i); //stockage de la sous-sous-diagonale
      row.push_back(i+Ny);
      val.push_back(ssdi);
      col.push_back(i+Ny); //stockage de la sur-sur-diagonale (symetrie)
      row.push_back(i);
      val.push_back(ssdi);
    }
  }


  // ---------------------------------------------------------
  //fonction qui remplit le vecteur second membre du probleme

  void secondMembre(std::vector<double> &S,std::vector<double> U, std::vector<double> V, int Nx, int Ny,int N, double dt,double t, double Lx, double Ly, double D, int mode, double alpha, double beta, int me){

    // Fonction qui remplit le second membre selon le domaine ou l'on se trouve

    double dx=Lx/(Nx+1),dy=Ly/(Ny+1); //pas d'espace
    double gamma=D*dt/(dx*dx);
    double eta=2*D*dt*beta/(dx*alpha);

    for(int j=0; j<Ny; j++){
      if (me==0){ //domaine 1
        for(int i=0; i<N; i++){

          S[i+j*Nu]=U[i+j*N]+dt*f((i+1)*dx,(j+1)*dy,t,Lx,Ly,mode);

          if (j==0) S[i+j*N]+=D*dt/dy*dy*g((i+1)*dx,0,mode);

          else if (j==Ny-1) S[i+j*N]+=D*dt/(dy*dy)*g((i+1)*dx,Ly,mode);

          else if (i==0) S[i+j*N]+=D*dt/dx*dx*h(0,(j+1)*dy,mode);

          else if (i==N-1) S[i+j*N]+=gamma*(V[Ny*2+j]-V[j]) + eta*V[Ny+j];
        }
      }
      else if (me==1){ //domaine 2
        // Invariant: Nu+Nv-h_part = Nx
        for(int i=0; i<N; i++){

          S[i+j*Nu]=U[i+j*N]+dt*f((Nx-N+i+1)*dx,(j+1)*dy,t,Lx,Ly,mode);

          if (j==0) S[i+j*N]+=D*dt/dy*dy*g((Nx-N+i+1)*dx,0,mode);

          else if (j==Ny-1) S[i+j*N]+=D*dt/(dy*dy)*g((Nx-N+i+1)*dx,Ly,mode);

          else if (i==0) S[i+j*N]+=gamma*(V[j]-V[Ny*2+j]) + eta*V[Ny+j];

          else if (i==N-1) S[i+j*N]+=D*dt/dx*dx*h(Lx,(j+1)*dy,mode);
        }
      }
    }

  }
