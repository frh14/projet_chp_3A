#include <vector>
#include <cstdio>
#include <cstdlib>
#include <math.h>

#include "fonctions.hpp"
#include "matrix.hpp"

//------------------------------------------------------------------------------
//fonction qui remplit la matrice de diffusion sur un sous domaine

void Matrix(std::vector<int> &row,std::vector<int> &col,std::vector<double> &val,int Nx, int Ny, int N, double Lx, double Ly, double D, double dt, double alpha, double beta, int me){

  double dx=Lx/(Nx+1),dy=Ly/(Ny+1); //pas d'espace

  double di(0),sdi(0),ssdi(0), gamma(0);
  sdi=-dt*D/(dx*dx);
  ssdi=-dt*D/(dy*dy);
  di=1-2*(sdi+ssdi);

  if (alpha != 0) gamma = di+2*D*dt*beta/(dx*alpha);
  else gamma = 1./(dx*dx);

  //printf("sdi=%f , ssdi=%f , di=%f , gamma=%f \n",sdi,ssdi,di,gamma);

  for(int i=0; i<(N*Ny); i++){

    if (me==0) {
      if((i+1)%N==0){}
      else{

        if ((i+1)%N==N-1) {
          if (alpha==0){}
          else{
            col.push_back(i); //stockage de la sous-diagonale premiere ligne
            row.push_back(i+1);
            val.push_back(2*sdi);
          }
        }

        else if((i+1)%N==N){}

        else{
          col.push_back(i); //stockage de la sous-diagonale
          row.push_back(i+1);
          val.push_back(sdi);
        }

        col.push_back(i+1); //stockage de la sur-diagonale (symetrie)
        row.push_back(i);
        val.push_back(sdi);

      }

      if(i+N<N*Ny){
        if (alpha==0 && (i+N)%N==(N-1)){}
        else{
          col.push_back(i); //stockage de la sous-sous-diagonale
          row.push_back(i+N);
          val.push_back(ssdi);
        }
        if (alpha==0 && i%N==(N-1)){}
        else{
          col.push_back(i+N); //stockage de la sur-sur-diagonale (symetrie)
          row.push_back(i);
          val.push_back(ssdi);
        }
      }
    }

    else if (me==1) {
      if((i+1)%N==0){}
      else{
        if (i%N==0) {
          if (alpha==0){}
          else{
            col.push_back(i+1); //stockage de la sur-diagonale derniere ligne
            row.push_back(i);
            val.push_back(2*sdi);
          }
        }

        else{
          col.push_back(i+1); //stockage de la sur-diagonale (symetrie)
          row.push_back(i);
          val.push_back(sdi);
        }

        col.push_back(i); //stockage de la sous-diagonale
        row.push_back(i+1);
        val.push_back(sdi);

      }

      if(i+N<N*Ny){
        if (alpha==0 && (i+N)%N==0){}
        else{
          col.push_back(i); //stockage de la sous-sous-diagonale
          row.push_back(i+N);
          val.push_back(ssdi);
        }
        if (alpha==0 && i%N==0){}
        else{
          col.push_back(i+N); //stockage de la sur-sur-diagonale (symetrie)
          row.push_back(i);
          val.push_back(ssdi);
        }
      }
    }

    if ((i+1)%N==0 && me==0) { //domaine 1 : modification de la derniere ligne
      col.push_back(i),row.push_back(i),val.push_back(gamma);
    }
    else if((i+1)%N==1 && me==1){ //domaine 2 : modification de la premiere ligne
      col.push_back(i),row.push_back(i),val.push_back(gamma);
    }
    else{col.push_back(i),row.push_back(i),val.push_back(di);}

  }
}

// -----------------------------------------------------------------------------
//fonction qui remplit le vecteur second membre du probleme selon le domaine ou l'on se trouve

void secondMembre(std::vector<double> &S,std::vector<double> U, std::vector<double> V, int Nx, int Ny,int N, double dt,double t, double Lx, double Ly, double D, int mode, double alpha, double beta, int me){

  double dx=Lx/(Nx+1),dy=Ly/(Ny+1); //pas d'espace
  double delta=D*dt/(dx*dx),eta(1);

  if (alpha != 0) eta=2*D*dt*beta/(dx*alpha);
  else eta=1./(dx*dx);

  for(int j=0; j<Ny; j++){
    for(int i=0; i<N; i++){

      S[i+j*N]=U[i+j*N];

      if (me==0){ //domaine 1

        S[i+j*N]+=dt*f((i+1)*dx,(j+1)*dy,t,Lx,Ly,mode);

        if (j==0) S[i+j*N]+=D*dt*g((i+1)*dx,0,mode)/(dy*dy);

        if (j==Ny-1) S[i+j*N]+=D*dt*g((i+1)*dx,Ly,mode)/(dy*dy);

        if (i==0) S[i+j*N]+=D*dt*h(0,(j+1)*dy,mode)/(dx*dx);

        if (i==N-1){
          if (alpha==0) S[i+j*N]=eta*V[Ny+j];
          else S[i+j*N]+=delta*(V[2*Ny+j]-V[j])+eta*V[Ny+j];
        }
      }

      else if (me==1){ //domaine 2

        S[i+j*N]+=dt*f((Nx-N+i+1)*dx,(j+1)*dy,t,Lx,Ly,mode);

        if (j==0) S[i+j*N]+=D*dt*g((Nx-N+i+1)*dx,0,mode)/(dy*dy);

        if (j==Ny-1) S[i+j*N]+=D*dt*g((Nx-N+i+1)*dx,Ly,mode)/(dy*dy);

        if (i==0){
          if (alpha==0) S[i+j*N]=eta*V[Ny+j];
          else S[i+j*N]+=delta*(V[j]-V[2*Ny+j])+eta*V[Ny+j];
        }

        if (i==N-1) S[i+j*N]+=D*dt*h(Lx,(j+1)*dy,mode)/(dx*dx);
      }
    }
  }
}

//------------------------------------------------------------------------------
//fonction qui remplit la matrice de diffusion sur un proc

void Matrix_(std::vector<int> &row,std::vector<int> &col,std::vector<double> &val,int Nx, int Ny, int N, int Nproc, double Lx, double Ly, double D, double dt, double alpha, double beta, int me){

  double dx=Lx/(Nx+1),dy=Ly/(Ny+1); //pas d'espace

  double di(0),sdi(0),ssdi(0), gamma(0);
  sdi=-dt*D/(dx*dx);
  ssdi=-dt*D/(dy*dy);
  di=1-2*(sdi+ssdi);

  if (alpha != 0) gamma = di+2*D*dt*beta/(dx*alpha);
  else gamma = 1./(dx*dx);

  //printf("sdi=%f , ssdi=%f , di=%f , di+2*D*dt*beta/(dx*alpha)=%f \n",sdi,ssdi,di,di+2*D*dt*beta/(dx*alpha));

  for(int i=0; i<(N*Ny); i++){

    if (me==0) { // Proc 0

      if((i+1)%N==0){}
      else{
        if ((i+1)%N==N-1) {
          if (alpha==0){}
          else{
            col.push_back(i); //stockage de la sous-diagonale
            row.push_back(i+1);
            val.push_back(2*sdi);
          }
        }
        else if((i+1)%N==N){}
        else{
          col.push_back(i); //stockage de la sous-diagonale
          row.push_back(i+1);
          val.push_back(sdi);
        }

        col.push_back(i+1); //stockage de la sur-diagonale (symetrie)
        row.push_back(i);
        val.push_back(sdi);
      }

      if(i+N<N*Ny){
        col.push_back(i); //stockage de la sous-sous-diagonale
        row.push_back(i+N);
        val.push_back(ssdi);
        col.push_back(i+N); //stockage de la sur-sur-diagonale (symetrie)
        row.push_back(i);
        val.push_back(ssdi);
      }

      if ((i+1)%N==0) {//domaine 1 //modification de la derniere ligne
        col.push_back(i),row.push_back(i),val.push_back(gamma);
      }

      else{col.push_back(i),row.push_back(i),val.push_back(di);}

    }

    else if (me==Nproc-1) { // Proc Nproc-1

      if((i+1)%N==1){//domaine 2 //modification de la premiere ligne
        col.push_back(i),row.push_back(i),val.push_back(gamma);
      }

      else{col.push_back(i),row.push_back(i),val.push_back(di);}

      if((i+1)%N==0){}
      else{
        if (i%N==0) {//stockage de la sur-diag
          col.push_back(i+1); //stockage de la sur-diagonale (symetrie)
          row.push_back(i);
          val.push_back(2*sdi);
        }

        else{
          col.push_back(i+1); //stockage de la sur-diagonale (symetrie)
          row.push_back(i);
          val.push_back(sdi);
        }

        col.push_back(i); //stockage de la sous-diagonale
        row.push_back(i+1);
        val.push_back(sdi);
      }
    }

    else{
      if ((i+1)%N==0) {//domaine 1 //modification de la derniere ligne
        col.push_back(i),row.push_back(i),val.push_back(gamma);
      }

      else if((i+1)%N==1){//domaine 2 //modification de la premiere ligne
        col.push_back(i),row.push_back(i),val.push_back(gamma);
      }

      else{col.push_back(i),row.push_back(i),val.push_back(di);}

      if((i+1)%N==0){}
      else{
        if (i%N==0) {//stockage de la sur-diag
          col.push_back(i+1); //sur-diagonale
          row.push_back(i);
          val.push_back(2*sdi);
        }

        else{
          col.push_back(i+1); //sur-diagonale
          row.push_back(i);
          val.push_back(sdi);
        }

        if ((i+1)%N==N-1) {
          col.push_back(i); //sous-diagonale
          row.push_back(i+1);
          val.push_back(2*sdi);
        }

        else if((i+1)%N==N){}

        else{
          col.push_back(i); //sous-diagonale
          row.push_back(i+1);
          val.push_back(sdi);
        }
      }
    }
  }
}

// -----------------------------------------------------------------------------
//fonction qui remplit le vecteur second membre selon le proc

void secondMembre_(std::vector<double> &S,std::vector<double> U, std::vector<double> V1,std::vector<double> V2, int Nx, int Ny, int N, int Nproc,double dt,double t, double Lx, double Ly, double D, int mode, double alpha, double beta, int me){

  double dx=Lx/(Nx+1),dy=Ly/(Ny+1); //pas d'espace
  double delta=D*dt/(dx*dx);
  double eta=2*D*dt*beta/(dx*alpha);

  for(int j=0; j<Ny; j++){
    for(int i=0; i<N; i++){

      S[i+j*N]=U[i+j*N];

      if (me==0){ // Proc 0

        S[i+j*N]+=dt*f((i+1)*dx,(j+1)*dy,t,Lx,Ly,mode);

        if (j==0) S[i+j*N]+=D*dt*g((i+1)*dx,0,mode)/(dy*dy);

        if (j==Ny-1) S[i+j*N]+=D*dt*g((i+1)*dx,Ly,mode)/(dy*dy);

        if (i==0) S[i+j*N]+=D*dt*h(0,(j+1)*dy,mode)/(dx*dx);

        if (i==N-1)S[i+j*N]+=delta*(V1[2*Ny+j]-V1[j])+eta*V1[Ny+j];
      }

      else if (me==Nproc-1){ // Proc Nproc-1

        S[i+j*N]+=dt*f((Nx-N+i+1)*dx,(j+1)*dy,t,Lx,Ly,mode);

        if (j==0) S[i+j*N]+=D*dt*g((Nx-N+i+1)*dx,0,mode)/(dy*dy);

        if (j==Ny-1) S[i+j*N]+=D*dt*g((Nx-N+i+1)*dx,Ly,mode)/(dy*dy);

        if (i==0) S[i+j*N]+=delta*(V2[j]-V2[2*Ny+j])+eta*V2[Ny+j];

        if (i==N-1) S[i+j*N]+=D*dt*h(Lx,(j+1)*dy,mode)/(dx*dx);}

        else{

          S[i+j*N]+=dt*f((Nx-N+i+1)*dx,(j+1)*dy,t,Lx,Ly,mode);

          if (j==0) S[i+j*N]+=D*dt*g((Nx-N+i+1)*dx,0,mode)/(dy*dy);

          if (j==Ny-1) S[i+j*N]+=D*dt*g((Nx-N+i+1)*dx,Ly,mode)/(dy*dy);

          if (i==0) S[i+j*N]+=delta*(V2[j]-V2[2*Ny+j])+eta*V2[Ny+j];

          if (i==N-1) S[i+j*N]+=delta*(V1[2*Ny+j]-V1[j])+eta*V1[Ny+j];
        }
      }
    }
  }
