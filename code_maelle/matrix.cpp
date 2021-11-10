#include <vector>
#include<cstdio>
#include<cstdlib>
#include"math.h"
#include "fonctions.hpp"
#include "matrix.hpp"


//------------------------------------------------------------------------------------------------------------------------
//fonction qui remplit la matrice de diffusion sur un sous domaine

void mul_Matrix(std::vector<double> &U, std::vector<double> &W, int Nx, int Ny, int N, double Lx, double Ly, double D, double dt, double alpha, double beta, int me){
  // Fonction qui remplit la matrice associée au problème avec un stockage coordonnées/valeur dans trois tableaux

  double dx=Lx/(Nx+1),dy=Ly/(Ny+1);

  double di(0),sdi(0),ssdi(0);
  di=1+2*dt*(1/pow(dx,2)+1/pow(dy,2));
  sdi=-dt*D/pow(dx,2);
  ssdi=-dt*D/pow(dy,2);

  if (me==0) {//domaine 1

    for (int j = 0; j < Ny; j++) {
      for (int i = 0; i < N; i++) {

        if (i==N-1) {//dernière ligne
          W[i+j*N]=(di+D*dt*beta/(alpha*dx))*U[i+j*N]+2*sdi*U[i-1+j*N]+ssdi*U[i+(j-1)*N]+ssdi*U[(j+1)*N];}

          else {//toute les autres lignes
            W[i+j*N]=di*U[i+j*N]+sdi*U[i+1+j*N]+sdi*U[i-1+j*N]+ssdi*U[i+(j-1)*N]+ssdi*U[(j+1)*N];}
          }
        }
      }

      else if (me==1) {//domaine 2

        for (int j = 0; j < Ny; j++) {
          for (int i = 0; i < N; i++) {

            if (i==0) {//première ligne
              W[i+j*N]=(di+D*dt*beta/(alpha*dx))*U[i+j*N]+2*sdi*U[i+1+j*N]+ssdi*U[i+(j-1)*N]+ssdi*U[(j+1)*N];}

              else {//toute les autres lignes
                W[i+j*N]=di*U[i+j*N]+sdi*U[i+1+j*N]+sdi*U[i-1+j*N]+ssdi*U[i+(j-1)*N]+ssdi*U[(j+1)*N];}
              }
            }
          }
        }


// ---------------------------------------------------------
//fonction qui remplit le vecteur second mambre du problème

void secondMembre_seq(std::vector<double> &S,std::vector<double> U, std::vector<double> V, int Nx, int Ny, int N, double dt,double t, double Lx, double Ly, double D, int mode, double alpha, double beta, int h_part, int me){

  // Fonction qui remplit le second membre selon le domaine où l'on se trouve

  double dx=Lx/(Nx+1),dy=Ly/(Ny+1); //pas d'espace

  for(int j=0; j<Ny; j++){
    for(int i=0; i<N; i++){

      S[i+j*Nx]=U[i+j*N]+dt*f((i+1)*dx,(j+1)*dy,t,Lx,Ly,mode);

      if (j==0) {
        S[i+j*N]+=D*dt/pow(dy,2)*g((i+1)*dx,0,mode);
      }

      else if (j==Ny-1) {
        S[i+j*N]+=D*dt/pow(dy,2)*g((i+1)*dx,Ny+1,mode);
      }

      else if (i==0) {
        S[i+j*N]+=D*dt/pow(dx,2)*h(0,(j+1)*dy,mode);
      }

      else if (i==Nx-1) {
        S[i+j*N]+=D*dt/pow(dx,2)*h(Nx+1,(j+1)*dy,mode);
      }
    }
  }
}
