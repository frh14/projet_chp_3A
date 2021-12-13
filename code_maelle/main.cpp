#include <cstdio>
#include <vector>
#include <cstdlib>
#include <fstream>
#include <math.h>
#include <iostream>
#include <string>
#include <sstream>
#include<sys/time.h>

#include "fonctions.hpp"
#include "matrix.hpp"
#include "charge.hpp"
#include "solveur.hpp"

using namespace std;

//résolution de l'equation de diffusion sur un domaine Omega partitionne avec en deux sous domaine Omega1 et Omega2
//conditions aux limites sur les bords Gamma1 et Gamma2 de type non homogene données par des fonctions
//le domaine est scinde en deux selon ses lignes

int main(int argc, char** argv){

  // Lecture des donnees dans le fichier parametres (parameters.txt)
  ifstream file("parameters.txt");
  int Nx(0), Ny(0);        // Nombre de noeuds en x et en y
  double Lx(0), Ly(0); // Taille du rectangle
  double D(0), dt(0); // Coefficient de diffusion et pas de temps
  int mode(0);  // Variable qui permet de choisir quelles fonctions f,g,h sont prises
  int h_part(0); // parametre de recouvrement: nombre de lignes partagees par chaque sous-domaine
  double alpha(1), beta(1); // Coefficients de la condition d'interface (Neumann-Dirichlet)
  file >> Nx >> Ny >> Lx >> Ly >> D >> dt >> mode >> h_part >> alpha >> beta ;
  file.close();

  //------------------------------------------------------------------------------------------------------------------------------------
  //resolution sequentielle du probleme de diffusion
  //----------------------------------------------------------------------------------------------------------------------------------------

  struct timeval t1,t2;
  gettimeofday(&t1, NULL);

  //construction de nos indices de separation de domaine selon ses lignes
  int Nu(0); //domaine 1
  int Nv(0); //domaine 2
  //test sur la parite de Nx
  if (Nx%2==0) Nu=Nx/2,Nv=Nx/2+h_part;

  else if (Nx%2==1) Nu=(Nx-1)/2,Nv=(Nx-1)/2+h_part;

  //construction des vecteurs inconnues sur chaque sous domaine au temps initial:
  std::vector<double> U(Nu*Ny,1); //domaine 1
  std::vector<double> V(Nv*Ny,1); //domaine 2

  //construction des stencils
  std::vector<double> U0(3*Ny); //domaine 1 pour le domaine 2
  std::vector<double> V0(3*Ny); //domaine 2 pour le domaine 1

  //construction des matrices de resolution sur chaque sous domaine
  std::vector<int> rowu, rowv;
  std::vector<int> colu, colv;
  std::vector<double> valu, valv;

  double dx=Lx/(Nx+1),dy=Ly/(Ny+1); //pas d'espace

  //construction matrices de résolution
  //domaine 1
  Matrix(rowu,colu,valu,Nx,Ny,Nu,Nv,Lx,Ly,D,dt,alpha,beta,0);

  double di(0),sdi(0),ssdi(0);
  sdi=-dt*D/(dx*dx);
  ssdi=-dt*D/(dy*dy);
  di=1-2*dt*(sdi+ssdi);

  //printf("di=%f,sdi=%f,ssdi=%f \n",di,sdi,ssdi);
  //for (int i = 0; i < rowu.size(); i++) printf("rowu[%d]=%d, colu[%d]=%d, valu[%d]=%f \n",i,rowu[i],i,colu[i],i,valu[i]);

  //domaine 2
  Matrix(rowv,colv,valv,Nx,Ny,Nu,Nv,Lx,Ly,D,dt,alpha,beta,1);
  //for (int i = 0; i < rowv.size(); i++) printf("rowv[%d]=%d, colv[%d]=%d, valv[%d]=%f \n",i,rowv[i],i,colv[i],i,valv[i]);

  int Nt=1; // nombre d'iterations en temps
  double e=0.00000000001; //tolerance pour le GC
  int kmax=Nx*Ny; //iteration max du BICG

  //Parametres d'arret pour Schwarz
  double errschwz = 0.000000001;
  int maxschwz = Nx*Ny;
  double error(1.); int iteschwz(0);

  //schema en temps
  for (int k = 0; k < Nt; k++) {

    double t = k*dt; //temps de l'experience

    //Initialisation des stencils
    for (int i = 0; i < U0.size(); i++) U0[i]=1.,V0[i]=1.;

    //Boucle de Schwarz
    while ((iteschwz <= maxschwz) && (error >= errschwz)){

      //construction du second membre sur chaque sous-domaine
      std::vector<double> Su(Nu*Ny,0);
      std::vector<double> Sv(Nv*Ny,0);

      secondMembre(Su,U,V0,Nx,Ny,Nu,dt,t,Lx,Ly,D,mode,alpha,beta,0);
      secondMembre(Sv,V,U0,Nx,Ny,Nv,dt,t,Lx,Ly,D,mode,alpha,beta,1);

/*      //test second membre domaine 1
      double gamma=D*dt/(dx*dx);
      double eta=2*D*dt*beta/(dx*alpha);
      std::vector<double> S(Nu*Ny,0);
      S[0]=U[0]+dt*f((0+1)*dx,(0+1)*dy,t,Lx,Ly,mode)+D*dt*h(0,(0+1)*dy,mode)/(dx*dx)+D*dt*g((0+1)*dx,0,mode)/(dy*dy);
      S[1]=U[1]+dt*f((1+1)*dx,(0+1)*dy,t,Lx,Ly,mode)+D*dt*h(0,(1+1)*dy,mode)/(dx*dx)+gamma*(V0[Ny*2+0]-V0[0]) + eta*V0[Ny+0];
      S[2]=U[2]+dt*f((0+1)*dx,(1+1)*dy,t,Lx,Ly,mode);
      S[3]=U[3]+dt*f((1+1)*dx,(1+1)*dy,t,Lx,Ly,mode)+gamma*(V0[Ny*2+1]-V0[1]) + eta*V0[Ny+1];
      S[4]=U[4]+dt*f((0+1)*dx,(2+1)*dy,t,Lx,Ly,mode);
      S[5]=U[5]+dt*f((1+1)*dx,(2+1)*dy,t,Lx,Ly,mode)+gamma*(V0[Ny*2+2]-V0[2]) + eta*V0[Ny+2];
      S[6]=U[6]+dt*f((0+1)*dx,(3+1)*dy,t,Lx,Ly,mode);
      S[7]=U[7]+dt*f((1+1)*dx,(3+1)*dy,t,Lx,Ly,mode)+gamma*(V0[Ny*2+3]-V0[3]) + eta*V0[Ny+3];
      S[8]=U[8]+dt*f((0+1)*dx,(4+1)*dy,t,Lx,Ly,mode)+D*dt*h(0,(0+1)*dy,mode)/(dx*dx);
      S[9]=U[9]+dt*f((1+1)*dx,(4+1)*dy,t,Lx,Ly,mode)+D*dt*h(0,(1+1)*dy,mode)/(dx*dx)+gamma*(V0[Ny*2+4]-V0[4]) + eta*V0[Ny+4];

      for (int i = 0; i < 10; i++) {
        printf("S[%d]=%f et Su[%d]=%f \n",i,S[i],i,Su[i]);
      }*/

/*
      S[i+j*N]+=dt*f((i+1)*dx,(j+1)*dy,t,Lx,Ly,mode);

      if (j==0) S[i+j*N]+=D*dt*g((i+1)*dx,0,mode)/(dy*dy);

      else if (j==Ny-1) S[i+j*N]+=D*dt*g((i+1)*dx,Ly,mode)/(dy*dy);

      if (i==0) S[i+j*N]+=D*dt*h(0,(j+1)*dy,mode)/(dx*dx);

      else if (i==N-1)S[i+j*N]+=gamma*(V[Ny*2+j]-V[j]) + eta*V[Ny+j];}
*/


    //  for (int i = 0; i < Su.size(); i++) {
        //printf("Su[%d]=%f \n",i,Su[i]);}

      /*for (int i = 0; i < Sv.size(); i++) {
        printf("Sv[%d]=%f \n",i,Sv[i]);}*/

      //resolution du systeme lineaire sur chaque sous-domaine
      BICG(rowu,colu,valu,U,Su,e,kmax,Nu,Ny);
      BICG(rowv,colv,valv,V,Sv,e,kmax,Nv,Ny);


  //    for (int i = 0; i < Su.size(); i++) {
        //printf("Su[%d]=%f \n",i,Su[i]);
    //  }

    /*  for (int i = 0; i < U.size(); i++) {
        printf("U[%d]=%f \n",i,U[i]);}*/

      //mise a jour des stencils

      for (int j = 0; j < Ny ; i++) {

        U0[j] = V[Nv*j];
        U0[Ny+j] = V[1+Nv*j];
        U0[2*Ny+j] = V[2+Nv*j];
        V0[j] = U[Nu-1+Nu*j];
        V0[Ny+j] =
        V0[2*Ny+j] = ;

      }


      for (int j = 0; j < Ny; j++){
        U0[j] = V[Ny*(h_part-1)+j];
        U0[Ny+j] = V[Ny*(h_part)+j];
        U0[2*Ny+j] = V[Ny*(h_part+1)+j];
        V0[j] = U[Ny*(Nu-h_part-1)+j];
        V0[Ny+j] = U[Ny*(Nu-h_part)+j];
        V0[2*Ny+j] = U[Ny*(Nu-h_part+1)+j];}

      // Evaluation de l'erreur pour la condition d'arret error
      error = maj_error(U,V,h_part,Ny,Nu);
      printf("erreur de schwartz=%f \n",error);

      iteschwz++;

    }

    //ecriture dans les  fichiers
    std::string prefixe="test";
    std::ostringstream oss;
    oss<<prefixe<<k<<".txt";
    std::ofstream NewFichier(oss.str().c_str(), std::ios::out | std::ios::trunc);

    std::vector<int> x_tab(Nx);
    std::vector<int> y_tab(Ny);

    for (int i = 0; i < Nx; i++) x_tab[i]=i*Lx/(Nx+1);
    for (int j = 0; j < Ny; j++) y_tab[j]=j*Ly/(Ny+1);

    for (int i = 0; i < Nu; i++)
      for (int j = 0; j < Ny; j++)
        NewFichier<< x_tab[i] <<" "<<  y_tab[j] << " " << U[i+j*Nu] << endl;

    for (int i = 0; i < Nv; i++)
      for (int j = 0; j < Ny; j++)
        NewFichier<< x_tab[i+Nu-h_part] <<" "<<  y_tab[j] << " " << V[i+j*Nv] << endl;

  }

  gettimeofday(&t2, NULL);

  printf("resolution sequentielle de parametres: \n");
  printf("Nx=%d, Ny=%d, Lx=%f , Ly=%f, D=%f, dt=%f , h_part=%d, alpha=%f, beta=%f \n",Nx,Ny,Lx,Ly,D,dt,h_part,alpha,beta);
  printf("temps ecoule =%lu \n",t2.tv_usec - t1.tv_usec);

  //-----------------------------------------------------------------------------------------------------------------------------------------
  //fin de la resolution sequentielle
  //------------------------------------

  return 0;

}
