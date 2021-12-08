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
  double alpha(0), beta(0); // Coefficients de la condition d'interface (Neumann-Dirichlet)
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
  if (Nx%2==0) {
    Nu=Nx/2;
    Nv=Nx/2+h_part;
  }
  else if (Nx%2==1) {
    Nu=(Nx-1)/2;
    Nv=(Nx-1)/2+h_part;
  }

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
  //domaine 1
  Matrix(rowu,colu,valu,Nx,Ny,Nu,Nv,Lx,Ly,D,dt,alpha,beta,0);
  //domaine 2
  Matrix(rowv,colv,valv,Nx,Ny,Nu,Nv,Lx,Ly,D,dt,alpha,beta,1);


  std::string prefixe = "solution_approchee_seq_t=";

  //creation des abcisses et ordonnees des points du maillage
  std::vector<double> x_tab(Nx), y_tab(Ny);
  for (int i = 0; i < Nx; i++) x_tab[i]=(i+1)*Lx/(Nx+1);
  for (int j = 0; j < Ny; j++) y_tab[j]=(j+1)*Ly/(Ny+1);

  int Nt=100; // nombre d'iterations en temps
  double e=0.00000001; //tolerance pour le GC
  int kmax=Nx*Ny; //iteration max du BICG

  //Parametres d'arret pour Schwarz
  double errschwz = 0.00000001;
  int maxschwz = Nx*Ny;
  double error(1); int iteschwz(0);

  //schema en temps
  for (int k = 0; k < Nt; k++) {

    double t = k*dt; //temps de l'experience

    //Initialisation des stencils
    for (int i = 0; i < U0.size(); i++) U0[i]=1.,V0[i]=1.;

    //Boucle de Schwarz
    while ((iteschwz <= maxschwz)||(error >= errschwz)){
      //construction du second membre sur chaque sous-domaine
      std::vector<double> Su(Nu*Ny,0);
      std::vector<double> Sv(Nv*Ny,0);
      secondMembre(Su,U,V0,Nx,Ny,Nu,dt,t,Lx,Ly,D,mode,alpha,beta,0);
      secondMembre(Sv,V,U0,Nx,Ny,Nv,dt,t,Lx,Ly,D,mode,alpha,beta,1);

      //resolution du systeme lineaire sur chaque sous-domaine
      BICG(rowu,colu,valu,U,Su,e,kmax,Nu,Ny);
      BICG(rowv,colv,valv,V,Sv,e,kmax,Nv,Ny);

      //mise a jour des stencils
      for (int j = 0; j < Ny; j++){
        U0[j] = V[Ny*(h_part-1)+j];
        U0[Ny+j] = V[Ny*(h_part)+j];
        U0[2*Ny+j] = V[Ny*(h_part+1)+j];
        V0[j] = U[Ny*(Nu-h_part-1)+j];
        V0[Ny+j] = U[Ny*(Nu-h_part)+j];
        V0[2*Ny+j] = U[Ny*(Nu-h_part+1)+j];
      }

      // Evaluation de l'erreur pour la condition d'arret error
      error = maj_error(U,V,h_part,Ny,Nu);
      printf("error = %f\n", error );

      iteschwz++;
    }

    //ecriture dans les  fichiers

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
