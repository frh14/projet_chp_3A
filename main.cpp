#include <cstdio>
#include <vector>
#include <cstdlib>
#include <fstream>
#include <math.h>
#include <iostream>
#include <string>
#include <sstream>
#include <sys/time.h>

#include "charge.hpp"
#include "update.hpp"

using namespace std;

//resolution de l'equation de diffusion sur un domaine Omega partitionne en deux sous domaine Omega1 et Omega2
//conditions aux limites sur les bords Gamma1 et Gamma2 de type non homogene donnees par des fonctions
//le domaine est scinde en deux (gauche et droite)

int main(int argc, char** argv){

  // Lecture des donnees dans le fichier parametres (parameters.txt)
  ifstream file("parameters.dat");
  int Nx(0), Ny(0);        // Nombre de noeuds en x et en y
  double Lx(0), Ly(0); // Taille du rectangle
  double D(0), dt(0); // Coefficient de diffusion et pas de temps
  int mode(0);  // Variable qui permet de choisir quelles fonctions f,g,h sont prises
  int h_part(0); // parametre de recouvrement: nombre de lignes partagees par chaque sous-domaine
  double alpha(1), beta(1); // Coefficients de la condition d'interface (Neumann-Dirichlet)
  file >> Nx >> Ny >> Lx >> Ly >> D >> dt >> mode >> h_part >> alpha >> beta ;
  file.close();

//----------------------------------------------------------------------
//resolution sequentielle du probleme de diffusion
//----------------------------------------------------------------------

  int Nt=1; // nombre d'iterations en temps
  double e=1e-10; //tolerance pour le GC
  int kmax=Nx*Ny; //iteration max du BICGStab

  //Parametres d'arret pour Schwarz
  double errschwz = 1e-8;
  int maxschwz = 100;
  
  //construction de nos indices de separation de domaine selon ses lignes
  int Nu(0); //domaine 1
  int Nv(0); //domaine 2
  //test sur la parite de Nx
  if (Nx%2==0) Nu=Nx/2,Nv=Nx/2+h_part;
  else if (Nx%2==1) Nu=(Nx-1)/2,Nv=(Nx+1)/2+h_part;

  printf("Nx=%d , Nu=%d , Nv=%d \n",Nx,Nu,Nv);

  //construction des vecteurs inconnues sur chaque sous domaine au temps initial:
  std::vector<double> U(Nu*Ny,1); //domaine 1
  std::vector<double> V(Nv*Ny,1); //domaine 2
  
  struct timeval t1,t2;
  gettimeofday(&t1, NULL);
  
  Update(U,V,Nx,Ny,Nu,Nv,Lx,Ly,D,dt,mode,h_part,alpha,beta,Nt,e,kmax,errschwz,maxschwz);

  gettimeofday(&t2, NULL);

  printf("resolution sequentielle de parametres: \n");
  printf("Nx=%d, Ny=%d, Lx=%f , Ly=%f, D=%f, dt=%f , h_part=%d, alpha=%f, beta=%f \n",Nx,Ny,Lx,Ly,D,dt,h_part,alpha,beta);
  printf("temps ecoule =%lu \n",t2.tv_usec - t1.tv_usec);

  //--------------------------------------------------------------------
  //fin de la resolution sequentielle
  //-------------------------------------------------------------------

  return 0;

}
