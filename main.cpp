#include <cstdio>
#include <vector>
#include <cstdlib>
#include <fstream>
#include <math.h>
#include <iostream>
#include <string>
#include <sstream>
#include <sys/time.h>
#include <mpi.h>

#include "fonctions.hpp"
#include "matrix.hpp"
#include "charge.hpp"
#include "solveur.hpp"
#include "update.hpp"

using namespace std;

//resolution de l'equation de diffusion sur un domaine Omega partitionne en deux sous domaine Omega1 et Omega2
//conditions aux limites sur les bords Gamma1 et Gamma2 de type non homogene donnees par des fonctions
//le domaine est scinde en deux (gauche et droite)

int main(int argc, char** argv){

  // Lecture des donnees dans le fichier parametres (parameters.txt)
  ifstream file("parameters.dat");
  int Nx(6), Ny(6);                // Nombre de noeuds en x et en y
  double Lx(1), Ly(1);             // Taille du rectangle
  double D(1), dt(0.01);              // Coefficient de diffusion et pas de temps
  int mode(1);                     // Mode pour le choix des fonctions f,g,h prises
  int h_part(2);                   // Parametre de recouvrement: nombre de lignes partagees par chaque sous-domaine
  double alpha(1), beta(1);        // Coefficients de la condition d'interface (Neumann-Dirichlet)
  file >> Nx >> Ny >> Lx >> Ly >> D >> dt >> mode >> h_part >> alpha >> beta ;
  file.close();

  int Nt=1; // nombre d'iterations en temps

  int sequential(1);
  if (argc>1) sequential = atoi(argv[1]);

  if (argc>2) h_part = atoi(argv[2]);
  if (argc>3) mode = atoi(argv[3]);
  if (argc>4) Nt = atoi(argv[4]);

  //Parametres d'arret pour le solveur
  double e=1e-10; //tolerance pour le CG
  int kmax=10*Nx*Ny; //iteration max du BICGStab

  //Parametres d'arret pour Schwarz
  double errschwz=1e-8; //tolerance pour Schwarz
  int maxschwz=10*Nx*Ny; //iteration max pour Schwarz

  struct timeval t1,t2;

  //----------------------------------------------------------------------
  // resolution sequentielle du probleme de diffusion par decomposition de domaine
  //----------------------------------------------------------------------

  if (sequential==1) {

    printf("\n!===== Resolution par decomposition de domaine =====! \n\n");
    printf("  Parametres maillage : Nx = %d, Ny = %d, Lx = %f , Ly = %f \n",Nx,Ny,Lx,Ly);
    printf("  Parametres simulation : D = %f, mode = %d, Nt = %d, dt = %f \n",D,mode,Nt,dt);
    printf("  Parametre decomposition domaine :  h_part = %d, alpha = %f, beta = %f \n",h_part,alpha,beta);

    gettimeofday(&t1, NULL);

    Update_dd(Nx,Ny,dt,Lx,Ly,D,mode,h_part,alpha,beta,Nt,e,kmax,errschwz,maxschwz);

    gettimeofday(&t2, NULL);

    printf("\n  temps d'execution : %lu \n\n",t2.tv_usec - t1.tv_usec);

  }

  //--------------------------------------------------------------------
  // fin de la resolution sequentielle
  //-------------------------------------------------------------------


  //--------------------------------------------------------------------
  // resolution parallele avec mpi
  //-------------------------------------------------------------------

  if (sequential==0) {

    printf("\n!===== Decomposition de domaines multithreads =====! \n\n");
    printf("  Parametres maillage : Nx = %d, Ny = %d, Lx = %f , Ly = %f \n",Nx,Ny,Lx,Ly);
    printf("  Parametres simulation : D = %f, mode = %d, Nt = %d, dt = %f \n",D,mode,Nt,dt);
    printf("  Parametre decomposition domaine :  h_part = %d, alpha = %f, beta = %f \n",h_part,alpha,beta);

    gettimeofday(&t1, NULL);

    Update_pll(argc, argv,Nx,Ny,dt,Lx,Ly,D,mode,h_part,alpha,beta,Nt,e,kmax,errschwz,maxschwz);

    gettimeofday(&t2, NULL);

    printf("\n  temps d'execution : %lu \n\n",t2.tv_usec - t1.tv_usec);

  }

  //--------------------------------------------------------------------
  // fin de la resolution parallele
  //-------------------------------------------------------------------

  return 0;

}
