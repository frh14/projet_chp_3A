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

//résolution de l'équation de diffusion sur un domaine Omega partitionné avec en deux sous domaine Omega1 et Omega2
//conditions aux limites sur les bords Gamma1 et Gamma2 de type non homogéne données par des fonctions
//le domaine est scindé en deux selon ses lignes



int main(int argc, char** argv) {

  // Lecture des données dans le fichier paramètres (parameters.txt)
  ifstream file("parameters.txt");
  int Nx(0), Ny(0);        // Nombre de noeuds en x et en y
  double Lx(0), Ly(0); //Taille du rectangle
  double D(0), dt(0); // Coefficient de diffusion et pas de temps
  int mode(0);  // Variable qui permet de choisir quelles fonctions f,g,h sont prises
  int h_part(0); // paramétre de recouvrement: nombre de lignes partagées par chaque sous-domaine
  double alpha(0), beta(0); // Coefficients de la condition d'interface (Neumann-Dirichlet)
  file >> Nx >> Ny >> Lx >> Ly >> D >> dt >> mode >> h_part >> alpha >> beta ;
  file.close();


//------------------------------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------------------------------

struct timeval t1,t2;
gettimeofday(&t1, NULL);

//résolution séquentielle du problème de diffusion

//construction de nos indices de séparation de domaine

int Nu(0); //domaine 1
int Nv(0); //domaine 2

//test sur la parité de Nx
if (Nx%2==0) {
  Nu=Nx/2;
  Nv=Nx/2+h_part;}
else if (Nx%2==1) {
  Nu=(Nx-1)/2;
  Nv=(Nx-1)/2+h_part;}

//construction des vecteurs inconnues sur chaque sous domaine:

std::vector<double> U(Nu*Ny,0); //domaine 1
std::vector<double> V(Nv*Ny,0); //domaine 2

//construction des vecteurs du produit matrice-vecteur
std::vector<double> Wu(Ny,0); //domaine 1
std::vector<double> Wv(Ny,0); //domaine 2








//initialisation du problème

  int Nt=100; // nombre d'itérations en temps
  double e=0.00000001; //tolerance pour le GC
  int kmax=Nx*Ny; //itération max du GC
  //std::vector<double> U(Nx*Ny,1); //vecteur inconnu initial // point de départ du GC

  std::string prefixe = "solution_approchée_seq_t=";

//schema en temps
//création des abcisses et ordonnées des points du maillage
std::vector<double> x_tab(Nx),y_tab(Ny);
for (int i = 0; i < Nx; i++) {
  x_tab[i]=(i+1)*Lx/(Nx+1);}
for (int j = 0; j < Ny; j++) {
  y_tab[j]=(j+1)*Ly/(Ny+1);}


  for (int k = 0; k < Nt; k++) {

    double t=k*dt; //temps de l'expérience

//construction du second membre
    std::vector<double> S(Nx*Ny,0);
    //secondMembre(S,U,Nx,Ny,dt,t,Lx,Ly,D,mode);

//résolution du systéme linéaire à l'aide du GC
    //GC(row,col,value,U,S,e,kmax,Nx,Ny);

//écriture dans les  fichiers
//A ENLEVER POUR CALCULER LES SPEED-UP
/*
 std::ostringstream oss;
 oss<<prefixe<<k<<".txt";
 std::ofstream NewFichier(oss.str().c_str(), std::ios::out | std::ios::trunc);

 for (int i = 0; i < Nx; i++) {
   for (int j = 0; j < Ny; j++) {
     NewFichier<< x_tab[i] <<" "<<  y_tab[j] << " " << U[i+j*Nx] << endl;}}*/

  }
  gettimeofday(&t2, NULL);

printf("résolution séquentielle de paramètres: \n");
printf("Nx=%d, Ny=%d, Lx=%f , Ly=%f, D=%f, dt=%f , h_part=%d, alpha=%f, beta=%f \n",Nx,Ny,Lx,Ly,D,dt,h_part,alpha,beta);
printf("temps écoulé =%lu \n",t2.tv_usec - t1.tv_usec);

  return 0;

}
