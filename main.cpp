#include <cstdio>
#include <vector>
#include <cstdlib>
#include <fstream>
#include <math.h>
#include <iostream>
#include <string>
#include <sstream>
#include<sys/time.h>
#include <mpi.h>

#include "fonctions.hpp"
#include "matrix.hpp"
#include "charge.hpp"
#include "solveur.hpp"

using namespace std;

int main(int argc, char** argv) {

  // Lecture des données dans le fichier paramètres (parameters.txt)
  ifstream file("parameters.txt");
  int    Nx(0), Ny(0);        // Nombre de noeuds en x et en y
  double Lx(0), Ly(0), D(0), dt(0); //taille du rectangle, coefficient de diffusion et pas de temps
  int    mode(0);             // Variable qui permet de choisir quelles fonctions f,g,h sont prises
  file >> Nx >> Ny >> Lx >> Ly >> D >> dt >> mode;
  file.close();


//------------------------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------------------

struct timeval t1,t2;
gettimeofday(&t1, NULL);

//résolution séquentielle du problème de diffusion

// Construction de la matrice de diffusion associée au problème
//La matrice est creuse // stockage adapté
//vecteur valeur->contient la valeur des termes non nuls de la matrice
//vecteur ligne et colonne->contiennent les indices des termes non nuls
  std::vector<int> row,col;
  std::vector<double> value;
  sparseMatrix(row,col,value,Nx,Ny,Lx,Ly,D,dt);

//initialisation du problème

  int N=100; // nombre d'itérations en temps
  double e=0.00000001; //tolerance pour le GC
  int kmax=Nx*Ny; //itération max du GC
  std::vector<double> U(Nx*Ny,1); //vecteur inconnu initial // point de départ du GC

  std::string prefixe = "solution_approchée_seq_t=";

//schema en temps
//création des abcisses et ordonnées des points du maillage
std::vector<double> x_tab(Nx),y_tab(Ny);
for (int i = 0; i < Nx; i++) {
  x_tab[i]=(i+1)*Lx/(Nx+1);}
for (int j = 0; j < Ny; j++) {
  y_tab[j]=(j+1)*Ly/(Ny+1);}


  for (int k = 0; k < N; k++) {

    double t=k*dt; //temps de l'expérience

//construction du second membre
    std::vector<double> S(Nx*Ny,0);
    secondMembre(S,U,Nx,Ny,dt,t,Lx,Ly,D,mode);

//résolution du systéme linéaire à l'aide du GC
    GC(row,col,value,U,S,e,kmax,Nx,Ny);

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
printf("Nx=%d, Ny=%d, Lx=%f , Ly=%f, D=%f, dt=%f \n",Nx,Ny,Lx,Ly,D,dt);
printf("temps écoulé =%lu \n",t2.tv_usec - t1.tv_usec);

  return 0;

}
