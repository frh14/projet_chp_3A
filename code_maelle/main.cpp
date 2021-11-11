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
//résolution séquentielle du problème de diffusion
//----------------------------------------------------------------------------------------------------------------------------------------

  struct timeval t1,t2;
  gettimeofday(&t1, NULL);

  //construction de nos indices de séparation de domaine selon ses lignes
  int Nu(0); //domaine 1
  int Nv(0); //domaine 2
  //test sur la parité de Nx
  if (Nx%2==0) {
    Nu=Nx/2;
    Nv=Nx/2+h_part;}
  else if (Nx%2==1) {
    Nu=(Nx-1)/2;
    Nv=(Nx-1)/2+h_part;}

  //construction des vecteurs inconnues sur chaque sous domaine au temps initial:
  std::vector<double> U(Nu*Ny,1); //domaine 1
  std::vector<double> V(Nv*Ny,1); //domaine 2

  //construction des matrices de résolution
  std::vector<int> rowu, rowv;
  std::vector<int> colu, colv;
  std::vector<double> valu, valv;
  //domaine 1
  int me=0;
  Matrix(rowu,colu,valu,Nx,Ny,Nu,Nv,Lx,Ly,D,dt,alpha,beta,me);
  //domaine 2
  me=1;
  Matrix(rowv,colv,valv,Nx,Ny,Nu,Nv,Lx,Ly,D,dt,alpha,beta,me);

  /*//affichage pour phase de test
printf("Nx=%d et Ny=%d \n",Nx,Ny);
printf("Nu=%d et Nv=%d \n",Nu,Nv);
double dx=Lx/(Nx+1),dy=Ly/(Ny+1);
double di(0),sdi(0),ssdi(0);
di=1+2*dt*(1/pow(dx,2)+1/pow(dy,2));
sdi=-dt*D/pow(dx,2);
ssdi=-dt*D/pow(dy,2);
printf("------------------------------------------------------------------------\n");
printf("di=%f \n",di);
printf("sdi=%f \n",sdi);
printf("ssdi=%f \n",ssdi);
printf("di+D*dt*beta/(dx*alpha))=%f \n",di+D*dt*beta/(dx*alpha));
printf("2*sdi=%f \n",2*sdi);

printf("------------------------------------------------------------------------\n");
printf("------------------------------------------------------------------------\n");
printf("matrice domaine 1 de taille %d \n",Nu*Ny);
for (int i = 0; i < Nu*Ny; i++) {
  for (int k = 0; k < rowu.size(); k++) {
    if (rowu[k]==i) {
      printf(" ligne=%d , colonne=%d , valeur=%f \n",rowu[k],colu[k],valu[k]);}}}

printf("------------------------------------------------------------------------\n");
printf("matrice domaine 2 de taille %d \n",Nv*Ny);
for (int i = 0; i < Nv*Ny; i++) {
  for (int k = 0; k < rowv.size(); k++) {
      if (rowv[k]==i) {
          printf(" ligne=%d , colonne=%d , valeur=%f \n",rowv[k],colv[k],valv[k]);}}}*/

  int Nt=100; // nombre d'itérations en temps
  double e=0.00000001; //tolerance pour le GC
  int kmax=Nx*Ny; //itération max du GC

  std::string prefixe = "solution_approchée_seq_t=";

  //création des abcisses et ordonnées des points du maillage
  std::vector<double> x_tab(Nx),y_tab(Ny);
  for (int i = 0; i < Nx; i++) {
    x_tab[i]=(i+1)*Lx/(Nx+1);}
  for (int j = 0; j < Ny; j++) {
    y_tab[j]=(j+1)*Ly/(Ny+1);}

  //schema en temps
  for (int k = 0; k < Nt; k++) {

    double t=k*dt; //temps de l'expérience

    //construction du second membre sur chaque sous-domaines
    std::vector<double> S(Nx*Ny,0);
    //secondMembre(S,U,Nx,Ny,dt,t,Lx,Ly,D,mode);

    //résolution du systéme linéaire sur chaque sous-domaine

//écriture dans les  fichiers
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

  //-----------------------------------------------------------------------------------------------------------------------------------------
  //fin de la résolution séquentielle
  //------------------------------------

  return 0;

}
