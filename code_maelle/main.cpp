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
  int    mode(0) , nb_shared_col(0);             // Variable qui permet de choisir quelles fonctions f,g,h sont prises , nb de colonne partagee
  file >> Nx >> Ny >> Lx >> Ly >> D >> dt >> mode >> nb_shared_col;
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

//--------------------------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------------------------

struct timeval t3,t4;

//résolution du probléme de diffusion par un code paralléle
  int N_loc=100;
  double e_loc=0.00000001; //tolerance pour le GC
  int kmax_loc=Nx*Ny; //nombre d'itération max du GC
  double t_loc(0); // temps initial
  std::vector<double> p_glob(Nx*Ny);
  std::vector<double> U_glob(Nx*Ny);




  MPI_Init(&argc,&argv);

  int me,Nproc;
  MPI_Comm_rank(MPI_COMM_WORLD,&me);  //processeur local
  MPI_Comm_size(MPI_COMM_WORLD,&Nproc); //nombre de processeur total utilisés

  if (me==0) {
    gettimeofday(&t3, NULL);
  }

//distribution de la charge entre les différents processeurs
  int iBeg,iEnd;
  charge(Nx*Ny,Nproc,me,iBeg,iEnd);

//construction de la matrice de diffusion par morceaux pour chaque proc
//matrice creuse // stockage adpaté
//vecteur valeur->contient la valeur des termes non nuls de la matrice
//vecteur ligne et colonne->contiennent les indices des termes non nuls
  std::vector<int> row_loc,col_loc;
  std::vector<double> value_loc;
  sparseMatrix_parallel(row_loc,col_loc,value_loc,Nx,Ny,Lx,Ly,D,dt,iBeg,iEnd);

  int loc_size=iEnd-iBeg+1; //taille locale des vecteurs

  std::vector<double> U_loc(loc_size,1); //vecteur inconnu initial // point de départ du GC

  //std::string prefixe2 = "par_solution_approchée_t=";

//schéma en temps

for (int j = 0; j < N_loc; j++) {

  t_loc=j*dt; //temps de l'expérience

  //construction du vecteur second membre local pour chaque proc
  std::vector<double> S_loc(loc_size,0);
  secondMembre_parallel(S_loc,U_loc,Nx,Ny,dt,t_loc,Lx,Ly,D,mode,iBeg,iEnd);

  //partie communication entre les processeurs:
  //transmission des vecteurs locaux U loc d'inconnues nécessaires aux autres processeurs
  for (int i = iBeg; i < iEnd+1; i++) {
    U_glob[i]=U_loc[i-iBeg];
  }

  for (int i = 0; i < Nproc; i++) {
      int iDeb,iFin;
      charge(Nx*Ny,Nproc,i,iDeb,iFin);
      MPI_Bcast(&U_glob[iDeb],iFin-iDeb+1,MPI_DOUBLE,i,MPI_COMM_WORLD);
  }

  std::vector<int> U_indice;

  for (int i = 0; i < Nx*Ny; i++) {
    U_indice.push_back(i);
  }

  //réalisation du produit matrice vecteur z=Ap
  std::vector<double> z(loc_size,0);

  mulSparseMatrix_2(row_loc,col_loc,value_loc,z,iBeg,iEnd,U_glob,U_indice);

//calul de la premiere direction de descente p
    std::vector<double> pprime(loc_size,0);

    mulSparseMatrix_2(row_loc,col_loc,value_loc,pprime,iBeg,iEnd,U_glob,U_indice);

    std::vector<double> p(loc_size,0);
    for (int i = 0; i < loc_size; i++) {
      p[i]=S_loc[i]-pprime[i];}

//remplissage du reste initial r0
    std::vector<double> r(loc_size,0);
    for (int i = 0; i < loc_size; i++) {
      r[i]=p[i];}

//calcul du beta global (norme du reste) global par réduction
   double beta_loc=pow(norme2(r),2);
   double beta;
   MPI_Allreduce(&beta_loc,&beta,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
   beta=sqrt(beta);

//initialisation
    int k=0;

//début de la méthode itérative

      while (beta>e_loc && k<=kmax_loc) {

//remplissage de p avec les valeurs utiles à un produit matrice/vecteur Az=p

        for (int i = iBeg; i < iEnd+1; i++) {
          p_glob[i]=p[i-iBeg];
        }

        for (int i = 0; i < Nproc; i++) {
            int iDeb,iFin;
            charge(Nx*Ny,Nproc,i,iDeb,iFin);
            MPI_Bcast(&p_glob[iDeb],iFin-iDeb+1,MPI_DOUBLE,i,MPI_COMM_WORLD);
        }

    //réalisation du produit matrice vecteur z=Ap
        std::vector<double> z(loc_size,0);
        std::vector<int> p_indice;

        for (int i = 0; i < Nx*Ny; i++) {
          p_indice.push_back(i);
        }

        mulSparseMatrix_2(row_loc,col_loc,value_loc,z,iBeg,iEnd,p_glob,p_indice);

 //calcul du pas de descente par réduction des produits scaliares locaux

        double w_loc=ps(r,r);
        double w;
        MPI_Allreduce(&w_loc,&w,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

        double y_loc=ps(z,p);
        double y;
        MPI_Allreduce(&y_loc,&y,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

        double alpha=w/y; //pas de descente global

//mise à jour des valeurs des composantes du vecteur x
        for (int i = iBeg; i < iEnd+1; i++) {
            U_glob[i]+=alpha*p_glob[i];}

//calcul du nouveau reste
        std::vector<double> rplus(loc_size,0);
        for (int i = 0; i < loc_size; i++) {
          rplus[i]=r[i]-alpha*z[i];}

//mise à jour de la nouvelle direction de descente
   //réduction des ps locaux
        double gamma_loc=ps(rplus,rplus);
        double gamma;
        MPI_Allreduce(&gamma_loc,&gamma,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        gamma=gamma/w;

   //calcul de la direction de descente
        for (int i=0; i < loc_size; i++) {
          p[i]=rplus[i]+gamma*p[i];}

//mise à jour du reste
        for (int i = 0; i < loc_size; i++) {
          r[i]=rplus[i];}

//mise à jour de beta global norme du reste par réduction des beta locaux
        double beta_loc2=pow(norme2(r),2);
        MPI_Allreduce(&beta_loc2,&beta,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        beta=sqrt(beta);

//incrémentation
        k+=1;}

//mise à jour du vecteur inconnu local avec le résultat du GC
    for (int i = 0; i < loc_size; i++) {
      U_loc[i]=U_glob[i+iBeg];}

    if (me==0) {
      gettimeofday(&t4, NULL);
    }


//écriture dans un fichier de ce vecteur local
//A ENLEVER POUR CALCULER LES SPEED-UP
   /*std::string prefixe2 = "solution_approchée_par_t=";
   std::ostringstream oss2;
   oss2<<prefixe2<<j<<"proc="<<me<<".txt";
   std::ofstream NewFichier2(oss2.str().c_str(), std::ios::out | std::ios::trunc);

  for (int k = 0; k < loc_size; k++) {
    int i=(k+iBeg)/Nx;
    int j=(k+iBeg)/Nx;

    NewFichier2<< x_tab[i] <<" "<<  y_tab[j] << " " << U_loc[k] << endl;}*/


  } //fin des itérations en temps

  MPI_Finalize();


//fin de la zone parallèle

  printf("résolution parallèle de paramètres: \n");
  printf("Nx=%d, Ny=%d, Lx=%f , Ly=%f, D=%f, dt=%f \n",Nx,Ny,Lx,Ly,D,dt);
  printf("temps de calcul=%lu \n",t4.tv_usec - t3.tv_usec);
  printf("speed up=%f \n",(double)(t2.tv_usec - t1.tv_usec)/(t4.tv_usec - t3.tv_usec) );

  return 0;

}
