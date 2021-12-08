#include "fonctions.hpp"
#include "solveur.hpp"

#include <math.h>

//module qui contient les fonctions utiles a la resolution du probleme
//notamment les fonctions de conditions de limite

// ---------------------------------------------------------------------------------------------

double f(double x,double y,double t, double Lx, double Ly, int mode){
  if(mode==1){
    return 2*(y-(y*y)+x-(x*x));
  }
  if(mode==2){
    return sin(x)+cos(y);
  }
  if(mode==3){
    return exp(-(x-(Lx/2))*(x-(Lx/2)))*exp(-(y-(Ly/2))*(y-(Ly/2))*cos((2*atan(1))*t));
  }
}

// ---------------------------------------------------------------------------------------------

double g(double x,double y,int mode){
  if(mode==1){
      return 0;
  }
  if(mode==2){
    return sin(x)+cos(y);
  }
  if(mode==3){
      return 0;
  }
}

// ---------------------------------------------------------------------------------------------

double h(double x,double y,int mode){
  if(mode==1){
      return 0;
  }
  if(mode==2){
    return sin(x)+cos(y);
  }
  if(mode==3){
      return 1;
  }
}

// ---------------------------------------------------------------------------------------------
//fonction qui calcule le max de l'erreur de la solution sur les points de recouvrement
//soit la difference de valeur de la solution sur les points communs aux deux sous domaines

double maj_error(std::vector<double> &U , std::vector<double> &V , int h , int Ny , int N){
  double error , error_max(0) ;
  std::vector<double> Vect_loc(Ny) ;

   for (int i=0 ; i< h ; i++ ){
     for (int j=0 ; j<Ny ; j++){
       Vect_loc[j] = U[(N-h+i)*Ny+j] - V[(i*Ny)+j];
       error = norm2(Vect_loc);
       if (error>error_max) error_max=error;
     }
   }
   return error_max ;
}
