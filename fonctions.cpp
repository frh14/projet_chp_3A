#include <vector>
#include <cstdio>
#include <cstdlib>
#include <math.h>

#include "fonctions.hpp"
#include "solveur.hpp"

//module qui contient les fonctions utiles a la resolution du probleme
//notamment les fonctions de conditions de limite

//------------------------------------------------------------------------------

double f(double x,double y,double t, double Lx, double Ly, int mode){
  if(mode==1) return 2*(y-(y*y)+x-(x*x));
  else if(mode==2) return sin(x)+cos(y);
  else if(mode==3) return exp(-(x-(Lx/2))*(x-(Lx/2)))*exp(-(y-(Ly/2))*(y-(Ly/2))*cos((2*atan(1))*t));
  else{
    printf("Le mode choisit n'existe pas!\n");
    return 0;
  }
}

// -----------------------------------------------------------------------------

double g(double x,double y,int mode){
  if(mode==1) return 0;
  else if(mode==2) return sin(x)+cos(y);
  else if(mode==3) return 0;
  else{
    printf("Le mode choisit n'existe pas!\n");
    return 0;
  }
}

// -----------------------------------------------------------------------------

double h(double x,double y,int mode){
  if(mode==1) return 0;
  else if(mode==2) return sin(x)+cos(y);
  else if(mode==3) return 1;
  else{
    printf("Le mode choisit n'existe pas!\n");
    return 0;
  }
}

// -----------------------------------------------------------------------------

double maj_error(std::vector<double> U, std::vector<double> V, int h, int Ny, int Nu, int Nv){
  double error , error_max(0) ;
  std::vector<double> Vect_loc(Ny) ;

  for (int j=0; j<Ny; j++){
    for (int i=0; i<h; i++ ){
      error = abs(U[j*Nu+(Nu-h-1)+i] - V[j*Nv+i]);
      if (error>error_max) error_max=error;
    }
  }
  return error_max;
}
