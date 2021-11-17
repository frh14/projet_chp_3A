#include "fonctions.hpp"

#include <math.h>

//module qui contient les fonctions utiles à la résolution du problème
//notamment les fonctions de conditions de limite


//----------------------------------------------------------------------------

double f(double x,double y,double t, double Lx, double Ly, int mode){
  if(mode==1){
    return 2*(y-(y*y)+x-(x*x));
  }
  if(mode==2){
    return sin(x)+cos(y);
  }
  if(mode==3){
    return exp(-(x-(Lx/2))*(x-(Lx/2)))*exp(-(y-(Ly/2))*(y-(Ly/2))*cos((2*atan(1))*t);
  }
}

// --------------------------------------------------------------------------

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

// -----------------------------------------------------------------------------

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
