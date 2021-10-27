#include "charge.hpp"

void charge(int n, int Np, int me, int& is, int& ie){
  int r = n % Np;
  if (me < r){
    is = me*(n/Np+1);
    ie = is + (n/Np+1) - 1;
  }
  else {
    is = r + me*(n/Np);
    ie = is + (n/Np) - 1;
  }
}
