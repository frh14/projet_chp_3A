#ifndef _SOLVEUR_HPP_
#define _SOLVEUR_HPP_
#include <vector>

void BICGStab(std::vector<int> &row,std::vector<int> &col,std::vector<double> &val,std::vector<double> &X,std::vector<double> &F, double e, int kmax,int Nx,int Ny);
double norm2(std::vector<double> &x);
double ps(std::vector<double> &x,std::vector<double> &y);
void mulSparseMatrix(std::vector<int> &row, std::vector<int> &col, std::vector<double> &val, std::vector<double> &y, std::vector<double> &x);

#endif
