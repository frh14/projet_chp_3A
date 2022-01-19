#ifndef _SOLVEUR_HPP_
#define _SOLVEUR_HPP_
#include <vector>

std::vector<double> BICGStabTest(std::vector<int> row,std::vector<int> col,std::vector<double> val,std::vector<double> &F,double e,int kmax,int Nx,int Ny);
void BICGStab(std::vector<int> row,std::vector<int> col,std::vector<double> val,std::vector<double> &X,std::vector<double> &F, double e, int kmax,int Nx,int Ny);
double norm2(std::vector<double> x);
double ps(std::vector<double> x,std::vector<double> y);
void mulSparseMatrix(std::vector<int> row, std::vector<int> col, std::vector<double> val, std::vector<double> &y, std::vector<double> x);
void transpSparseMatrix(std::vector<int> row, std::vector<int> col, std::vector<double> val,std::vector<int> &rowT, std::vector<int> &colT, std::vector<double> &valT );

#endif
