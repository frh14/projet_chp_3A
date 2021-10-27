#ifndef _SOLVEUR_HPP_
#define _SOLVEUR_HPP_

void GC(std::vector<int> row,std::vector<int> col, std::vector<double> value, std::vector<double> &Uo, std::vector<double> F,double e,int kmax,int Nx,int Ny);
double norme2(std::vector<double> x);
double ps(std::vector<double> x,std::vector<double> y);
void mulSparseMatrix(std::vector<int> row, std::vector<int> col, std::vector<double> value, std::vector<double> &y, std::vector<double> x);
void mulSparseMatrix_2(std::vector<int> row, std::vector<int> col, std::vector<double> value, std::vector<double> &y,int iBeg,int iEnd,std::vector<double> x_val,std::vector<int> x_indice);

#endif
