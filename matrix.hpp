#ifndef _MATRIX_HPP_
#define _MATRIX_HPP_

void sparseMatrix(std::vector<int> &row,std::vector<int> &col,std::vector<double> &value, int Nx, int Ny, double Lx, double Ly, double D, double dt);
void secondMembre(std::vector<double> &S,std::vector<double> U, int Nx, int Ny, double dt,double t, double Lx, double Ly, double D,int mode);
void sparseMatrix_parallel(std::vector<int> &row,std::vector<int> &col,std::vector<double> &value, int Nx, int Ny, double Lx, double Ly, double D, double dt,int iBeg,int iEnd);
void secondMembre_parallel(std::vector<double> &S,std::vector<double> U, int Nx, int Ny,double dt, double t, double Lx, double Ly, double D, int mode,int iBeg,int iEnd);

#endif // _MATRIX_HPP_
