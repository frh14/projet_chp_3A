#ifndef _MATRIX_HPP_
#define _MATRIX_HPP_

void Matrix(std::vector<int> &row,std::vector<int> &col,std::vector<double> &value,int Nx, int Ny, int N, double Lx, double Ly, double D, double dt, double alpha, double beta, int me);
void secondMembre(std::vector<double> &S,std::vector<double> U, std::vector<double> V, int Nx, int Ny, int N, double dt,double t, double Lx, double Ly, double D, int mode, double alpha, double beta, int me);

void Matrix_(std::vector<int> &row,std::vector<int> &col,std::vector<double> &value,int Nx, int Ny, int N, int Nproc, double Lx, double Ly, double D, double dt, double alpha, double beta, int me);
void secondMembre_(std::vector<double> &S,std::vector<double> U, std::vector<double> V1, std::vector<double> V2, int Nx, int Ny, int N, int Nproc, double dt,double t, double Lx, double Ly, double D, int mode, double alpha, double beta, int me);


#endif // _MATRIX_HPP_
