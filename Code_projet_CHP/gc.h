#include <vector>
#include <cmath>
#include <fstream>
#include <string>
#include <mpi.h>
#include <iomanip>
#include <sstream>

#include "Data_File.h"
#include "Build.h"
#include "Fonctions.h"

//Fonctions de caculs surchargées
double norme(const std::vector<double> &x);
double ps(const std::vector<double> &x, const std::vector<double> &y);
std::vector<double> operator+(const std::vector<double> &x, const std::vector<double> &y);
std::vector<double> operator-(const std::vector<double> &x, const std::vector<double> &y);
std::vector<double> operator*(const double &real, const std::vector<double> &x);
void mpi_share_vec(const std::vector<double> &x, std::vector<double> &x_prod, int me, int Np, int ibeg, int iend, int Ny);
std::vector<double> prodmat(SparseMatrix &A, const std::vector<double> &x, int ibeg);
std::vector<double> operator*(SparseMatrix &A, const std::vector<double> &x);

//Sauvegarde de la solution u
void save_sol_gnuplot(std::vector<double> &u, double t, Data_File *df, int me, int ibeg, int iend);
//Gradient conjugué
void gc(SparseMatrix &A, std::vector<double> &u, std::vector<double> b, double t, int me, int Np, int ibeg, int iend, Data_File *df);
