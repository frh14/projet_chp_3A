#ifndef _FONCTIONS_HPP_
#define _FONCTIONS_HPP_

double f (double x, double y, double t, double Lx, double Ly, int mode);
double g (double x, double y, int mode);
double h (double x, double y, int mode);
double maj_error(std::vector<double> U, std::vector<double> V, int h ,int Ny ,int Nu, int Nv);

#endif // _FONCTIONS_HPP_
