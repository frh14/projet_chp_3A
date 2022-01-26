#ifndef _UPDATE_HPP_
#define _UPDATE_HPP_

void Update(std::vector<double> U,std::vector<double> V,int Nx,int Ny,int Nu,int Nv,double dt,double Lx,double Ly,double D,int mode,int h_part,double alpha,double beta,int Nt,double e,int kmax,double errschwz,int maxschwz);
void Write(std::vector<double> U,std::vector<double> V,int Nx,int Ny,int Nu,int Nv,double Lx,double Ly,int h_part,int k);
void Update_(int argc, char** argv,int Nx,int Ny,double dt,double Lx,double Ly,double D,int mode,int h_part,double alpha,double beta,int Nt,double e,int kmax,double errschwz,int maxschwz);
void Write_(std::vector<double> U_loc,int IBeg,int IEnd,int Nx,int Ny,int N,double Lx,double Ly,int h_part,int k,int me);

#endif // _UPDATE_HPP_
