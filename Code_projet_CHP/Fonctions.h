#pragma once

#include "Data_File.h"
#include <cmath>
#include <vector>

void charge(int me, int n, int Np, int &ibeg, int &iend);
void print(const std::vector<double> &v, int me);
double f(const double x, const double y, const double t, Data_File *df);
double g(const double x, const double y, const double t, Data_File *df);
double h(const double x, const double y, const double t, Data_File *df);