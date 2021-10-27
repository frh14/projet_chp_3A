#include "Fonctions.h"
void charge(int me, int n, int Np, int &ibeg, int &iend)
{
    int size(n / Np);
    int r(n % Np);

    ibeg = size * me;
    iend = ibeg + size - 1;

    if (me < r)
    {
        ibeg += me;
        iend += me + 1;
    }
    else
    {
        ibeg += r;
        iend += r;
    }
}

void print(const std::vector<double> &v, int me)
{
    for (int i = 0; i < v.size(); i++)
        std::cout << "vec " << me << " " << v[i] << std::endl;
}

double f(const double x, const double y, const double t, Data_File *df)
{
    switch (df->get_conditions())
    {
    case 1:
        return 2. * (y - y * y + x - x * x);
        break;
    case 2:
        return sin(x) + cos(y);
        break;
    case 3:
        return exp(-pow((x - df->get_Lx() / 2), 2)) * exp(-pow((y - df->get_Ly() / 2), 2)) * cos(M_PI / 2 * t);
        break;
    default:
        return 0.;
        break;
    }
}

double g(const double x, const double y, const double t, Data_File *df)
{
    switch (df->get_conditions())
    {
    case 1:
        return 0;
        break;
    case 2:
        return sin(x) + cos(y);
        break;
    case 3:
        return 0.;
        break;
    default:
        return 0.;
        break;
    }
}

double h(const double x, const double y, const double t, Data_File *df)
{
    switch (df->get_conditions())
    {
    case 1:
        return 0;
        break;
    case 2:
        return sin(x) + cos(y);
        break;
    case 3:
        return 1.;
        break;
    default:
        return 0.;
        break;
    }
}

