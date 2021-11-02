#include "gc.h"

double norme(const std::vector<double> &x)
{
    return sqrt(ps(x, x));
}

double ps(const std::vector<double> &x, const std::vector<double> &y)
{
    double sum(0.), sum_local(0.);

    for (int k = 0; k < x.size(); k++)
        sum_local += x[k] * y[k];
    int Np;
    MPI_Comm_size(MPI_COMM_WORLD, &Np); // Nombre de Procs

    if (Np != 1)
        MPI_Allreduce(&sum_local, &sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    if (Np == 1)
        return sum_local;

    return sum;
}

std::vector<double> operator+(const std::vector<double> &x, const std::vector<double> &y)
{
    std::vector<double> c(x.size(), -12.);
    for (int i = 0; i < x.size(); i++)
        c[i] = x[i] + y[i];

    return c;
}

std::vector<double> operator-(const std::vector<double> &x, const std::vector<double> &y)
{
    std::vector<double> c(x.size(), -12.);
    for (int i = 0; i < x.size(); i++)
        c[i] = x[i] - y[i];

    return c;
}

std::vector<double> operator*(const double &real, const std::vector<double> &x)
{
    std::vector<double> c(x.size(), -12.);
    for (int i = 0; i < x.size(); i++)
        c[i] = real * x[i];

    return c;
}

std::vector<double> operator*(SparseMatrix &A, const std::vector<double> &x)
{
    std::vector<double> y(A.rows());

    for (int i = 0; i < A.rows(); i++)
        for (int k = A._AI[i]; k < A._AI[i + 1]; k++)
        {
            y[i] += A._AA[k] * x[A._AJ[k]];
        }
    return y;
}

std::vector<double> prodmat(SparseMatrix &A, const std::vector<double> &x, int ibeg)
{
    std::vector<double> y(A.rows());

    for (int i = 0; i < A.rows(); i++)
        for (int k = A._AI[i]; k < A._AI[i + 1]; k++)
        {
            y[i] += A._AA[k] * x[A._AJ[k] - ibeg];
        }
    return y;
}

void save_sol_gnuplot(std::vector<double> &u, double t, Data_File *df, int me, int ibeg, int iend)
{
    int Nx(df->get_Nx()), Ny(df->get_Ny());
    double Lx(df->get_Lx()), Ly(df->get_Ly());
    double dx(Lx / (Nx + 1)), dy(Ly / (Ny + 1));

    std::stringstream stream;
    stream << std::fixed << std::setprecision(2) << t;
    std::string s = stream.str();

    std::ofstream _file_u;
    _file_u.open("Results/sol_" + s + "_" + std::to_string(me) + ".txt", std::ios::out);

    int i, j;
    for (int k = ibeg; k <= iend; k++)
    {
        i = k / Ny;
        j = k - i * Ny;
        //std::cout<< i<< " " << j << " " << k-ibeg<< " " << u.size() << std::endl;
        _file_u << i * dx << " " << j * dy << " " << u[k - ibeg] << " " << f(i * dx, j * dy, t, df) << std::endl;
    }
}

//Echange les morceau de vecteurs manquants avec les proc voisins
//Si me pair , envoie dabord, si me impair envoie second
void mpi_share_vec(const std::vector<double> &x, std::vector<double> &x_prod, int me, int Np, int ibeg, int iend, int Ny)
{
    int size_charge = iend - ibeg + 1;
    if (me == 0)
    {
        //il manque les Ny valeurs suivantes de u pour completer le produit mat vec
        for (int i = Ny; i < size_charge + Ny; i++)
            x_prod[i] = x[i - Ny]; //Copie de la partie connue par me

        //Envoi du morceau de me pour me +1
        MPI_Send(&x[size_charge - Ny], Ny, MPI_DOUBLE, 1, 100, MPI_COMM_WORLD);

        //Copie de la partie connue par me+1
        MPI_Status Status;
        MPI_Recv(&x_prod[Ny + size_charge], Ny, MPI_DOUBLE, 1, 100, MPI_COMM_WORLD, &Status);
    }
    if (me == Np - 1)
    {
        //il manque les Ny valeurs precedentes de u pour completer le produit mat vec
        for (int i = Ny; i < size_charge + Ny; i++)
            x_prod[i] = x[i - Ny]; //Copie de la partie connue par me

        if (me % 2 == 0)
        {
            //Envoi du morceau de me pour me - 1
            MPI_Send(&x[0], Ny, MPI_DOUBLE, Np - 2, 100, MPI_COMM_WORLD);

            //Copie de la partie connue par me-1
            MPI_Status Status;
            MPI_Recv(&x_prod[0], Ny, MPI_DOUBLE, Np - 2, 100, MPI_COMM_WORLD, &Status);
        }
        else
        {
            //Copie de la partie connue par me-1
            MPI_Status Status;
            MPI_Recv(&x_prod[0], Ny, MPI_DOUBLE, Np - 2, 100, MPI_COMM_WORLD, &Status);

            //Envoi du morceau de me pour me - 1
            MPI_Send(&x[0], Ny, MPI_DOUBLE, Np - 2, 100, MPI_COMM_WORLD);
        }
    }
    if ((me > 0) and (me < Np - 1))
    {
        //il manque les Ny valeurs precedentes et suivantes de u pour completer le produit mat vec
        for (int i = Ny; i < size_charge + Ny; i++)
            x_prod[i] = x[i - Ny]; //Copie de la partie connue par me

        if (me % 2 == 0)
        {
            //Envoi du morceau de me pour me - 1 et me +1
            MPI_Send(&x[0], Ny, MPI_DOUBLE, me - 1, 100, MPI_COMM_WORLD);
            MPI_Send(&x[size_charge - Ny], Ny, MPI_DOUBLE, me + 1, 100, MPI_COMM_WORLD);

            //Copie de la partie connue par me-1 et me +1
            MPI_Status Status;
            MPI_Recv(&x_prod[0], Ny, MPI_DOUBLE, me - 1, 100, MPI_COMM_WORLD, &Status);
            MPI_Recv(&x_prod[Ny + size_charge], Ny, MPI_DOUBLE, me + 1, 100, MPI_COMM_WORLD, &Status);
        }
        else
        {
            //Copie de la partie connue par me-1 et me +1
            MPI_Status Status;
            MPI_Recv(&x_prod[0], Ny, MPI_DOUBLE, me - 1, 100, MPI_COMM_WORLD, &Status);
            MPI_Recv(&x_prod[Ny + size_charge], Ny, MPI_DOUBLE, me + 1, 100, MPI_COMM_WORLD, &Status);

            //Envoi du morceau de me pour me - 1 et me +1
            MPI_Send(&x[0], Ny, MPI_DOUBLE, me - 1, 100, MPI_COMM_WORLD);
            MPI_Send(&x[size_charge - Ny], Ny, MPI_DOUBLE, me + 1, 100, MPI_COMM_WORLD);
        }
    }
}

void gc(SparseMatrix &A, std::vector<double> &u, std::vector<double> b, double t, int me, int Np, int ibeg, int iend, Data_File *df)
{
    int kmax(df->get_kmax()), Ny(df->get_Ny());
    double epsilon(df->get_epsilon());

    int size_charge = iend - ibeg + 1; //Nombre de lignes dans la charge

    std::vector<double> rk, rk1, p, z; //Vecteurs de taille iend - ibeg + 1

    //les vecteurs utilisés pour les produits matriciels sont plus long de 2*Ny
    std::vector<double> u_prod(Ny + size_charge + Ny, -12.12); //-12.12 pour debug
    std::vector<double> p_prod(Ny + size_charge + Ny, -12.12); //-12.12 pour debug

    //Si //
    if (Np > 1)
    {
        //Echange les morceau de vecteurs manquants avec les proc voisins
        mpi_share_vec(u, u_prod, me, Np, ibeg, iend, Ny);
        rk = b - prodmat(A, u_prod, ibeg - Ny);
    }
    //Si seq
    if (Np == 1)
        rk = b - A * u;

    double beta = norme(rk);
    p = rk;

    double alpha = 0;
    double gamma = 0;
    double ps_rk = 0;

    int k = 0;
    while ((beta > epsilon) and (k < kmax))
    {
        //Si //
        if (Np > 1)
        {
        //Echange les morceau de vecteurs manquants avec les proc voisins
            mpi_share_vec(p, p_prod, me, Np, ibeg, iend, Ny);
            z = prodmat(A, p_prod, ibeg - Ny);
        }
        //Si Seq
        if (Np == 1)
            z = A * p;
        ps_rk = ps(rk, rk); //calculé une seule fois
        alpha = ps_rk / ps(z, p);
        u = u + alpha * p;
        rk1 = rk - alpha * z;
        gamma = ps(rk1, rk1) / ps_rk;
        p = rk1 + gamma * p;
        beta = sqrt(ps_rk); //Norme(rk) = sqrt(ps(rk,rk))
        rk = rk1;
        k++;
    }
    if ((me == 0) and (df->get_affichage()))
    {
        std::cout << std::endl
                  << std::endl
                  << "Iteration t=" << t << std::endl;
        std::cout << "k : " << k << std::endl;
        std::cout << "residu : " << beta << std::endl;
        std::cout << me << " " << t << " " << alpha << " " << beta << " " << gamma << std::endl;
        if (k >= kmax)
        {
            std::cout << "me " << me << "tolérence non atteinte" << std::endl;
        }
    }
}
