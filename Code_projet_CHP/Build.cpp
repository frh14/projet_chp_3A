#include "Build.h"

//Construction du second membre
std::vector<double> build_F(const double t, int ibeg, int iend, Data_File *df)
{
    std::vector<double> F; //second membre
    //Données du problème
    int Nx(df->get_Nx()), Ny(df->get_Ny());
    double Lx(df->get_Lx()), Ly(df->get_Ly());
    double dx(Lx / (Nx + 1)), dy(Ly / (Ny + 1)), dt(df->get_dt());
    double D(df->get_D());

    F.resize(0);   //taille nulle avant pushback
    int count = 0; //compteur de ligne
    for (int i = 0; i < Nx; i++)
        for (int j = 0; j < Ny; j++)
        {
            if ((count >= ibeg) and (count <= iend)) //si dans les lignes de la charge
            {
                double val = 0;
                if ((i == 0) or (i == Nx - 1))
                    val += (dt * D) / (dx * dx) * g(i * dx, j * dy, t, df);
                if ((j == 0) or (j == Ny - 1))
                    val += (dt * D) / (dy * dy) * h(i * dx, j * dy, t, df);
                if ((i > 0) and (j > 0) and (i < Nx - 1) and (j < Ny - 1))
                    val += dt * f(i * dx, j * dy, t, df);
                F.push_back(val);
            }
            count++;
        }

    return F;
}

//Classe SparseMatrix
SparseMatrix::SparseMatrix()
{
}

SparseMatrix::~SparseMatrix()
{
}

//Fonctions pour seq

void SparseMatrix::set(int &i_creux, double v, int i, int j)
{
    _AA[i_creux] = v;
    _AJ[i_creux] = j;
    _AI[i + 1] += 1;
    i_creux++;
}

void SparseMatrix::resize(int size, int rows, int cols)
{
    _AA.resize(size);
    _AI.resize(rows + 1);
    _AJ.resize(size);

    _n_nnuls = size;
    _rows = rows;
    _cols = cols;
}

void SparseMatrix::print()
{
    //Des zeros partout dans A
    double A[_rows][_cols];
    for (int i = 0; i < _rows; i++)
        for (int j = 0; j < _cols; j++)
            A[i][j] = 0;

    //Sauf dans les valeurs non nulles
    int ligne = 0;
    for (int k = 0; k < size(); k++)
    {
        A[ligne][_AJ[k]] = _AA[k];
        if (k + 1 == _AI[ligne + 1])
            ligne++;
    }

    //Affichage en matrice
    for (int i = 0; i < _rows; i++)
    {
        for (int j = 0; j < _cols; j++)
            std::cout << A[i][j] << " ";
        std::cout << std::endl;
    }
}

void SparseMatrix::print_csr(int me)
{
    std::cout << me << " AA [ ";
    for (int k = 0; k < _AA.size(); k++)
        std::cout << _AA[k] << " ";
    std::cout << " ]" << std::endl;

    std::cout << me << " AJ [ ";
    for (int k = 0; k < _AJ.size(); k++)
        std::cout << _AJ[k] << " ";
    std::cout << " ]" << std::endl;

    std::cout << me << " tableau ligne CSR [ ";
    for (int k = 0; k < _AI.size(); k++)
        std::cout << _AI[k] << " ";
    std::cout << " ]" << std::endl;
}

// ========================================
//Fonctions pour para
//La matrice doit etre push back dans l'ordre en lignes
void SparseMatrix::push(double v, int i, int j)
{
    _AA.push_back(v);
    _AJ.push_back(j);
    int n_lignes = _AI.size() - 1;
    //Si la l'indice ligne est egal au dernier indice de la matrice
    //on est tjr sur la meme ligne donc +=1
    if (i == n_lignes - 1)
        _AI[n_lignes]++;
    else
    {
        // ! Si une ligne est vide , ne marche pas 
        _AI.push_back(_AI[_AI.size() - 1] + 1);
    }
    _n_nnuls++;
}

void SparseMatrix::reserve(int size)
{
    _AA.reserve(size);
    _AI.reserve(size);
    _AJ.reserve(size);
}

SparseMatrix build_A(int ibeg, int iend, Data_File *df)
{
    int Nx(df->get_Nx()), Ny(df->get_Ny());
    double Lx(df->get_Lx()), Ly(df->get_Ly());
    double dx(Lx / (Nx + 1)), dy(Ly / (Ny + 1));

    double D(df->get_D()), dt(df->get_dt());

    double tc = 1 + 2. * D * dt * (1 / (dx * dx) + 1 / (dy * dy)); // valeur trigdiag centre
    double ts = -D * dt / (dy * dy);                               // valeur tridiag décentrée
    double d = -D * dt / (dx * dx);                                //valeur bloc diag

    //Nombre de ligne de la matrice locale à me
    int rows = iend - ibeg + 1;
    int cols = Nx * Ny;
    SparseMatrix A; // Matrice creuse

    //Initialisation du vecteur avant les pushback
    A._AA.resize(0);
    A._AJ.resize(0);
    A._AI.resize(1);
    A._AI[0] = 0;

    // pour optimiser la mémoire , on reserve 5* le nombre de ligne( un peu trop ) puis on supprimera les cases reservées en trop
    A.reserve(rows * 5);

    //definition des lignes et colones
    A.set_cols(cols);
    A.set_rows(rows);

    //Construction de la matrice A creuse
    //ligne de bloc par ligne de bloc
    //puis ligne par ligne
    A._AI[0] = 0;
    int i_creux = 0;
    //pour chaque bloc de ligne
    for (int i_bloc = 0; i_bloc < Nx; i_bloc++)
    {
        //Pour chaque ligne dans le bloc
        for (int i_smat = 0; i_smat < Ny; i_smat++)
        {
            int i = i_bloc * Ny + i_smat; //vrai indice ligne
            if ((i >= ibeg) and (i <= iend))
                //haut de la matrice
                if (i_bloc == 0)
                {
                    //haut dy bloc
                    if (i_smat == 0)
                    {
                        A.push(tc, i - ibeg, i);
                        A.push(ts, i - ibeg, i + 1);
                        A.push(d, i - ibeg, i + Ny);
                    }
                    //bas du bloc
                    else if (i_smat == Ny - 1)
                    {
                        A.push(ts, i - ibeg, i - 1);
                        A.push(tc, i - ibeg, i);
                        A.push(d, i - ibeg, i + Ny);
                    }
                    //mileu du bloc
                    else
                    {
                        A.push(ts, i - ibeg, i - 1);
                        A.push(tc, i - ibeg, i);
                        A.push(ts, i - ibeg, i + 1);
                        A.push(d, i - ibeg, i + Ny);
                    }
                }
                else if (i_bloc == Nx - 1)
                {
                    if (i_smat == 0)
                    {
                        A.push(d, i - ibeg, i - Ny);
                        A.push(tc, i - ibeg, i);
                        A.push(ts, i - ibeg, i + 1);
                    }
                    else if (i_smat == Ny - 1)
                    {
                        A.push(d, i - ibeg, i - Ny);
                        A.push(ts, i - ibeg, i - 1);
                        A.push(tc, i - ibeg, i);
                    }
                    else
                    {
                        A.push(d, i - ibeg, i - Ny);
                        A.push(ts, i - ibeg, i - 1);
                        A.push(tc, i - ibeg, i);
                        A.push(ts, i - ibeg, i + 1);
                    }
                }
                else if ((i_bloc > 0) and (i_bloc < Nx - 1))
                {
                    if (i_smat == 0)
                    {
                        A.push(d, i - ibeg, i - Ny);
                        A.push(tc, i - ibeg, i);
                        A.push(ts, i - ibeg, i + 1);
                        A.push(d, i - ibeg, i + Ny);
                    }
                    else if (i_smat == Ny - 1)
                    {
                        A.push(d, i - ibeg, i - Ny);
                        A.push(ts, i - ibeg, i - 1);
                        A.push(tc, i - ibeg, i);
                        A.push(d, i - ibeg, i + Ny);
                    }
                    else
                    {
                        A.push(d, i - ibeg, i - Ny);
                        A.push(ts, i - ibeg, i - 1);
                        A.push(tc, i - ibeg, i);
                        A.push(ts, i - ibeg, i + 1);
                        A.push(d, i - ibeg, i + Ny);
                    }
                }
        }
    }

    //Reduction du tableau , Optimisation de la memoire
    A._AA.shrink_to_fit();
    A._AI.shrink_to_fit();
    A._AJ.shrink_to_fit();

    return A;
}
