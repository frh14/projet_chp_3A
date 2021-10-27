#include <iostream>
#include <vector>
#include <mpi.h>

#include "Data_File.h"
#include "Fonctions.h"
#include "Build.h"
#include "gc.h"

int main(int argc, char *argv[])
{
    //Lacement des processus MPI
    MPI_Init(&argc, &argv);
    std::vector<double> time;
    time.push_back(MPI_Wtime());

    int me, Np;
    MPI_Comm_rank(MPI_COMM_WORLD, &me); //Numéro du proc
    MPI_Comm_size(MPI_COMM_WORLD, &Np); // Nombre de Procs

    //Verif arguments à l'exc du programe
    //Il faut mettre un datafile en arg
    if (argc < 2)
    {
        std::cout << "Entrer le nom d'un data file en argument" << std::endl;
        abort();
    }

    //Fichier datafiles avec les parametres
    Data_File *df = new Data_File;
    df->read_parameter_file(argv[1]); //Lecture du datafile

    //Si on modifie Nx et Ny avec les arguments
    if (argc == 3)
       df->set_Nx(std::stoi(argv[2]));
    
    if (argc == 4)
    {
       df->set_Nx(std::stoi(argv[2]));
       df->set_Ny(std::stoi(argv[3]));
    }
    

    if ((me == 0) and (df->get_affichage()))
        df->print_parameters();

    int Nx(df->get_Nx()), Ny(df->get_Ny());
    double dt(df->get_dt()), t_final(df->get_tfinal());

    //Calcul de la charge
    //la charge est sur le nombre de lignes
    int ibeg, iend;
    charge(me, Nx * Ny, Np, ibeg, iend);

    //Garde fou si trop de procs
    if (iend - ibeg+1<Ny)
    {
        std::cout << "! Erreur ! : Taille trop petite par rapport au nombre de proc " << std::endl;
        std::cout << "! Nombre de proc max est : " << Nx << std::endl;
        abort();
    }

    double t = 0;

    //Création du morceau de matrice creuse de la ligne ibeg à iend
    SparseMatrix A;
    A = build_A(ibeg, iend, df);
    if (Nx * Ny <= 20)
        if ((me == 0) and (df->get_affichage()))
            A.print_csr(me);

    std::vector<double> F;                      //Second membre
    std::vector<double> u(iend - ibeg + 1, 1.); //Solution
    int kfinal = t_final / dt;
    for (int k = 0; k <= kfinal; k++)
    {
        //Construction de F
        F = build_F(t + dt, ibeg, iend, df);

        time.push_back(MPI_Wtime());

        //Résoution de Au=B avec le gc
        gc(A, u, u + F, t, me, Np, ibeg, iend, df);

        time.push_back(MPI_Wtime());
        if ((me == 0) and (df->get_affichage()))
            std::cout << "temps gc " << time[time.size() - 1] - time[time.size() - 2] << std::endl;

        time.push_back(MPI_Wtime());

        //sauve sol
        if (df->get_stockage())
            save_sol_gnuplot(u, t, df,me,ibeg,iend);

        time.push_back(MPI_Wtime());
        if ((me == 0) and (df->get_affichage()))
            std::cout << "temps ecriture " << time[time.size() - 1] - time[time.size() - 2] << std::endl;

        t += dt;
    }
    time.push_back(MPI_Wtime());
    if ((me == 0) and (df->get_affichage()))
        std::cout << "temps total " << time[time.size() - 1] - time[0] << std::endl;

    //ecriture du fichier de temps
    if ((me == 0) and (not(df->get_affichage())))
        std::cout << time[time.size() - 1] - time[0] << std::endl;

    MPI_Finalize();
    return 0;
}
