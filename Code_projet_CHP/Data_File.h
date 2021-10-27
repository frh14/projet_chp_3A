#pragma once

#include <string>
#include <fstream>
#include <iostream>

class Data_File
{
private:
    //Fichiers
    std::string _parameter_file;

    //Parametres lu
    int _Nx;
    int _Ny;
    double _Lx;
    double _Ly;
    double _D;
    double _dt;
    int _conditions;
    int _kmax;
    double _epsilon;
    double _tfinal;
    bool _affichage;
    bool _stockage;

public:
    Data_File();

    void read_parameter_file(std::string parameter_file_name);
    void print_parameters();


    void set_Nx(int Nx) {_Nx = Nx;};
    void set_Ny(int Ny) {_Ny = Ny;};

    int get_Nx() { return _Nx; };
    int get_Ny() { return _Ny; };
    double get_Lx() { return _Lx; };
    double get_Ly() { return _Ly; };
    double get_D() { return _D; };
    double get_dt() { return _dt; };
    int get_conditions() { return _conditions; };
    int get_kmax() { return _kmax; };
    double get_epsilon() { return _epsilon; };
    double get_tfinal() { return _tfinal; };
    bool get_stockage() { return _stockage; };
    bool get_affichage() { return _affichage; };
};