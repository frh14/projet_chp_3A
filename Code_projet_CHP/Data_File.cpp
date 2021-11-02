#include "Data_File.h"

using namespace std;

Data_File::Data_File()
{
}


void Data_File::read_parameter_file(string parameter_file_name)
{
    _parameter_file = parameter_file_name;
    ifstream dataFile(parameter_file_name.data());
    if (!dataFile.is_open())
    {
        cout << "Unable to open parameter file " << parameter_file_name << endl;
        abort();
    }

    string line;

    while (!dataFile.eof())
    {
        getline(dataFile, line);

        if (line.find("Nx") != std::string::npos)
        {
            getline(dataFile, line);
            _Nx = stoi(line);
        }

        if (line.find("Ny") != std::string::npos)
        {
            getline(dataFile, line);
            _Ny = stoi(line);
        }

        if (line.find("Lx") != std::string::npos)
        {
            getline(dataFile, line);
            _Lx = stod(line);
        }

        if (line.find("Ly") != std::string::npos)
        {
            getline(dataFile, line);
            _Ly = stod(line);
        }

        if (line.find("D") != std::string::npos)
        {
            getline(dataFile, line);
            _D = stod(line);
        }

        if (line.find("dt") != std::string::npos)
        {
            getline(dataFile, line);
            _dt = stod(line);
        }

        if (line.find("conditions") != std::string::npos)
        {
            getline(dataFile, line);
            _conditions = stoi(line);
        }

        if (line.find("kmax") != std::string::npos)
        {
            getline(dataFile, line);
            _kmax = stoi(line);
        }

        if (line.find("epsilon") != std::string::npos)
        {
            getline(dataFile, line);
            _epsilon = stod(line);
        }

        if (line.find("tfinal") != std::string::npos)
        {
            getline(dataFile, line);
            _tfinal = stod(line);
        }

        if (line.find("affichage") != std::string::npos)
        {
            getline(dataFile, line);
            _affichage = stoi(line);
        }

        if (line.find("stockage") != std::string::npos)
        {
            getline(dataFile, line);
            _stockage = stoi(line);
        }
    }
}

void Data_File::print_parameters()
{
    cout << "=================================================" << endl;
    cout << "Affichage des parametres : " << _parameter_file << endl;
    cout << "Nx " << _Nx << endl;
    cout << "Ny " << _Ny << endl;
    cout << "Lx " << _Lx << endl;
    cout << "Ly " << _Ly << endl;
    cout << "D " << _D << endl;
    cout << "dt " << _dt << endl;
    cout << "conditions " << _conditions << endl;
    cout << "kmax " << _kmax << endl;
    cout << "epsilon " << _epsilon << endl;
    cout << "tfinal " << _tfinal << endl;
    cout << "affichage " << _affichage << endl;
    cout << "stockage " << _stockage << endl;
    cout << "=================================================" << endl;
}