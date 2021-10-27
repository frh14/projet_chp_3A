#pragma once

#include <iostream>
#include <vector>
#include "Fonctions.h"
#include "Data_File.h"

//Au format csr
class SparseMatrix
{
public:
    std::vector<double> _AA;
    std::vector<int> _AI;
    std::vector<int> _AJ;

private:
    int _n_nnuls;
    int _rows;
    int _cols;

public:
    SparseMatrix();
    ~SparseMatrix();

    int rows() { return _rows; };
    int cols() { return _cols; };

    //Utile en seq
    void set(int &i_creux, double v, int i, int j);
    void resize(int size, int rows, int cols);
    void print();
    void print_csr(int me);

    //Utile en parallele
    int size() { return _AJ.size(); };
    void set_rows(int rows) { _rows = rows; };
    void set_cols(int cols) { _cols = cols; };
    void push(double v, int i, int j);
    void reserve(int size);
};

//Build le second membre
std::vector<double> build_F(const double t, int ibeg, int iend, Data_File *df);
SparseMatrix build_A(int ibeg, int iend, Data_File *df);