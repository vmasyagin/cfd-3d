//
// Created by victor on 11.03.2020.
//
#include <grid.h>
#include <algorithm>

#ifndef CFD_2D_INITIALIZATOR_H
#define CFD_2D_INITIALIZATOR_H

#endif //CFD_2D_INITIALIZATOR_H

//typedef std::map<int, std::vector<double>> block_row;
typedef std::vector<int> row;

class Initializator{
public:
    static void init(Grid* g, int matr_dim, int block_dimx, int& n, int& nnz);
    static void initMatrixParams(Grid* g, int matr_dim, int block_dimx, int* h_row_ptrs, int* h_col_indices, row* rows);
};
