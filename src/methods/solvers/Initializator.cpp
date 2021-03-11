//
// Created by victor on 11.03.2020.
//

#include "Initializator.h"

void Initializator::init(Grid* g, int matr_dim, int block_dimx, int& n, int& nnz) {
    try {
        int block_num = (int)(matr_dim / block_dimx); // blocks per matrix dimension
        n = g->cCount * block_num;
        nnz = 0;
        for (int i = 0; i < g->cCount; i++){ 
            Cell& c = g->cells[i];
            nnz+=block_num * block_num;
            for (int j = 0; j < c.fCount; j++) {
                if (c.neigh[j] > -1) {
                    nnz += block_num * block_num;
                }
            }
        }
    }
    catch (Exception e){
        printf(e.getMessage());
    }
}

void Initializator::initMatrixParams(Grid* g, int matr_dim, int block_dimx, int* row_ptrs, int* col_indices, row* rows) {
    try {
        int block_num = (int)(matr_dim / block_dimx); // blocks per matrix dimension
        int block_count = block_num * block_num;
        int n = g->cCount * block_num;
        int nnz = 0;
        for (int i = 0; i < g->cCount; i++) {
            Cell& c = g->cells[i];
            nnz += block_count;
            rows[i].push_back(i);
            for (int j = 0; j < c.fCount; j++) {
                if (c.neigh[j] > -1) {
                    nnz += block_count;
                    rows[i].push_back(c.neigh[j]);
                }
            }
            std::sort(rows[i].begin(), rows[i].end());
        }
        int ind = 0;
        for (int i = 0; i < g->cCount; i++) {
            for (int j = 0; j < block_num; j++) {
                row_ptrs[i * block_num + j] = ind;
                ind += rows[i].size() * block_num;
            }
        }
        row_ptrs[n] = ind;

        ind = 0;
        for (int i = 0; i < g->cCount; i++) {
            for (int j = 0; j < block_num; j++) {
                for (int k = 0; k < rows[i].size(); k++) {
                    for (int l = 0; l < block_num; l++) {
                        col_indices[ind] = rows[i][k] * block_num + l;
                        ind++;
                    }
                }
            }
        }
    }
    catch (Exception e) {
        printf(e.getMessage());
    }
}