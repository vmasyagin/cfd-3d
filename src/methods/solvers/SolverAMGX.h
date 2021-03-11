#pragma once
#include "MatrixSolver.h"
#include "amgx_c.h"
#include <map>
#include <string>
#include <algorithm>
#include <mpi.h>
//#include "cuda_runtime.h"

///* CUDA error macro */
//#define CUDA_SAFE_CALL(call) do {                                 \
//  cudaError_t err = call;                                         \
//  if(cudaSuccess != err) {                                        \
//    fprintf(stderr, "Cuda error in file '%s' in line %i : %s.\n", \
//            __FILE__, __LINE__, cudaGetErrorString( err) );       \
//    exit(EXIT_FAILURE);                                           \
//  } } while (0)

typedef std::vector<int> row;

class SolverAMGX :
	public MatrixSolver
{
public:
	~SolverAMGX();
	virtual void init(Grid* g, int matrDimension, int blockDimension);

	virtual void zero();

	virtual void setMatrElement(int i, int j, double** matrDim);
	virtual void setRightElement(int i, double* vectDim);
	virtual void addMatrElement(int i, int j, double** matrDim);
	virtual void addRightElement(int i, double* vectDim);
	virtual void createMatrElement(int i, int j) {}
	virtual void printToFile(const char* fileName);

protected:
	void initMatrVectors();

protected:
	AMGX_Mode mode = AMGX_mode_dDDI; // using dDDI by default
    AMGX_config_handle config;
    AMGX_resources_handle rsrc;
    AMGX_solver_handle solver;
    AMGX_matrix_handle A;
    AMGX_vector_handle rhs;
    AMGX_vector_handle sln;
	AMGX_SOLVE_STATUS status;

	////MPI (with CUDA GPUs)
	//int rank = 0;
	//int lrank = 0;
	//int nranks = 0;
	//int nx, ny;
	//int gpu_count = 0;
	//MPI_Comm amgx_mpi_comm = MPI_COMM_WORLD;

    int n, nnz, block_dimx, block_dimy, block_size, num_neighbors, matr_dim, block_num;

    int *row_ptrs = NULL, *col_indices = NULL, *neighbors = NULL;
    double *values = NULL, *diag = NULL, *dh_x = NULL, *dh_b = NULL;
    int *h_row_ptrs = NULL, *h_col_indices = NULL;
    double *h_values = NULL, *h_diag = NULL, *h_x = NULL, *h_b = NULL;

	row* rows;

	char* solver_name;

	int PRINT_LEVEL = 0;
};

