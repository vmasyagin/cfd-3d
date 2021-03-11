#include "SolverAMGX.h"
#include "config.h"


void SolverAMGX::initMatrVectors()
{
    AMGX_initialize();
    AMGX_initialize_plugins();
    std::string config_name = CONFIGS_DIR;
    config_name.append(solver_name + 5);
    config_name += ".json";
    AMGX_config_create_from_file(&config, config_name.c_str());
    AMGX_resources_create_simple(&rsrc, config);

    AMGX_solver_create(&solver, rsrc, mode, config);
    AMGX_matrix_create(&A, rsrc, mode);
    AMGX_vector_create(&rhs, rsrc, mode);
    AMGX_vector_create(&sln, rsrc, mode);
}

void SolverAMGX::zero() {
    AMGX_matrix_destroy(A);
    AMGX_vector_destroy(rhs);
    AMGX_vector_destroy(sln);

    AMGX_matrix_create(&A, rsrc, mode);
    AMGX_vector_create(&rhs, rsrc, mode);
    AMGX_vector_create(&sln, rsrc, mode);

    memset(h_values, 0, sizeof(double)*nnz*block_size);
    memset(h_b, 0, sizeof(double)*n*block_dimx);
    memset(h_x, 0, sizeof(double)*n*block_dimx);
    //memset(h_col_indices, 0, nnz * sizeof(int));
    //memset(h_row_ptrs, 0, (n + 1)*sizeof(int));
}

SolverAMGX::~SolverAMGX()
{
    AMGX_matrix_destroy(A);
    AMGX_vector_destroy(rhs);
    AMGX_vector_destroy(sln);
    AMGX_solver_destroy(solver);
    AMGX_resources_destroy(rsrc);

    AMGX_finalize_plugins();
    AMGX_finalize();
	delete[] x;
}

void SolverAMGX::setMatrElement(int i, int j, double** matrDim)
{
    // TODO: реализовать позже при необходимости
}

void SolverAMGX::addMatrElement(int i, int j, double** matrDim)
{
    std::vector<double> row;
    row.resize(block_size, 0);
    auto it = std::find(rows[i].begin(), rows[i].end(), j);
    int pos = std::distance(rows[i].begin(), it);
    for (int i_block = 0; i_block < block_num; i_block++){
        for (int j_block = 0; j_block < block_num; j_block++){
            int index = h_row_ptrs[i * block_num] + pos * block_num + i_block * rows[i].size() * block_num + j_block; // current block id
            memcpy(&row[0], &h_values[index * block_size], block_size * sizeof(double));
            for (int ii = 0; ii < block_dimx; ii++) {
                for (int jj = 0; jj < block_dimy; jj++) {       
                    row[ii * block_dimx + jj] += matrDim[i_block * block_dimx + ii][j_block * block_dimy + jj];
                }
            }
            memcpy(&h_values[index * block_size], &row[0], block_size * sizeof(double));
        }
    }
}


void SolverAMGX::setRightElement(int i, double* vectDim)
{
    memcpy(&h_b[i*matr_dim], vectDim, matr_dim * sizeof(double));
}


void SolverAMGX::addRightElement(int i, double* vectDim)
{
    for (int ii = 0; ii < matr_dim; ii++) {
        h_b[ii + i * matr_dim] += vectDim[ii];
    }
}

void SolverAMGX::init(Grid* g, int matrDimension, int blockDimension) {

    ///* MPI init (with CUDA GPUs) */
    ////MPI
    //MPI_Comm_size(amgx_mpi_comm, &nranks);
    //MPI_Comm_rank(amgx_mpi_comm, &rank);

    ////CUDA GPUs
    //CUDA_SAFE_CALL(cudaGetDeviceCount(&gpu_count));
    //lrank = rank % gpu_count;
    //CUDA_SAFE_CALL(cudaSetDevice(lrank));
    //printf("Process %d selecting device %d\n", rank, lrank);

    Initializator::init(g, matrDimension, blockDimension, n, nnz);
    block_dimx = blockDimension;
    block_dimy = blockDimension;
    block_size = block_dimx * block_dimy;
    matr_dim = matrDimension;
    block_num = (int)matr_dim / block_dimx;

    h_col_indices = new int[nnz];
    h_row_ptrs = new int[n + 1];
    rows = new row[g->cCount];
    Initializator::initMatrixParams(g, matrDimension, blockDimension, h_row_ptrs, h_col_indices, rows);

    h_values = new double[nnz * block_size];
    h_b = new double[n * block_dimx];
    h_x = new double[n * block_dimx];

    x = new double[n * block_dimx];

    initMatrVectors();
}

void SolverAMGX::printToFile(const char* fileName)
{
    FILE * fp = fopen(fileName, "w");
	for (int i = 0; i < n*block_dimx; i++) {
	    fprintf(fp, "%25.16e  ", values[h_row_ptrs[i]]);
    }
	fprintf(fp, "\n\n=============================================================================================\n\n\n");
	for (int i = 0; i < n*block_dimx; i++) {
		fprintf(fp, "%25.16e  ", x[i]);
	}

	fclose(fp);
}





