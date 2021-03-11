#include "SolverAMGXImpl.h"
#include "global.h"
#include <cmath>
#include <math.h>

int SolverAMGXImpl::solve(double eps, int& maxIter)
{
	int result = MatrixSolver::RESULT_OK;

	/* AMG */
	{
		double final_res_norm;
        AMGX_SAFE_CALL(AMGX_pin_memory(h_x, n * block_dimx * sizeof(double)));
        AMGX_SAFE_CALL(AMGX_pin_memory(h_b, n * block_dimx * sizeof(double)));
        AMGX_SAFE_CALL(AMGX_pin_memory(h_col_indices, nnz * sizeof(int)));
        AMGX_SAFE_CALL(AMGX_pin_memory(h_row_ptrs, (n + 1)*sizeof(int)));
        AMGX_SAFE_CALL(AMGX_pin_memory(h_values, nnz * block_size * sizeof(double)));

        if (h_diag != NULL)
        {
            AMGX_SAFE_CALL(AMGX_pin_memory(h_diag, n * block_size * sizeof(double)));
        }
        /* set pointers to point to CPU (host) memory */
        row_ptrs = h_row_ptrs;
        col_indices = h_col_indices;
        values = h_values;
        diag = h_diag;
        dh_x = h_x;
        dh_b = h_b;

        //int nrings; //=1; //=2;
        //AMGX_config_get_default_number_of_rings(config, &nrings);
        //printf("nrings=%d\n",nrings);

        /* set the connectivity information (for the vector) */
        AMGX_SAFE_CALL(AMGX_vector_bind(sln, A));
        AMGX_SAFE_CALL(AMGX_vector_bind(rhs, A));
        /* upload the matrix (and the connectivity information) */
        AMGX_SAFE_CALL(AMGX_matrix_upload_all(A, n, nnz, block_dimx, block_dimy, row_ptrs, col_indices, values, diag));
        /* upload the vector (and the connectivity information) */
        AMGX_SAFE_CALL(AMGX_vector_upload(sln, n, block_dimx, dh_x));
        AMGX_SAFE_CALL(AMGX_vector_upload(rhs, n, block_dimx, dh_b));

        AMGX_SAFE_CALL(AMGX_solver_setup(solver, A));

        AMGX_SAFE_CALL(AMGX_solver_solve(solver, rhs, sln));

        AMGX_SAFE_CALL(AMGX_unpin_memory(h_x));
        AMGX_SAFE_CALL(AMGX_unpin_memory(h_b));
        AMGX_SAFE_CALL(AMGX_unpin_memory(h_col_indices));
        AMGX_SAFE_CALL(AMGX_unpin_memory(h_row_ptrs));
        AMGX_SAFE_CALL(AMGX_unpin_memory(h_values));
		/* Run info - needed logging turned on */
		int initMaxIter = maxIter;
        AMGX_SAFE_CALL(AMGX_solver_get_iterations_number(solver, &maxIter));
        AMGX_SAFE_CALL(AMGX_solver_get_iteration_residual(solver, maxIter-1, 0, &final_res_norm));
		if (initMaxIter <= maxIter) {
			result |= MatrixSolver::RESULT_ERR_MAX_ITER;
		}
		if (/*final_res_norm >= eps || */!std::isfinite(final_res_norm)) {
			result |= MatrixSolver::RESULT_ERR_CONVERG;
		}
	}

	if (result == MatrixSolver::RESULT_OK) {
		/* get the local solution */
        AMGX_SAFE_CALL(AMGX_pin_memory(x, n * block_dimx * sizeof(double)));
        AMGX_SAFE_CALL(AMGX_vector_download(sln, x));
        AMGX_SAFE_CALL(AMGX_unpin_memory(x));
        AMGX_SAFE_CALL(AMGX_vector_upload(sln, n, block_dimx, dh_x));
        AMGX_SAFE_CALL(AMGX_vector_upload(rhs, n, block_dimx, dh_b));
	}


	return result;
}


