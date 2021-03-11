#ifndef _MatrixSolver_
#define _MatrixSolver_

#include "CSR.h"
#include "_hypre_utilities.h"
#include "HYPRE_krylov.h"
#include "HYPRE.h"
#include "HYPRE_parcsr_ls.h"
#include "grid.h"
#include "Initializator.h"

class MatrixSolver
{
public:

	static const int RESULT_OK					= 0x0000;
	static const int RESULT_ERR_ZERO_DIAG		= 0x0001;
	static const int RESULT_ERR_MAX_ITER		= 0x0002;
	static const int RESULT_ERR_CONVERG			= 0x0004;

	static MatrixSolver* create(const char* solverName);
	
	
	virtual ~MatrixSolver();

	virtual void init(int cellsCount, int blockDimension);
    virtual void init(Grid* g, int matrDimension, int blockDimension) {init(g->cCount, matrDimension);};
	
	virtual void zero();
	
	virtual void setMatrElement(int i, int j, double** matrDim);
	virtual void setRightElement(int i, double* vectDim);
	virtual void addMatrElement(int i, int j, double** matrDim);
	virtual void addRightElement(int i, double* vectDim);
	virtual void createMatrElement(int i, int j);
	virtual void initCSR();

	virtual int solve(double eps, int& maxIter) = 0;
	virtual char* getName() = 0;

	virtual void setParameter(const char* name, int val) {}
	virtual void setParameter(const char* name, double val) {}
	virtual void setParameter(const char* name, const char* val) {}

	virtual void printToFile(const char* fileName);

	CSRMatrix	*a;
	HYPRE_Int    blockDim;
	double		*b;
	double		*x;
	Initializator *initializator;
};

class SolverJacobi: public MatrixSolver
{
public:
	SolverJacobi() {
		tempXAlloc = false;
	}
	~SolverJacobi() {
		if (tempXAlloc) {
			delete [] tempX;
			tempXAlloc = false;
		}
	}
	virtual int	solve(double eps, int& maxIter);
	double			*tempX;
	bool			tempXAlloc;
};




#endif