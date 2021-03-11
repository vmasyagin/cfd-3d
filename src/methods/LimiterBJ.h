#pragma once

#include "Limiter.h"
#include "grid.h"

class LimiterBJ : public Limiter
{
public:
	LimiterBJ(FEM_DG *mthd, Grid *grid, double **ro, double **ru, double **rv, double **rw, double **re, int fCount);
	virtual ~LimiterBJ();
	virtual void run();
	virtual const char* getName() { return "Barth-Jespersen"; }
private:
	void initLimiterParameters();
	//void chooseDirection(int* n, Point pm, double* alpha, int* neigh);
	//void getMatrLR(double** L, double** R, Param& par, const Point& n);
	//void calcPositivityLimiter();
private:
	Grid *grid;
	int cellsCount;

	/*double*** limAlpha;
	int*** limNeigh;
	Vector** limLm;
	Vector** limLmN;
	Point** limPm;

	double* deltaS1;
	double* deltaS2;
	double** deltaU1;
	double** deltaU2;

	double** matrL;
	double** matrR;*/

	FEM_DG *method;

	int funcCount;

	double **fRO;
	double **fRU;
	double **fRV;
	double **fRW;
	double **fRE;

	double *fROlim;
	double *fRUlim;
	double *fRVlim;
	double *fRWlim;
	double *fRElim;

	const double LIMITER_ALPHA = 1.9;

	const double GAM = 1.4;
	const double AGAM = GAM - 1.0;

	const double EPS = 1.0e-10;
};

