#pragma once

#include "Limiter.h"
#include "grid.h"

class LimiterEntropy : public Limiter
{
public:
	LimiterEntropy(FEM_DG *mthd, Grid *grid, double **ro, double **ru, double **rv, double **rw, double **re, int fCount, LimiterParam* par);
	virtual ~LimiterEntropy();
	virtual void run();
	virtual const char* getName() { return "Entropy"; }
private:
	void initLimiterParameters();
	double* getField(int id, int iCell);
	double calcPsi(int iCell, double lambda, double mu);
	double calcMinPhi(int iCell, double lambda);
	double getDsDu(int iCell, int id, Point pt);
	double getDpDu(int iCell, int id, Point pt);
	double calcDpsiDmu(int iCell, double lambda);
	double calcMinDphiDlambda(int iCell);
	void calcCs();
private:
	Grid *grid;
	int cellsCount;

	double PRef;
	double RRef;

	double Pmin;

	double* Cs;

	FEM_DG *method;

	int funcCount;

	double **fRO;
	double **fRU;
	double **fRV;
	double **fRW;
	double **fRE;

	const double GAM = 1.4;
	const double AGAM = GAM - 1.0;

	const double EPS = 1.0e-10;
};

