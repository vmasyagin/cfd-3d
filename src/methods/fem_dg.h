#pragma once

#include "method.h"

class FEM_DG : public Method
{
protected:
	const static int BASE_FUNC_COUNT = 4;
	const static int GP_CELL_COUNT = 4;
	const static int GP_FACE_COUNT = 3;

	double	TAU;

	double **ro;
	double **ru;
	double **rv;
	double **rw;
	double **re;
	double **tau_xx;
	double **tau_yy;
	double **tau_zz;
	double **tau_xy;
	double **tau_xz;
	double **tau_yz;

	Point** cellGP;
	Point** faceGP;
	double** cellGW;
	double* cellJ;
	double** faceGW;
	double* faceJ;

	double* Y;			// массив расто€ний от цента €чейки до стенки

	LimiterParam limParam;

public:
	virtual void getFields(double &fRO, double &fRU, double &fRV, double& fRW, double &fRE, int iCell, Point p) = 0;
	virtual void getFields(double &fRO, double &fRU, double &fRV, double& fRW, double &fRE, int iCell, double x, double y, double z) = 0;

	virtual void getFieldsLinear(double& fRO, double& fRU, double& fRV, double& fRW, double& fRE, int iCell, Point p) {};
	virtual void getFieldsLinear(double& fRO, double& fRU, double& fRV, double& fRW, double& fRE, int iCell, double x, double y, double z) {};

	virtual double getField(int fld, int iCell, Point p) = 0;
	virtual double getField(int fld, int iCell, double x, double y, double z) = 0;

	virtual double getFieldLinear(int fld, int iCell, Point p) { return 0.0; };
	virtual double getFieldLinear(int fld, int iCell, double x, double y, double z) { return 0.0; };

	virtual double getF(int id, int iCell, Point p) = 0;
	virtual double getF(int id, int iCell, double x, double y, double z) = 0;

	int getFaceGPCount() { return GP_FACE_COUNT; }

	Point getFaceGP(int iFace, int i) { return faceGP[iFace][i]; }

	double getFaceGW(int iFace, int i) { return faceGW[iFace][i]; }

	double getFaceJ(int iFace) { return faceJ[iFace]; }

	int getCellGPCount() { return GP_CELL_COUNT; }

	Point getCellGP(int iCell, int i) { return cellGP[iCell][i]; }

	double getCellGW(int iCell, int i) { return cellGW[iCell][i]; }

	double getCellJ(int iCell) { return cellJ[iCell]; }

	virtual Material& getMaterial(int iCell) = 0;

	virtual void calcFlux(double& fr, double& u, double& v, double& w, double& fe, Param pL, Param pR, Vector n, double GAM) {};

	void calcFluxFields(Param& pF, Param pL, Param pR, Vector n, double GAM) {};

	virtual void boundaryCond(int iFace, Param& pL, Param& pR) {};

	double getTau() { return TAU; }

	virtual void calcRoeAverage(Param& average, Param pL, Param pR, double GAM, Vector n) = 0;

	friend class Limiter;
};