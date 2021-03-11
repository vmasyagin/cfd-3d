#include "LimiterCockburn.h"

#define MIN(a,b) ((a)<(b)?(a):(b))
#define MAX(a,b) ((a)>(b)?(a):(b))

double _det_(double a11, double a12, double a13,
	double a21, double a22, double a23,
	double a31, double a32, double a33)
{
	return a11 * a22 * a33 + a12 * a23 * a31 + a21 * a32 * a13 - a31 * a22 * a13 - a21 * a12 * a33 - a11 * a32 * a23;
}

double _det_(double a11, double a12, double a13, double a14,
	double a21, double a22, double a23, double a24,
	double a31, double a32, double a33, double a34,
	double a41, double a42, double a43, double a44)
{
	double d1 = _det_(a21, a22, a23,
		a31, a32, a33,
		a41, a42, a43);
	double d2 = _det_(a11, a12, a13,
		a31, a32, a33,
		a41, a42, a43);
	double d3 = _det_(a11, a12, a13,
		a21, a22, a23,
		a41, a42, a43);
	double d4 = _det_(a11, a12, a13,
		a21, a22, a23,
		a31, a32, a33);
	return -a14 * d1 + a24 * d2 - a34 * d3 + a44 * d4;
}

bool _solve_SLE_3_(double  a11, double  a12, double  a13,
	double  a21, double  a22, double  a23,
	double  a31, double  a32, double  a33,
	double  b1, double  b2, double  b3,
	double& x1, double& x2, double& x3)
{
	double det = _det_(a11, a12, a13,
		a21, a22, a23,
		a31, a32, a33);
	if (fabs(det) <= 1.e-12)
	{
		log("WARNING: matrix determinat is nearly zero,\n         when choosing direction for Cockburn limiter!\n");
		x1 = 0.0;
		x2 = 0.0;
		x3 = 0.0;
		return false;
	}
	x1 = _det_(b1, a12, a13,
		b2, a22, a23,
		b3, a32, a33) / det;
	x2 = _det_(a11, b1, a13,
		a21, b2, a23,
		a31, b3, a33) / det;
	x3 = _det_(a11, a12, b1,
		a21, a22, b2,
		a31, a32, b3) / det;
	return true;
}

void _solve_SLE_4_(double  a11, double  a12, double  a13, double a14,
	double  a21, double  a22, double  a23, double a24,
	double  a31, double  a32, double  a33, double a34,
	double  a41, double  a42, double  a43, double a44,
	double  b1, double  b2, double  b3, double b4,
	double& x1, double& x2, double& x3, double& x4)
{
	double det = _det_(a11, a12, a13, a14,
		a21, a22, a23, a24,
		a31, a32, a33, a34,
		a41, a42, a43, a44);

	x1 = _det_(b1, a12, a13, a14,
		b2, a22, a23, a24,
		b3, a32, a33, a34,
		b4, a42, a43, a44) / det;

	x2 = _det_(a11, b1, a13, a14,
		a21, b2, a23, a24,
		a31, b3, a33, a34,
		a41, b4, a43, a44) / det;

	x3 = _det_(a11, a12, b1, a14,
		a21, a22, b2, a24,
		a31, a32, b3, a34,
		a41, a42, b4, a44) / det;

	x4 = _det_(a11, a12, a13, b1,
		a21, a22, a23, b2,
		a31, a32, a33, b3,
		a41, a42, a43, b4) / det;

}

LimiterCockburn::LimiterCockburn(FEM_DG *mthd, Grid *grid, double **ro, double **ru, double **rv, double **rw, double **re, int fCount)
{
	this->method = mthd;
	this->grid = grid;
	cellsCount = grid->cCount;
	funcCount = fCount;
	initLimiterParameters();

	fRO = ro;
	fRU = ru;
	fRV = rv;
	fRW = rw;
	fRE = re;
}

LimiterCockburn::~LimiterCockburn()
{
	delete[] grid;

	for (int i = 0; i < cellsCount; i++)
	{
		for (int j = 0; j < 4; j++) {
			delete[] limAlpha[i][j];
			delete[] limNeigh[i][j];
		}
		delete[] limAlpha[i];
		delete[] limNeigh[i];

		delete[] limLm[i];
		delete[] limLmN[i];
		delete[] limPm[i];

		delete[] fROlim[i];
		delete[] fRUlim[i];
		delete[] fRVlim[i];
		delete[] fRWlim[i];
		delete[] fRElim[i];

		delete[] fRO[i];
		delete[] fRU[i];
		delete[] fRV[i];
		delete[] fRW[i];
		delete[] fRE[i];

	}

	delete[] fRO;
	delete[] fRU;
	delete[] fRV;
	delete[] fRW;
	delete[] fRE;

	delete[] fROlim;
	delete[] fRUlim;
	delete[] fRVlim;
	delete[] fRWlim;
	delete[] fRElim;

	delete[] limAlpha;
	delete[] limNeigh;
	delete[] limLm;
	delete[] limLmN;
	delete[] limPm;

	delete[] deltaS1;
	delete[] deltaS2;

	for (int m = 0; m < 4; m++)
	{
		delete[] deltaU1[m];
		delete[] deltaU2[m];
	}
	delete[] deltaU1;
	delete[] deltaU2;

	for (int i = 0; i < 5; i++)
	{
		delete[] matrL[i];
		delete[] matrR[i];
	}
	delete[] matrL;
	delete[] matrR;

	delete[] method;
}

void LimiterCockburn::initLimiterParameters()
{
	limAlpha = new double** [cellsCount];
	limNeigh = new int** [cellsCount];
	limLm = new Vector* [cellsCount];
	limLmN = new Vector* [cellsCount];
	limPm = new Point* [cellsCount];

	fROlim = new double* [cellsCount];
	fRUlim = new double* [cellsCount];
	fRVlim = new double* [cellsCount];
	fRWlim = new double* [cellsCount];
	fRElim = new double* [cellsCount];

	for (int i = 0; i < cellsCount; i++)
	{
		limAlpha[i] = new double*[4];
		limAlpha[i][0] = new double[3];
		limAlpha[i][1] = new double[3];
		limAlpha[i][2] = new double[3];
		limAlpha[i][3] = new double[3];

		limNeigh[i] = new int*[4];
		limNeigh[i][0] = new int[3];
		limNeigh[i][1] = new int[3];
		limNeigh[i][2] = new int[3];
		limNeigh[i][3] = new int[3];

		limLm[i] = new Vector[4];
		limLmN[i] = new Vector[4];
		limPm[i] = new Point[4];

		fROlim[i] = new double[funcCount];
		fRUlim[i] = new double[funcCount];
		fRVlim[i] = new double[funcCount];
		fRWlim[i] = new double[funcCount];
		fRElim[i] = new double[funcCount];

	}

	deltaS1 = new double[5];
	deltaS2 = new double[5];
	deltaU1 = new double*[4];
	deltaU2 = new double*[4];
	for (int m = 0; m < 4; m++)
	{
		deltaU1[m] = new double[5];
		deltaU2[m] = new double[5];
	}

	matrL = new double*[5];
	matrR = new double*[5];
	for (int i = 0; i < 5; i++)
	{
		matrL[i] = new double[5];
		matrR[i] = new double[5];
	}

	// находим угол и коэффициенты разложения
	for (int iCell = 0; iCell < cellsCount; iCell++)
	{
		int n[5]; // 0 - текущая ячейка, 1-4 - номера соседних ячеек
		n[0] = iCell;
		for (int i = 0; i < 4; i++)
		{
			n[i + 1] = grid->cells[iCell].neigh[i];
		}

		if (n[1] > 0 && n[2] > 0 && n[3] > 0 && n[4] > 0) // не лимитируем граничные ячейки
		{
			for (int m = 0; m < 4; m++)
			{
				limPm[iCell][m] = grid->faces[grid->cells[iCell].facesInd[m]].c;
				limLm[iCell][m] = limPm[iCell][m] - grid->cells[iCell].c;
				double tmp = sqrt(limLm[iCell][m].x * limLm[iCell][m].x + limLm[iCell][m].y * limLm[iCell][m].y + limLm[iCell][m].z * limLm[iCell][m].z);
				limLmN[iCell][m] /= tmp;
				chooseDirection(n, limPm[iCell][m], limAlpha[iCell][m], limNeigh[iCell][m]);
			}
		}
	}
}

inline double MINMOD(double a, double b)
{
	if (a*b>0.0)
	{
		if (abs(a)<abs(b))
			return a;
		else
			return b;
	}
	else {
		return 0;
	}
}

inline double MINMOD_B(double a, double b)
{
	return MINMOD(a, b);
}

void LimiterCockburn::getMatrLR(double** L, double** R, Param& par, const Point& n)
{
	double nx = n.x;
	double ny = n.y;
	double nz = n.z;

	double fU = par.u;
	double fV = par.v;
	double fW = par.w;

	double un = fU * nx + fV * ny + fW * nz;

	double cz = par.cz;
	double cz2 = cz * cz;

	double q2 = 0.5 * (fU * fU + fV * fV + fW * fW);

	double fH = par.e + q2 + par.p / par.r;

	R[0][0] = nx;		
	R[0][1] = ny;		
	R[0][2] = nz;	
	R[0][3] = 1.0;
	R[0][4] = 1.0;

	R[1][0] = fU * nx;
	R[1][1] = fU * ny - cz * nz;
	R[1][2] = fU * nz + cz * ny;
	R[1][3] = fU + cz * nx;
	R[1][4] = fU - cz * nx;

	R[2][0] = fV * nx + cz * nz;
	R[2][1] = fV * ny;
	R[2][2] = fV * nz - cz * nx;
	R[2][3] = fV + cz * ny;
	R[2][4] = fV - cz * ny;

	R[3][0] = fW * nx - cz * ny;
	R[3][1] = fW * ny + cz * nx;
	R[3][2] = fW * nz;
	R[3][3] = fW + cz * nz;
	R[3][4] = fW - cz * nz;

	R[4][0] = q2 * nx + cz * fV * nz - cz * fW * ny;
	R[4][1] = q2 * ny + cz * fW * nx - cz * fU * nz;
	R[4][2] = q2 * nz + cz * fU * ny - cz * fV * nx;
	R[4][3] = fH + cz * un;
	R[4][4] = fH - cz * un;


	L[0][0] = (nx * (cz2 - AGAM * q2) + cz * (fW * ny - fV * nz)) / cz2;
	L[0][1] = nx * AGAM * fU / cz2;
	L[0][2] = (nx * AGAM * fV + cz * nz) / cz2;
	L[0][3] = (nx * AGAM * fW - cz * ny) / cz2;
	L[0][4] = -nx * AGAM / cz2;

	L[1][0] = (ny * (cz2 - AGAM * q2) + cz * (fU * nz - fW * nx)) / cz2;
	L[1][1] = (ny * AGAM * fU - cz * nz) / cz2;
	L[1][2] = ny * AGAM * fV / cz2;
	L[1][3] = (ny * AGAM * fW + cz * nx) / cz2;
	L[1][4] = -ny * AGAM / cz2;

	L[2][0] = (nz * (cz2 - AGAM * q2) + cz * (fV * nx - fU * ny)) / cz2;
	L[2][1] = (nz * AGAM * fU + cz * ny) / cz2;
	L[2][2] = (nz * AGAM * fV - cz * nx) / cz2;
	L[2][3] = nz * AGAM * fW / cz2;
	L[2][4] = -nz * AGAM / cz2;
	
	L[3][0] = 0.5 * (AGAM * q2 - cz * un) / cz2;
	L[3][1] = 0.5 * (cz * nx - AGAM * fU) / cz2;
	L[3][2] = 0.5 * (cz * ny - AGAM * fV) / cz2;
	L[3][3] = 0.5 * (cz * nz - AGAM * fW) / cz2;
	L[3][4] = 0.5 * AGAM / cz2;

	L[4][0] = 0.5 * (AGAM * q2 + cz * un) / cz2;
	L[4][1] = 0.5 * (-cz * nx - AGAM * fU) / cz2;
	L[4][2] = 0.5 * (-cz * ny - AGAM * fV) / cz2;
	L[4][3] = 0.5 * (-cz * nz - AGAM * fW) / cz2;
	L[4][4] = 0.5 * AGAM / cz2;
}

void LimiterCockburn::run()
{
	for (int i = 0; i < cellsCount; i++)
	{
		memcpy(fROlim[i], fRO[i], funcCount * sizeof(double));
		memcpy(fRUlim[i], fRU[i], funcCount * sizeof(double));
		memcpy(fRVlim[i], fRV[i], funcCount * sizeof(double));
		memcpy(fRWlim[i], fRW[i], funcCount * sizeof(double));
		memcpy(fRElim[i], fRE[i], funcCount * sizeof(double));
	}

	for (int iCell = 0; iCell < cellsCount; iCell++) // цикл по ячейкам
	{
		if (grid->cells[iCell].flag & CELL_FLAG_BAD) continue;
		int n0 = grid->cells[iCell].neigh[0];
		int n1 = grid->cells[iCell].neigh[1];
		int n2 = grid->cells[iCell].neigh[2];
		int n3 = grid->cells[iCell].neigh[3];
		if (!(n0 > 0 && n1 > 0 && n2 > 0 && n3 > 0)) // в граничных ячейках зануляем все коэффициенты кроме нулевого
		{
			for (int j = 1; j < 4; j++) {
				fROlim[iCell][j] = 0.0;
				fRUlim[iCell][j] = 0.0;
				fRVlim[iCell][j] = 0.0;
				fRWlim[iCell][j] = 0.0;
				fRElim[iCell][j] = 0.0;
			}
			memcpy(fRO[iCell], fROlim[iCell], funcCount * sizeof(double));
			memcpy(fRU[iCell], fRUlim[iCell], funcCount * sizeof(double));
			memcpy(fRV[iCell], fRVlim[iCell], funcCount * sizeof(double));
			memcpy(fRW[iCell], fRWlim[iCell], funcCount * sizeof(double));
			memcpy(fRE[iCell], fRElim[iCell], funcCount * sizeof(double));
			continue;
		}

		double ROc, RUc, RVc, RWc, REc;
		method->getFields(ROc, RUc, RVc, RWc, REc, iCell, grid->cells[iCell].c);

		for (int m = 0; m < 4; m++) // цикл по направлениям
		{
			double	a1, a2, a3;
			int	n1, n2, n3;
			a1 = limAlpha[iCell][m][0];
			a2 = limAlpha[iCell][m][1];
			a3 = limAlpha[iCell][m][2];
			n1 = limNeigh[iCell][m][0];
			n2 = limNeigh[iCell][m][1];
			n3 = limNeigh[iCell][m][2];

			// вычисляем приращения
			double ROm, RUm, RVm, RWm, REm;
			Point pm = limPm[iCell][m];
			method->getFields(ROm, RUm, RVm, RWm, REm, iCell, pm);

			double RO1, RU1, RV1, RW1, RE1;
			Point p1 = grid->cells[n1].c;
			method->getFields(RO1, RU1, RV1, RW1, RE1, iCell, p1);

			double RO2, RU2, RV2, RW2, RE2;
			Point p2 = grid->cells[n2].c;
			method->getFields(RO2, RU2, RV2, RW2, RE2, iCell, p2);

			double RO3, RU3, RV3, RW3, RE3;
			Point p3 = grid->cells[n3].c;
			method->getFields(RO3, RU3, RV3, RW3, RE3, iCell, p3);

			deltaU1[m][0] = ROm - ROc;
			deltaU1[m][1] = RUm - RUc;
			deltaU1[m][2] = RVm - RVc;
			deltaU1[m][3] = RWm - RWc;
			deltaU1[m][4] = REm - REc;

			deltaU2[m][0] = a1 * (RO1 - ROc) + a2 * (RO2 - ROc) + a3 * (RO3 - ROc);
			deltaU2[m][1] = a1 * (RU1 - RUc) + a2 * (RU2 - RUc) + a3 * (RO3 - RUc);
			deltaU2[m][2] = a1 * (RV1 - RVc) + a2 * (RV2 - RVc) + a3 * (RV3 - RVc);
			deltaU2[m][3] = a1 * (RW1 - RWc) + a2 * (RW2 - RWc) + a3 * (RW3 - RWc);
			deltaU2[m][4] = a1 * (RE1 - REc) + a2 * (RE2 - REc) + a3 * (RE3 - REc);

			Point lm = limLmN[iCell][m];
			Param par;
			par.r = ROc;
			par.u = RUc / ROc;
			par.v = RVc / ROc;
			par.w = RWc / ROc;
			par.e = REc / ROc - 0.5 * par.U2();
			par.p = par.r * par.e * AGAM;
			par.cz = sqrt(GAM * par.p / par.r);

			// переходим к инвариантам
			getMatrLR(matrL, matrR, par, lm);
			memset(deltaS1, 0, 5 * sizeof(double));
			memset(deltaS2, 0, 5 * sizeof(double));
			for (int k = 0; k < 5; k++)
			{
				deltaS1[0] += matrL[0][k] * deltaU1[m][k];
				deltaS1[1] += matrL[1][k] * deltaU1[m][k];
				deltaS1[2] += matrL[2][k] * deltaU1[m][k];
				deltaS1[3] += matrL[3][k] * deltaU1[m][k];
				deltaS1[4] += matrL[4][k] * deltaU1[m][k];

				deltaS2[0] += matrL[0][k] * deltaU2[m][k];
				deltaS2[1] += matrL[1][k] * deltaU2[m][k];
				deltaS2[2] += matrL[2][k] * deltaU2[m][k];
				deltaS2[3] += matrL[3][k] * deltaU2[m][k];
				deltaS2[4] += matrL[4][k] * deltaU2[m][k];
			}

			// лимитируем инварианты
			deltaS1[0] = MINMOD_B(deltaS1[0], LIMITER_ALPHA * deltaS2[0]);
			deltaS1[1] = MINMOD_B(deltaS1[1], LIMITER_ALPHA * deltaS2[1]);
			deltaS1[2] = MINMOD_B(deltaS1[2], LIMITER_ALPHA * deltaS2[2]);
			deltaS1[3] = MINMOD_B(deltaS1[3], LIMITER_ALPHA * deltaS2[3]);
			deltaS1[4] = MINMOD_B(deltaS1[4], LIMITER_ALPHA * deltaS2[4]);

			// переходим к консервативным
			memset(deltaU1[m], 0, 5 * sizeof(double));
			for (int k = 0; k < 5; k++)
			{
				deltaU1[m][0] += matrR[0][k] * deltaS1[k];
				deltaU1[m][1] += matrR[1][k] * deltaS1[k];
				deltaU1[m][2] += matrR[2][k] * deltaS1[k];
				deltaU1[m][3] += matrR[3][k] * deltaS1[k];
				deltaU1[m][4] += matrR[4][k] * deltaS1[k];
			}

			// TODO: нужно ли поворачивать скорости обратно? 
		}

		double pos, neg;
		for (int k = 0; k < 4; k++) // // бежим по компонентам RO, RU, RV, RW, RE
		{

			double& d1n = deltaU1[0][k];
			double& d2n = deltaU1[1][k];
			double& d3n = deltaU1[2][k];
			double& d4n = deltaU1[3][k];
			// делаем поправки для приращений
			if (fabs(d1n + d2n + d3n + d4n) > EPS)
			{
				pos = 0; neg = 0;
				for (int m = 0; m < 3; m++)
				{
					pos += _max_(0.0, deltaU1[m][k]);
					neg += _max_(0.0, -deltaU1[m][k]);
				}
				double thetaP, thetaM;
				if (fabs(pos) < EPS)
				{
					thetaP = 1.0;
					thetaM = 0.0;
				}
				else if (fabs(neg) < EPS)
				{
					thetaP = 0.0;
					thetaM = 1.0;
				}
				else
				{
					thetaP = _min_(1.0, neg / pos);
					thetaM = _min_(1.0, pos / neg);
				}
				d1n = thetaP * _max_(0.0, d1n) - thetaM * _max_(0.0, -d1n);
				d2n = thetaP * _max_(0.0, d2n) - thetaM * _max_(0.0, -d2n);
				d3n = thetaP * _max_(0.0, d3n) - thetaM * _max_(0.0, -d3n);
				d4n = thetaP * _max_(0.0, d4n) - thetaM * _max_(0.0, -d4n);
			}


			// вычисляем отлимитированные коэффициенты разложения
			double& xb0 = grid->cells[iCell].c.x;
			double& yb0 = grid->cells[iCell].c.y;
			double& zb0 = grid->cells[iCell].c.z;

			double& xm1 = limPm[iCell][0].x;
			double& ym1 = limPm[iCell][0].y;
			double& zm1 = limPm[iCell][0].z;

			double& xm2 = limPm[iCell][1].x;
			double& ym2 = limPm[iCell][1].y;
			double& zm2 = limPm[iCell][1].z;

			double& xm3 = limPm[iCell][2].x;
			double& ym3 = limPm[iCell][2].y;
			double& zm3 = limPm[iCell][2].z;

			double& xm4 = limPm[iCell][3].x;
			double& ym4 = limPm[iCell][3].y;
			double& zm4 = limPm[iCell][3].z;

			// коэффициенты функций принимающих значение 1.0 в центре соответствующей грани и 0.0 в центрах остальных граней
			double a1, b1, c1, d1;
			double a2, b2, c2, d2;
			double a3, b3, c3, d3;
			double a4, b4, c4, d4;

			_solve_SLE_4_(	xm1, ym1, zm1, 1.0,
							xm2, ym2, zm2, 1.0,
							xm3, ym3, zm3, 1.0,
							xm4, ym4, zm4, 1.0,
							1.0, 0.0, 0.0, 0.0,
							a1, b1, c1, d1);

			_solve_SLE_4_(	xm1, ym1, zm1, 1.0,
							xm2, ym2, zm2, 1.0,
							xm3, ym3, zm3, 1.0,
							xm4, ym4, zm4, 1.0,
							0.0, 1.0, 0.0, 0.0,
							a2, b2, c2, d2);

			_solve_SLE_4_(	xm1, ym1, zm1, 1.0,
							xm2, ym2, zm2, 1.0,
							xm3, ym3, zm3, 1.0,
							xm4, ym4, zm4, 1.0,
							0.0, 0.0, 1.0, 0.0,
							a3, b3, c3, d3);

			_solve_SLE_4_(	xm1, ym1, zm1, 1.0,
							xm2, ym2, zm2, 1.0,
							xm3, ym3, zm3, 1.0,
							xm4, ym4, zm4, 1.0,
							0.0, 0.0, 0.0, 1.0,
							a4, b4, c4, d4);

			double *pRR, *pRC;
			switch (k)
			{
			case 0:
				pRR = fROlim[iCell];
				pRC = &ROc;
				break;
			case 1:
				pRR = fRUlim[iCell];
				pRC = &RUc;
				break;
			case 2:
				pRR = fRVlim[iCell];
				pRC = &RVc;
				break;
			case 3:
				pRR = fRWlim[iCell];
				pRC = &RWc;
				break;
			case 4:
				pRR = fRElim[iCell];
				pRC = &REc;
				break;
			}

			pRR[1] = (d1n * a1 + d2n * a2 + d3n * a3 + d4n * a4) * grid->cells[iCell].HX;
			pRR[2] = (d1n * b1 + d2n * b2 + d3n * b3 + d4n * b4) * grid->cells[iCell].HY;
			pRR[3] = (d1n * c1 + d2n * c2 + d3n * c3 + d4n * c4) * grid->cells[iCell].HZ;

			if (fabs(pRR[1]) < EPS) pRR[1] = 0.0;
			if (fabs(pRR[2]) < EPS) pRR[2] = 0.0;
			if (fabs(pRR[3]) < EPS) pRR[3] = 0.0;
		}
	}

	for (int i = 0; i < cellsCount; i++)
	{
		memcpy(fRO[i], fROlim[i], funcCount * sizeof(double));
		memcpy(fRU[i], fRUlim[i], funcCount * sizeof(double));
		memcpy(fRV[i], fRVlim[i], funcCount * sizeof(double));
		memcpy(fRW[i], fRWlim[i], funcCount * sizeof(double));
		memcpy(fRE[i], fRElim[i], funcCount * sizeof(double));
	}

	calcPositivityLimiter();
}

void LimiterCockburn::chooseDirection(int* n, Point pm, double* alpha, int* neigh)
{
	Point pc = grid->cells[n[0]].c;
	Point lm = pm;
	Point l1 = grid->cells[n[1]].c;
	Point l2 = grid->cells[n[2]].c;
	Point l3 = grid->cells[n[3]].c;
	Point l4 = grid->cells[n[4]].c;

	lm -= pc;
	l1 -= pc;
	l2 -= pc;
	l3 -= pc;
	l4 -= pc;

	double a1, a2, a3;

	// перебираем комбинации векторов
	// (1,2,3)
	if (!_solve_SLE_3_(	l1.x, l2.x, l3.x,
						l1.y, l2.y, l3.y,
						l1.z, l2.z, l3.z,
						lm.x, lm.y, lm.z,
						a1, a2, a3))
	{
		grid->cells[n[0]].flag |= CELL_FLAG_BAD;
		log("WARNING: cell #%d is marked BAD...\n\n", n[0]);
		return;
	}
	if (a1 >= 0.0 && a2 >= 0.0 && a3 >= 0.0)
	{
		alpha[0] = a1; neigh[0] = n[1];
		alpha[1] = a2; neigh[1] = n[2];
		alpha[2] = a3; neigh[2] = n[3];
		return;
	}

	// (1,2,4)
	if (!_solve_SLE_3_(	l1.x, l2.x, l4.x,
						l1.y, l2.y, l4.y,
						l1.z, l2.z, l4.z,
						lm.x, lm.y, lm.z,
						a1, a2, a3))
	{
		grid->cells[n[0]].flag |= CELL_FLAG_BAD;
		log("WARNING: cell #%d is marked BAD...\n\n", n[0]);
		return;
	}
	if (a1 >= 0.0 && a2 >= 0.0 && a3 >= 0.0)
	{
		alpha[0] = a1; neigh[0] = n[1];
		alpha[1] = a2; neigh[1] = n[2];
		alpha[2] = a3; neigh[2] = n[4];
		return;
	}

	// (1,3,4)
	if (!_solve_SLE_3_(	l1.x, l3.x, l4.x,
						l1.y, l3.y, l4.y,
						l1.z, l3.z, l4.z,
						lm.x, lm.y, lm.z,
						a1, a2, a3))
	{
		grid->cells[n[0]].flag |= CELL_FLAG_BAD;
		log("WARNING: cell #%d is marked BAD...\n\n", n[0]);
		return;
	}
	if (a1 >= 0.0 && a2 >= 0.0 && a3 >= 0.0)
	{
		alpha[0] = a1; neigh[0] = n[1];
		alpha[1] = a2; neigh[1] = n[3];
		alpha[2] = a3; neigh[2] = n[4];
		return;
	}

	// (2,3,4)
	if (!_solve_SLE_3_(	l2.x, l3.x, l4.x,
						l2.y, l3.y, l4.y,
						l2.z, l3.z, l4.z,
						lm.x, lm.y, lm.z,
						a1, a2, a3))
	{
		grid->cells[n[0]].flag |= CELL_FLAG_BAD;
		log("WARNING: cell #%d is marked BAD...\n\n", n[0]);
		return;
	}
	if (a1 >= 0.0 && a2 >= 0.0 && a3 >= 0.0)
	{
		alpha[0] = a1; neigh[0] = n[2];
		alpha[1] = a2; neigh[1] = n[3];
		alpha[2] = a3; neigh[2] = n[4];
		return;
	}

	alpha[0] = 0; neigh[0] = n[0];
	alpha[1] = 0; neigh[1] = n[0];
	alpha[2] = 0; neigh[2] = n[0];
}

/*!
	Вычисление лимитера, обеспечивающего положительность плотности и давления
*/
void LimiterCockburn::calcPositivityLimiter()
{
	//double *xtt, *ytt, *ztt;
	Point *ptt;
	double Rcur1, RUcur1, RVcur1, RWcur1, Ecur1;
	double Rcur2, RUcur2, RVcur2, RWcur2, Ecur2;
	double Rcur3, RUcur3, RVcur3, RWcur3, Ecur3;
	double Rcur4, RUcur4, RVcur4, RWcur4, Ecur4;
	double sR, sRU, sRV, sRW, sE;
	double Pmin, Amin, Amax, Acur, fcur, fmax, ffcur, ffmax;
	int iPmin;
	double gamma = GAM;



	for (int k = 0; k < cellsCount; k++)
	{
		double x1 = grid->getNode(grid->cells[k].nodesInd[0]).x;
		double x2 = grid->getNode(grid->cells[k].nodesInd[1]).x;
		double x3 = grid->getNode(grid->cells[k].nodesInd[2]).x;
		double x4 = grid->getNode(grid->cells[k].nodesInd[3]).x;

		double y1 = grid->getNode(grid->cells[k].nodesInd[0]).y;
		double y2 = grid->getNode(grid->cells[k].nodesInd[1]).y;
		double y3 = grid->getNode(grid->cells[k].nodesInd[2]).y;
		double y4 = grid->getNode(grid->cells[k].nodesInd[3]).y;

		double z1 = grid->getNode(grid->cells[k].nodesInd[0]).z;
		double z2 = grid->getNode(grid->cells[k].nodesInd[1]).z;
		double z3 = grid->getNode(grid->cells[k].nodesInd[2]).z;
		double z4 = grid->getNode(grid->cells[k].nodesInd[3]).z;

		double xc = (x1 + x2 + x3 + x4) / 3.;
		double yc = (y1 + y2 + y3 + y4) / 3.;
		double zc = (z1 + z2 + z3 + z4) / 3.;

		method->getFields(Rcur1, RUcur1, RVcur1, RWcur1, Ecur1, k, x1, y1, z1);
		method->getFields(Rcur2, RUcur2, RVcur2, RWcur2, Ecur2, k, x2, y2, z2);
		method->getFields(Rcur3, RUcur3, RVcur3, RWcur3, Ecur3, k, x3, y3, z3);
		method->getFields(Rcur4, RUcur4, RVcur4, RWcur4, Ecur4, k, x4, y4, z4);

		if ((Ecur1 - (RUcur1*RUcur1 + RVcur1*RVcur1 + RWcur1*RWcur1) / (2.*Rcur1) > 0.) && (Ecur2 - (RUcur2*RUcur2 + RVcur2*RVcur2 + RWcur2*RWcur2) / (2.*Rcur2) > 0.) && (Ecur3 - (RUcur3*RUcur3 + RVcur3*RVcur3 + RWcur3*RWcur3) / (2.*Rcur3) > 0.) && (Ecur4 - (RUcur4*RUcur4 + RVcur4*RVcur4 + RWcur4*RWcur4) / (2.*Rcur4) > 0.)){
			goto lbl_Lim_1;

		}
		else{
			ptt = new Point[method->getFaceGPCount() * 4];

			// Находим все узлы квадратур Гаусса на границе ячейки
			int iGP = 0;
			for (int iF = 0; iF < 4; iF++)
			{
				int iFace = grid->cells[k].facesInd[iF];
				for (int iFaceGP = 0; iFaceGP < method->getFaceGPCount(); ++iFaceGP)
				{
					ptt[iGP] = method->getFaceGP(k, iFaceGP);
					iGP++;
				}
			}

			double eps = 1.e-13;

			double Pcur, Pmin, Rcur, Rmin;

			// корректируем плотность
			for (int iGP = 0; iGP < method->getFaceGPCount() * 4; iGP++)
			{
				Rcur = fRO[k][0] + fRO[k][1] * method->getF(1, k, ptt[iGP]) + fRO[k][2] * method->getF(2, k, ptt[iGP]) + fRO[k][3] * method->getF(3, k, ptt[iGP]);
				if (iGP == 0)
				{
					Rmin = Rcur;
				}
				else
				{
					if (Rcur < Rmin) Rmin = Rcur;
				}
			}

			double alfR = (fRO[k][0] == Rmin) ? 1.0 : MIN((fRO[k][0] - eps) / (fRO[k][0] - Rmin), 1.);
			for (int i = 1; i < 4; i++) fRO[k][i] *= alfR;

			// корректируем давление
			for (int iGP = 0; iGP < method->getFaceGPCount() * 4; iGP++)
			{
				method->getFields(sR, sRU, sRV, sRW, sE, k, ptt[iGP]);

				double se = sE / sR - 0.5 * (sRU * sRU + sRV * sRV + sRW * sRW) / (sR * sR);
				double Pcur = (gamma - 1.) * sR * se;

				if (iGP == 0) {
					Pmin = Pcur;
					iPmin = iGP;
				}
				else {
					if (Pcur < Pmin) {
						Pmin = Pcur;
						iPmin = iGP;
					}
				}
				
			}

			if (Pmin >= eps){
				goto lbl_Lim_1;
			}
			else {
				Amin = 0.;
				Amax = 1.;
				Acur = (Amax + Amin) * 0.5;

				sR = 0.;
				sRU = 0.;
				sRV = 0.;
				sRW = 0.;
				sE = 0.;

				for (int i = 1; i < 4; i++)
				{
					sR += fRO[k][i] * method->getF(i, k, ptt[iPmin]);
					sRU += fRU[k][i] * method->getF(i, k, ptt[iPmin]);
					sRV += fRV[k][i] * method->getF(i, k, ptt[iPmin]);
					sRW += fRW[k][i] * method->getF(i, k, ptt[iPmin]);
					sE += fRE[k][i] * method->getF(i, k, ptt[iPmin]);
				}

				fcur = fRE[k][0] + Acur*sE - ((fRU[k][0] + Acur*sRU)*(fRU[k][0] + Acur*sRU) + (fRV[k][0] + Acur*sRV)*(fRV[k][0] + Acur*sRV) + (fRW[k][0] + Acur * sRW) * (fRW[k][0] + Acur * sRW)) / (2.*(fRO[k][0] + Acur*sR));
				fmax = fRE[k][0] + Amax*sE - ((fRU[k][0] + Amax*sRU)*(fRU[k][0] + Amax*sRU) + (fRV[k][0] + Amax*sRV)*(fRV[k][0] + Amax*sRV) + (fRW[k][0] + Amax * sRW) * (fRW[k][0] + Amax * sRW)) / (2.*(fRO[k][0] + Amax*sR));
				ffmax = (gamma - 1.)*fmax - eps;
				ffcur = (gamma - 1.)*fcur - eps;

				double tfl = 0.;
				// метод половинного деления
				while (abs(ffmax)>eps) {
					if (tfl > 100) goto lbl_Lim_2;
					if (ffmax * ffcur > 0.)
						Amax = Acur;
					else
						Amin = Acur;
					Acur = (Amin + Amax) * 0.5;

					sR = 0.;
					sRU = 0.;
					sRV = 0.;
					sRW = 0.;
					sE = 0.;

					for (int i = 1; i < funcCount; i++)
					{
						sR += fRO[k][i] * method->getF(i, k, ptt[iPmin]);
						sRU += fRU[k][i] * method->getF(i, k, ptt[iPmin]);
						sRV += fRV[k][i] * method->getF(i, k, ptt[iPmin]);
						sRW += fRW[k][i] * method->getF(i, k, ptt[iPmin]);
						sE += fRE[k][i] * method->getF(i, k, ptt[iPmin]);
					}

					fcur = fRE[k][0] + Acur*sE - ((fRU[k][0] + Acur*sRU)*(fRU[k][0] + Acur*sRU) + (fRV[k][0] + Acur*sRV)*(fRV[k][0] + Acur*sRV) + (fRW[k][0] + Acur * sRW) * (fRW[k][0] + Acur * sRW)) / (2.*(fRO[k][0] + Acur*sR));
					fmax = fRE[k][0] + Amax*sE - ((fRU[k][0] + Amax*sRU)*(fRU[k][0] + Amax*sRU) + (fRV[k][0] + Amax*sRV)*(fRV[k][0] + Amax*sRV) + (fRW[k][0] + Amax * sRW) * (fRW[k][0] + Amax * sRW)) / (2.*(fRO[k][0] + Amax*sR));
					ffmax = (gamma - 1.)*fmax - eps;
					ffcur = (gamma - 1.)*fcur - eps;


					tfl++;
				}


			lbl_Lim_2:
				for (int i = 1; i< funcCount; i++)
				{
					fRO[k][i] *= Acur;
					fRU[k][i] *= Acur;
					fRV[k][i] *= Acur;
					fRW[k][i] *= Acur;
					fRE[k][i] *= Acur;
				}
			}
		}
	lbl_Lim_1:;
	}

}

