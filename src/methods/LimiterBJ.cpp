#include "LimiterBJ.h"

#define MIN(a,b) ((a)<(b)?(a):(b))
#define MAX(a,b) ((a)>(b)?(a):(b))

//double _det_(double a11, double a12, double a13,
//	double a21, double a22, double a23,
//	double a31, double a32, double a33)
//{
//	return a11 * a22 * a33 + a12 * a23 * a31 + a21 * a32 * a13 - a31 * a22 * a13 - a21 * a12 * a33 - a11 * a32 * a23;
//}
//
//double _det_(double a11, double a12, double a13, double a14,
//	double a21, double a22, double a23, double a24,
//	double a31, double a32, double a33, double a34,
//	double a41, double a42, double a43, double a44)
//{
//	double d1 = _det_(a21, a22, a23,
//		a31, a32, a33,
//		a41, a42, a43);
//	double d2 = _det_(a11, a12, a13,
//		a31, a32, a33,
//		a41, a42, a43);
//	double d3 = _det_(a11, a12, a13,
//		a21, a22, a23,
//		a41, a42, a43);
//	double d4 = _det_(a11, a12, a13,
//		a21, a22, a23,
//		a31, a32, a33);
//	return -a14 * d1 + a24 * d2 - a34 * d3 + a44 * d4;
//}
//
//bool _solve_SLE_3_(double  a11, double  a12, double  a13,
//	double  a21, double  a22, double  a23,
//	double  a31, double  a32, double  a33,
//	double  b1, double  b2, double  b3,
//	double& x1, double& x2, double& x3)
//{
//	double det = _det_(a11, a12, a13,
//		a21, a22, a23,
//		a31, a32, a33);
//	if (fabs(det) <= 1.e-12)
//	{
//		log("WARNING: matrix determinat is nearly zero,\n         when choosing direction for Cockburn limiter!\n");
//		x1 = 0.0;
//		x2 = 0.0;
//		x3 = 0.0;
//		return false;
//	}
//	x1 = _det_(b1, a12, a13,
//		b2, a22, a23,
//		b3, a32, a33) / det;
//	x2 = _det_(a11, b1, a13,
//		a21, b2, a23,
//		a31, b3, a33) / det;
//	x3 = _det_(a11, a12, b1,
//		a21, a22, b2,
//		a31, a32, b3) / det;
//	return true;
//}
//
//void _solve_SLE_4_(double  a11, double  a12, double  a13, double a14,
//	double  a21, double  a22, double  a23, double a24,
//	double  a31, double  a32, double  a33, double a34,
//	double  a41, double  a42, double  a43, double a44,
//	double  b1, double  b2, double  b3, double b4,
//	double& x1, double& x2, double& x3, double& x4)
//{
//	double det = _det_(a11, a12, a13, a14,
//		a21, a22, a23, a24,
//		a31, a32, a33, a34,
//		a41, a42, a43, a44);
//
//	x1 = _det_(b1, a12, a13, a14,
//		b2, a22, a23, a24,
//		b3, a32, a33, a34,
//		b4, a42, a43, a44) / det;
//
//	x2 = _det_(a11, b1, a13, a14,
//		a21, b2, a23, a24,
//		a31, b3, a33, a34,
//		a41, b4, a43, a44) / det;
//
//	x3 = _det_(a11, a12, b1, a14,
//		a21, a22, b2, a24,
//		a31, a32, b3, a34,
//		a41, a42, b4, a44) / det;
//
//	x4 = _det_(a11, a12, a13, b1,
//		a21, a22, a23, b2,
//		a31, a32, a33, b3,
//		a41, a42, a43, b4) / det;
//
//}

LimiterBJ::LimiterBJ(FEM_DG *mthd, Grid *grid, double **ro, double **ru, double **rv, double **rw, double **re, int fCount)
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

LimiterBJ::~LimiterBJ()
{
	delete[] grid;

	for (int i = 0; i < cellsCount; i++)
	{
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

	delete[] method;
}

void LimiterBJ::initLimiterParameters()
{
	/*limAlpha = new double** [cellsCount];
	limNeigh = new int** [cellsCount];
	limLm = new Vector* [cellsCount];
	limLmN = new Vector* [cellsCount];
	limPm = new Point* [cellsCount];*/

	fROlim = new double [cellsCount];
	fRUlim = new double [cellsCount];
	fRVlim = new double [cellsCount];
	fRWlim = new double [cellsCount];
	fRElim = new double [cellsCount];

	//for (int i = 0; i < cellsCount; i++)
	//{
	//	limAlpha[i] = new double*[4];
	//	limAlpha[i][0] = new double[3];
	//	limAlpha[i][1] = new double[3];
	//	limAlpha[i][2] = new double[3];
	//	limAlpha[i][3] = new double[3];

	//	limNeigh[i] = new int*[4];
	//	limNeigh[i][0] = new int[3];
	//	limNeigh[i][1] = new int[3];
	//	limNeigh[i][2] = new int[3];
	//	limNeigh[i][3] = new int[3];

	//	limLm[i] = new Vector[4];
	//	limLmN[i] = new Vector[4];
	//	limPm[i] = new Point[4];

	//	fROlim[i] = new double[funcCount];
	//	fRUlim[i] = new double[funcCount];
	//	fRVlim[i] = new double[funcCount];
	//	fRWlim[i] = new double[funcCount];
	//	fRElim[i] = new double[funcCount];

	//}

	/*deltaS1 = new double[5];
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
	}*/

	//// находим угол и коэффициенты разложения
	//for (int iCell = 0; iCell < cellsCount; iCell++)
	//{
	//	int n[5]; // 0 - текущая ячейка, 1-4 - номера соседних ячеек
	//	n[0] = iCell;
	//	for (int i = 0; i < 4; i++)
	//	{
	//		n[i + 1] = grid->cells[iCell].neigh[i];
	//	}

	//	if (n[1] > 0 && n[2] > 0 && n[3] > 0 && n[4] > 0) // не лимитируем граничные ячейки
	//	{
	//		for (int m = 0; m < 4; m++)
	//		{
	//			limPm[iCell][m] = grid->faces[grid->cells[iCell].facesInd[m]].c;
	//			limLm[iCell][m] = limPm[iCell][m] - grid->cells[iCell].c;
	//			double tmp = sqrt(limLm[iCell][m].x * limLm[iCell][m].x + limLm[iCell][m].y * limLm[iCell][m].y + limLm[iCell][m].z * limLm[iCell][m].z);
	//			limLmN[iCell][m] /= tmp;
	//			//chooseDirection(n, limPm[iCell][m], limAlpha[iCell][m], limNeigh[iCell][m]);
	//		}
	//	}
	//}
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

void LimiterBJ::run()
{
	for (int iCell = 0; iCell < cellsCount; iCell++)
	{
		/*INDEX_LIST rFaceList;
		INDEX_LIST rNodeList;
		INDEX_LIST rNodeCellList;*/

		Point pC = grid->cells[iCell].c;

		double fROm, fRUm, fRVm, fRWm, fREm;
		double fROmin, fRUmin, fRVmin, fRWmin, fREmin;
		double fROmax, fRUmax, fRVmax, fRWmax, fREmax;
		
		method->getFieldsLinear(fROm, fRUm, fRVm, fRWm, fREm, iCell, pC);
		fROmin = fROmax = fROm;
		fRUmin = fRUmax = fRUm;
		fRVmin = fRVmax = fRVm; 
		fRWmin = fRWmax = fRWm; 
		fREmin = fREmax = fREm;

		for (int iF = 0; iF < grid->cells[iCell].fCount; iF++)			// цикл по списку граней ячейки
		{
			int iFace = grid->cells[iCell].facesInd[iF];
			int iP = grid->faces[iFace].c2 == iCell ? grid->faces[iFace].c1 : grid->faces[iFace].c2;
			if (iP > -1)
			{
				Point pCp = grid->cells[iP].c;
				double fROp, fRUp, fRVp, fRWp, fREp;
				method->getFieldsLinear(fROp, fRUp, fRVp, fRWp, fREp, iCell, pCp);

				if (fROp < fROmin) fROmin = fROp;
				if (fRUp < fRUmin) fRUmin = fRUp;
				if (fRVp < fRVmin) fRVmin = fRVp;
				if (fRWp < fRWmin) fRWmin = fRWp;
				if (fREp < fREmin) fREmin = fREp;

				if (fROp > fROmax) fROmax = fROp;
				if (fRUp > fRUmax) fRUmax = fRUp;
				if (fRVp > fRVmax) fRVmax = fRVp;
				if (fRWp > fRWmax) fRWmax = fRWp;
				if (fREp > fREmax) fREmax = fREp;
			}
		}

		fROlim[iCell] = 1.0;
		fRUlim[iCell] = 1.0;
		fRVlim[iCell] = 1.0;
		fRWlim[iCell] = 1.0;
		fRElim[iCell] = 1.0;
		for (int i = 0; i < grid->cells[iCell].nCount; i++)
		{
			Point pt = grid->nodes[grid->cells[iCell].nodesInd[i]];
			double fROf, fRUf, fRVf, fRWf, fREf;
			method->getFieldsLinear(fROf, fRUf, fRVf, fRWf, fREf, iCell, pt);

			double fLR = 1.0;
			double fLU = 1.0;
			double fLV = 1.0;
			double fLW = 1.0;
			double fLE = 1.0;

			if (fROf != fROm)
			{
				if (fROf < fROm)
				{
					fLR = _min_(1.0, (fROmin - fROm) / (fROf - fROm));
				}
				else
				{
					fLR = _min_(1.0, (fROmax - fROm) / (fROf - fROm));
				}
			}
			if (fRUf != fRUm)
			{
				if (fRUf < fRUm)
				{
					fLU = _min_(1.0, (fRUmin - fRUm) / (fRUf - fRUm));
				}
				else
				{
					fLU = _min_(1.0, (fRUmax - fRUm) / (fRUf - fRUm));
				}
			}
			if (fRVf != fRVm)
			{
				if (fRVf < fRVm)
				{
					fLV = _min_(1.0, (fRVmin - fRVm) / (fRVf - fRVm));
				}
				else
				{
					fLV = _min_(1.0, (fRVmax - fRVm) / (fRVf - fRVm));
				}
			}
			if (fRWf != fRWm)
			{
				if (fRWf < fRWm)
				{
					fLW = _min_(1.0, (fRWmin - fRWm) / (fRWf - fRWm));
				}
				else
				{
					fLW = _min_(1.0, (fRWmax - fRWm) / (fRWf - fRWm));
				}
			}
			if (fREf != fREm)
			{
				if (fREf < fREm)
				{
					fLE = _min_(1.0, (fREmin - fREm) / (fREf - fREm));
				}
				else
				{
					fLE = _min_(1.0, (fREmax - fREm) / (fREf - fREm));
				}
			}

			if (fLR < fROlim[iCell])	fROlim[iCell] = fLR;
			if (fLU < fRUlim[iCell])	fRUlim[iCell] = fLU;
			if (fLV < fRVlim[iCell])	fRVlim[iCell] = fLV;
			if (fLW < fRWlim[iCell])	fRWlim[iCell] = fLW;
			if (fLE < fRElim[iCell])	fRElim[iCell] = fLE;
		}
		for (int j = 1; j < funcCount; j++) {
			fRO[iCell][j] *= fROlim[iCell];
			fRU[iCell][j] *= fRUlim[iCell];
			fRV[iCell][j] *= fRVlim[iCell];
			fRW[iCell][j] *= fRWlim[iCell];
			fRE[iCell][j] *= fRElim[iCell];
		}
	}
}

