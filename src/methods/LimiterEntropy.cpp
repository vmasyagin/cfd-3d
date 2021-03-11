#include "LimiterEntropy.h"

#define MIN(a,b) ((a)<(b)?(a):(b))
#define MAX(a,b) ((a)>(b)?(a):(b))

LimiterEntropy::LimiterEntropy(FEM_DG *mthd, Grid *grid, double **ro, double **ru, double **rv, double **rw, double **re, int fCount, LimiterParam* par)
{
	this->method = mthd;
	this->grid = grid;
	cellsCount = grid->cCount;
	funcCount = fCount;
	Pmin = par->Pmin;
	PRef = par->PRef;
	RRef = par->RRef;
	initLimiterParameters();

	fRO = ro;
	fRU = ru;
	fRV = rv;
	fRW = rw;
	fRE = re;
}

double LimiterEntropy::getDsDu(int iCell, int id, Point pt) {
	double sRO, sRU, sRV, sRW, sRE;
	method->getFields(sRO, sRU, sRV, sRW, sRE, iCell, pt);
	double qP = sRE - 0.5 * (_pow_2_(sRU) + _pow_2_(sRV) + _pow_2_(sRW)) / sRO;
	switch (id) {
	case 0: 
		return log((GAM - 1.) * qP / PRef) + sRO * ((0.5 * (_pow_2_(sRU) + _pow_2_(sRV) + _pow_2_(sRW)) / _pow_2_(sRO) / qP) - GAM / sRO) - GAM * log(sRO / RRef);
		break;
	case 1:
		return -sRU / qP;
		break;
	case 2:
		return -sRV / qP;
		break;
	case 3:
		return -sRW / qP;
		break;
	case 4:
		return sRO / qP;
		break;
	default:
		return 0.0;
		break;
	}
}

double LimiterEntropy::getDpDu(int iCell, int id, Point pt)
{
	double sRO, sRU, sRV, sRW, sRE;
	method->getFields(sRO, sRU, sRV, sRW, sRE, iCell, pt);
	double qR = (GAM - 1.) / sRO;
	switch (id) {
	case 0:
		return 0.5 * qR * (_pow_2_(sRU) + _pow_2_(sRV) + _pow_2_(sRW)) / sRO;
		break;
	case 1:
		return -sRU * qR;
		break;
	case 2:
		return -sRV * qR;
		break;
	case 3:
		return -sRW * qR;
		break;
	case 4:
		return (GAM - 1.);
		break;
	default:
		return 0.0;
		break;
	}
}

double LimiterEntropy::calcDpsiDmu(int iCell, double lambda)
{
	double intDs = 0.;
	for (int iGP = 0; iGP < method->getCellGPCount(); iGP++) {
		double gw = method->getCellGW(iCell, iGP);
		Point gp = method->getCellGP(iCell, iGP);
		double dS = 0., dUi;
		for (int i = 0; i < 5; i++) {
			dUi = 0.;
			double* fld = getField(i, iCell);
			for (int j = 1; j < funcCount; j++) {
				dUi += fld[j] * method->getF(j, iCell, gp);
			}
			dUi *= lambda;
			dS += getDsDu(iCell, i, gp) * dUi;
		}
		intDs += dS * gw;
	}
	intDs *= method->getCellJ(iCell);
	return intDs;
}

double LimiterEntropy::calcMinDphiDlambda(int iCell)
{
	int sign[3];
	double min = DBL_MAX;
	for (int iApex = 0; iApex < 4; iApex++) {
		for (int i = 0; i < 2; i++) {
			sign[0] = i;
			for (int j = 0; j < 2; j++) {
				sign[1] = j;
				for (int k = 0; k < 2; k++) {
					sign[2] = k;
					double fP = 0., dUi;
					for (int iField = 0; iField < 5; iField++) {
						dUi = 0.;
						double* fld = getField(iField, iCell);
						for (int z = 1; z < funcCount; z++) {
							dUi += pow(-1., sign[z - 1]) * fld[z] * method->getF(z, iCell, grid->getNode(grid->cells[iCell].nodesInd[iApex]));
						}
						
						fP += getDpDu(iCell, iField, grid->getNode(grid->cells[iCell].nodesInd[iApex])) * dUi;
					}
					if (min > fP) min = fP;
				}
			}
		}
	}
	return min;
}

double* LimiterEntropy::getField(int id, int iCell)
{
	switch (id) {
	case 0:
		return fRO[iCell];
		break;
	case 1:
		return fRU[iCell];
		break;
	case 2:
		return fRV[iCell];
		break;
	case 3:
		return fRW[iCell];
		break;
	case 4:
		return fRE[iCell];
		break;
	default:
		return nullptr;
		break;
	}
}

LimiterEntropy::~LimiterEntropy()
{
	delete[] grid;
	delete[] Cs;

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

	delete[] method;
}

void LimiterEntropy::initLimiterParameters()
{
	Cs = new double[cellsCount];
}

double LimiterEntropy::calcPsi(int iCell, double lambda, double mu)
{
	double Rcur, RUcur, RVcur, RWcur, Ecur, Pcur;
	double res = 0.;
	for (int iGP = 0; iGP < method->getCellGPCount(); iGP++) {
		double gw = method->getCellGW(iCell, iGP);
		Point gp = method->getCellGP(iCell, iGP);
		Rcur = fRO[iCell][0];
		RUcur = fRU[iCell][0];
		RVcur = fRV[iCell][0];
		RWcur = fRW[iCell][0];
		Ecur = fRE[iCell][0];
		for (int i = 1; i < funcCount; i++) {
			Rcur += lambda * fRO[iCell][i] * method->getF(i, iCell, gp);
			RUcur += lambda * mu * fRU[iCell][i] * method->getF(i, iCell, gp);
			RVcur += lambda * mu * fRV[iCell][i] * method->getF(i, iCell, gp);
			RWcur += lambda * mu * fRW[iCell][i] * method->getF(i, iCell, gp);
			Ecur += lambda * mu * fRE[iCell][i] * method->getF(i, iCell, gp);
		}
		Pcur = (GAM - 1.) * (Ecur - (RUcur * RUcur + RVcur * RVcur + RWcur * RWcur) / (2. * Rcur));
		res += Rcur * (log(Pcur / PRef) - GAM * log(Rcur / RRef)) * gw;
	}
	res *= method->getCellJ(iCell);
	return res;
}

double LimiterEntropy::calcMinPhi(int iCell, double lambda)
{
	double Rcur, RUcur, RVcur, RWcur, Ecur, Pcur;
	int sign[3];
	double min = DBL_MAX;
	for (int iApex = 0; iApex < 4; iApex++) {
		for (int i = 0; i < 2; i++) {
			sign[0] = i;
			for (int j = 0; j < 2; j++) {
				sign[1] = j;
				for (int k = 0; k < 2; k++) {
					sign[2] = k;
					Rcur = fRO[iCell][0];
					RUcur = fRU[iCell][0];
					RVcur = fRV[iCell][0];
					RWcur = fRW[iCell][0];
					Ecur = fRE[iCell][0];
					for (int z = 1; z < funcCount; z++) {
						Rcur += pow(-1., sign[z - 1]) * lambda * fRO[iCell][z] * method->getF(z, iCell, grid->getNode(grid->cells[iCell].nodesInd[iApex]));
						RUcur += pow(-1., sign[z - 1]) * lambda * fRU[iCell][z] * method->getF(z, iCell, grid->getNode(grid->cells[iCell].nodesInd[iApex]));
						RVcur += pow(-1., sign[z - 1]) * lambda * fRV[iCell][z] * method->getF(z, iCell, grid->getNode(grid->cells[iCell].nodesInd[iApex]));
						RWcur += pow(-1., sign[z - 1]) * lambda * fRW[iCell][z] * method->getF(z, iCell, grid->getNode(grid->cells[iCell].nodesInd[iApex]));
						Ecur += pow(-1., sign[z - 1]) * lambda * fRE[iCell][z] * method->getF(z, iCell, grid->getNode(grid->cells[iCell].nodesInd[iApex]));
					}
					Pcur = (GAM - 1.) * (Ecur - (RUcur * RUcur + RVcur * RVcur + RWcur * RWcur) / (2. * Rcur));
					if (min > Pcur) min = Pcur;
				}
			}
		}
	}
	return min;
}

void LimiterEntropy::calcCs()
{
	memset(Cs, 0, cellsCount * sizeof(double));
	/*volume integrals*/
	double Rcur, RUcur, RVcur, RWcur, Ecur, Pcur;
	for (int iCell = 0; iCell < cellsCount; iCell++) {
		double fCs = 0.;
		for (int iGP = 0; iGP < method->getCellGPCount(); iGP++) {
			double gw = method->getCellGW(iCell, iGP);
			Point gp = method->getCellGP(iCell, iGP);
			method->getFields(Rcur, RUcur, RVcur, RWcur, Ecur, iCell, gp);
			Pcur = (GAM - 1.) * (Ecur - (RUcur * RUcur + RVcur * RVcur + RWcur * RWcur) / (2. * Rcur));
			fCs += Rcur * (log(Pcur / PRef) - GAM * log(Rcur / RRef)) * gw;
		}
		fCs *= method->getCellJ(iCell);
		Cs[iCell] += fCs;
	}

	/*surf integrals*/
	for (int iFace = 0; iFace < grid->fCount; iFace++) {
		int c1 = grid->faces[iFace].c1;
		int c2 = grid->faces[iFace].c2;
		Vector n = grid->faces[iFace].n;
		double FCs = 0.;
		if (c2 > -1) {
			
			for (int iGP = 0; iGP < method->getFaceGPCount(); iGP++) {
				double fRO, fRU, fRV, fRW, fRE;
				double FR, FU, FV, FW;// , FE;

				Point& gp = method->getFaceGP(iFace, iGP);
				double gw = method->getFaceGW(iFace, iGP);

				method->getFields(fRO, fRU, fRV, fRW, fRE, c1, gp);
				Param par1;
				par1.r = fRO;
				par1.u = fRU / fRO;
				par1.v = fRV / fRO;
				par1.w = fRW / fRO;
				par1.E = fRE / fRO;
				par1.e = par1.E - par1.U2() * 0.5;
				Material& mat1 = method->getMaterial(c1);
				mat1.URS(par1, 0); // p=p(r,e)

				method->getFields(fRO, fRU, fRV, fRW, fRE, c2, gp);
				Param par2;
				par2.r = fRO;
				par2.u = fRU / fRO;
				par2.v = fRV / fRO;
				par2.w = fRW / fRO;
				par2.E = fRE / fRO;
				par2.e = par2.E - par2.U2() * 0.5;
				Material& mat2 = method->getMaterial(c2);
				mat2.URS(par2, 0); // p=p(r,e)

				//method->calcFlux(FR, FU, FV, FW, FE, par1, par2, n, GAM);
				Param pF;
				//method->calcRoeAverage(average, par1, par2, GAM, n);
				method->calcFluxFields(pF, par1, par2, n, GAM);
				double FP = pF.p;// (GAM - 1.)* (FE - (FU * FU + FV * FV * FW * FW) / (2. * FR));
				FR = pF.r;
				FU = FR * pF.u;
				FV = FR * pF.v;
				FW = FR * pF.w;
				double FHx = FU * (log(FP / PRef) - GAM * log(FR / RRef));
				double FHy = FV * (log(FP / PRef) - GAM * log(FR / RRef));
				double FHz = FW * (log(FP / PRef) - GAM * log(FR / RRef));
				FCs += method->getTau() * (FHx + FHy + FHz) * gw;
			}
			FCs *= method->getFaceJ(iFace);
			Cs[c1] += FCs;
			Cs[c2] -= FCs;
		}
		else {
			for (int iGP = 0; iGP < method->getFaceGPCount(); iGP++) {
				double fRO, fRU, fRV, fRW, fRE;
				double FR, FU, FV, FW;// , FE;

				Point& gp = method->getFaceGP(iFace, iGP);
				double gw = method->getFaceGW(iFace, iGP);

				method->getFields(fRO, fRU, fRV, fRW, fRE, c1, gp);
				Param par1;
				par1.r = fRO;
				par1.u = fRU / fRO;
				par1.v = fRV / fRO;
				par1.w = fRW / fRO;
				par1.E = fRE / fRO;
				par1.e = par1.E - par1.U2() * 0.5;
				Material& mat = method->getMaterial(c1);
				mat.URS(par1, 0); // p=p(r,e)
				mat.URS(par1, 1);

				Param par2;
				method->boundaryCond(iFace, par1, par2);

				//method->calcFlux(FR, FU, FV, FW, FE, par1, par2, grid->faces[iFace].n, GAM);
				Param pF;
				//method->calcRoeAverage(average, par1, par2, GAM, n);
				method->calcFluxFields(pF, par1, par2, n, GAM);
				double FP = pF.p;// (GAM - 1.)* (FE - (FU * FU + FV * FV * FW * FW) / (2. * FR));
				FR = pF.r;
				FU = FR * pF.u;
				FV = FR * pF.v;
				FW = FR * pF.w;
				double FHx = FU * (log(FP / PRef) - GAM * log(FR / RRef));
				double FHy = FV * (log(FP / PRef) - GAM * log(FR / RRef));
				double FHz = FW * (log(FP / PRef) - GAM * log(FR / RRef));
				FCs += method->getTau() * (FHx + FHy + FHz) * gw;
			}
			FCs *= method->getFaceJ(iFace);
			Cs[c1] += FCs;
		}
	}
}

void LimiterEntropy::run()
{
	//double Rcur, RUcur, RVcur, RWcur, Ecur, Pcur;
	double Lp, Ms;
	double eps = 1.e-5;
	double gamma = GAM;

	calcCs();

	for (int k = 0; k < cellsCount; k++)
	{
		int n0 = grid->cells[k].neigh[0];
		int n1 = grid->cells[k].neigh[1];
		int n2 = grid->cells[k].neigh[2];
		int n3 = grid->cells[k].neigh[3];
		if (!(n0 > -1 && n1 > -1 && n2 > -1 && n3 > -1)) // в граничных ячейках зануляем все коэффициенты кроме нулевого
		{
			for (int j = 1; j < 4; j++) {
				fRO[k][j] = 0.0;
				fRU[k][j] = 0.0;
				fRV[k][j] = 0.0;
				fRW[k][j] = 0.0;
				fRE[k][j] = 0.0;
			}
			continue;
		}
		bool limSkip = false;
		// calc lambda_p
		if (calcMinPhi(k, 1.) >= Pmin) {
			Lp = 1.0;
			limSkip = true;
		}
		else {
			// Newton method
			double lambda0 = 1.0, lambda1;
			do
			{
				double corr = (calcMinPhi(k, lambda0) - Pmin) / (calcMinDphiDlambda(k) - Pmin);
				lambda1 = lambda0 - corr;
				if ((lambda1 > 1.) || (lambda1 < 0.)) {
					int tfl = 0;
					// метод половинного деления
					while ((lambda1 > 1.) || (lambda1 < 0.)) {
						if (tfl > 100) break;
						lambda1 = (lambda0 + lambda1) * 0.5;
						tfl++;
					}
				}
				lambda0 = lambda1;
			} while (fabs(calcMinPhi(k, lambda1) - Pmin) > eps);
			Lp = lambda0;
		}

		// calc mu_s

		if (calcPsi(k, Lp, 1.) >= Cs[k]) {
			Ms = 1.0;
		}
		else {
			limSkip = false;
			// Newton method
			double mu0 = 1.0, mu1;
			do
			{
				double corr = (calcPsi(k, Lp, mu0) - Cs[k]) / (calcDpsiDmu(k, Lp) - Cs[k]);
				mu1 = mu0 - corr;
				if ((mu1 > 1.) || (mu1 < 0.)) {
					int tfl = 0;
					// метод половинного деления
					while ((mu1 > 1.) || (mu1 < 0.)) {
						if (tfl > 100) break;
						mu1 = (mu0 + mu1) * 0.5;
						tfl++;
					}
				}
				mu0 = mu1;
			} while (fabs(calcPsi(k, Lp, mu1) - Cs[k]) > eps);
			Ms = mu0;
		}
		if (limSkip) goto lbl_Lim_1;

	lbl_Lim_2:
		for (int i = 1; i < funcCount; i++)
		{
			fRO[k][i] *= Lp * Ms;
			fRU[k][i] *= Lp * Ms;
			fRV[k][i] *= Lp * Ms;
			fRW[k][i] *= Lp * Ms;
			fRE[k][i] *= Lp * Ms;
		}
	lbl_Lim_1:;
	}
}

