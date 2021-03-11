#include "fvm_tvd_implicit.h"
#include "tinyxml.h"
#include <string>
#include <time.h>
#include "global.h"
#include "MatrixSolver.h"

const char* FLUX_NAMES[2] = {
		"GODUNOV",
		"LAX"
};

void FVM_TVD_IMPLICIT::init(char * xmlFileName)
{
	TiXmlDocument doc( xmlFileName );
	bool loadOkay = doc.LoadFile( TIXML_ENCODING_UTF8 );
	if (!loadOkay)
	{
		log("ERROR: %s\n", doc.ErrorDesc());
		exit(doc.ErrorId());
	}
	
	TiXmlNode* task = 0;
	TiXmlElement* el = 0;
	TiXmlNode* node0 = 0;
	TiXmlNode* node1 = 0;
	task = doc.FirstChild( "task" );

	int steadyVal = 1;
	node0 = task->FirstChild("control");
	node0->FirstChild("STEADY")->ToElement()->Attribute("value", &steadyVal);
	node0->FirstChild("TAU")->ToElement()->Attribute("value", &TAU);
	node0->FirstChild("TMAX")->ToElement()->Attribute("value", &TMAX);
	node0->FirstChild("STEP_MAX")->ToElement()->Attribute("value", &STEP_MAX);
	node0->FirstChild("FILE_OUTPUT_STEP")->ToElement()->Attribute("value", &FILE_SAVE_STEP);
	node0->FirstChild("LOG_OUTPUT_STEP")->ToElement()->Attribute("value", &PRINT_STEP);
	const char * flxStr = node0->FirstChild("FLUX")->ToElement()->Attribute("value");
	if (strcmp(flxStr, "GODUNOV") == 0) {
		FLUX = FLUX_GODUNOV;
	}
	else if (strcmp(flxStr, "LAX") == 0) {
		FLUX = FLUX_LAX;
	}
	else {
		FLUX = FLUX_GODUNOV;
	}	

	if (steadyVal == 0) {
		STEADY = false;
	} else {
		STEADY = true;
		node1 = node0->FirstChild("CFL");
		node1->FirstChild("start")->ToElement()->Attribute("value", &CFL);
		node1->FirstChild("scale")->ToElement()->Attribute("value", &scaleCFL);
		node1->FirstChild("max")->ToElement()->Attribute("value", &maxCFL);
		node1->FirstChild("step")->ToElement()->Attribute("value", &stepCFL);
		node1->FirstChild("max_limited_cells")->ToElement()->Attribute("value", &maxLimCells);
	}

	// сглаживание нев€зок
	int smUsing = 1;
	node0 = task->FirstChild("smoothing");
	node0->FirstChild("using")->ToElement()->Attribute("value", &smUsing);
	node0->FirstChild("coefficient")->ToElement()->Attribute("value", &SMOOTHING_PAR);
	SMOOTHING = (smUsing == 1);

	// чтение параметров о ѕ–≈ƒ≈Ћ№Ќџ’ «Ќј„≈Ќ»я’
	node0 = task->FirstChild("limits");
	node0->FirstChild("ro")->ToElement()->Attribute("min", &limitRmin);
	node0->FirstChild("ro")->ToElement()->Attribute("max", &limitRmax);
	node0->FirstChild("p")->ToElement()->Attribute( "min", &limitPmin);
	node0->FirstChild("p")->ToElement()->Attribute( "max", &limitPmax);
	node0->FirstChild("u")->ToElement()->Attribute( "max", &limitUmax);

	// чтение параметров о ћј“≈–»јЋј’
	node0 = task->FirstChild("materials");
	node0->ToElement()->Attribute("count", &matCount);;
	materials = new Material[matCount];
	TiXmlNode* matNode = node0->FirstChild("material"); 
	for (int i = 0; i < matCount; i++)
	{
		Material & mat = materials[i];
		matNode->ToElement()->Attribute("id", &mat.id);
		node1 = matNode->FirstChild("name");
		el = node1->ToElement();
		mat.name = el->GetText();
		node1 = matNode->FirstChild("parameters");
		node1->FirstChild( "M"  )->ToElement()->Attribute( "value", &mat.M  );
		node1->FirstChild( "Cp" )->ToElement()->Attribute( "value", &mat.Cp );
		node1->FirstChild( "K"  )->ToElement()->Attribute( "value", &mat.K  );
		node1->FirstChild( "ML" )->ToElement()->Attribute( "value", &mat.ML );
		matNode = matNode->NextSibling("material");
	}

	// чтение параметров о –≈√»ќЌј’
	node0 = task->FirstChild("regions");
	node0->ToElement()->Attribute("count", &regCount);
	regions = new Region[regCount];
	TiXmlNode* regNode = node0->FirstChild("region");
	for (int i = 0; i < regCount; i++)
	{
		Region & reg = regions[i];
		regNode->ToElement()->Attribute("id", &reg.id);
		regNode->FirstChild("material")->ToElement()->Attribute("id", &reg.matId);
		regNode->FirstChild("cell")->ToElement()->Attribute("type", &reg.cellType);
		
		node1 = regNode->FirstChild("parameters");
		node1->FirstChild( "Vx" )->ToElement()->Attribute( "value", &reg.par.u );
		node1->FirstChild( "Vy" )->ToElement()->Attribute( "value", &reg.par.v );
		node1->FirstChild( "T"  )->ToElement()->Attribute( "value", &reg.par.T );
		node1->FirstChild( "P"  )->ToElement()->Attribute( "value", &reg.par.p );
		
		Material& mat = materials[reg.matId];
		mat.URS(reg.par, 2);	// r=r(p,T)
		mat.URS(reg.par, 1);	// e=e(p,r)
		
		regNode = regNode->NextSibling("region");
	}

	// чтение параметров о √–јЌ»„Ќџ’ ”—Ћќ¬»я’
	node0 = task->FirstChild("boundaries");
	node0->ToElement()->Attribute("count", &bCount);
	boundaries = new Boundary[bCount];
	TiXmlNode* bNode = node0->FirstChild("boundCond");
	for (int i = 0; i < bCount; i++)
	{
		Boundary & b = boundaries[i];
		bNode->ToElement()->Attribute("edgeType", &b.edgeType);
		const char * str = bNode->FirstChild("type")->ToElement()->GetText();
		if (strcmp(str, "BOUND_WALL") == 0) 
		{
			b.parCount = 0;
			b.par = NULL;
			b.type = Boundary::BOUND_WALL;
		} else
		if (strcmp(str, "BOUND_OUTLET") == 0) 
		{
			b.parCount = 0;
			b.par = NULL;
			b.type = Boundary::BOUND_OUTLET;
		} else
		if (strcmp(str, "BOUND_INLET") == 0) 
		{
			b.parCount = 4;
			b.par = new double[4];
			b.type = Boundary::BOUND_INLET;

			node1 = bNode->FirstChild("parameters");
			node1->FirstChild( "T"  )->ToElement()->Attribute( "value", &b.par[0] );
			node1->FirstChild( "P"  )->ToElement()->Attribute( "value", &b.par[1] );
			node1->FirstChild( "Vx" )->ToElement()->Attribute( "value", &b.par[2] );
			node1->FirstChild( "Vy" )->ToElement()->Attribute( "value", &b.par[3] );
		} else {
			log("ERROR: unsupported boundary condition type '%s'", str);
			EXIT(1);
		}

		bNode = bNode->NextSibling("boundCond");
	}


	node0 = task->FirstChild("mesh");
	const char* fName = task->FirstChild("mesh")->FirstChild("name")->ToElement()->Attribute("value");
	grid.initFromFiles((char*)fName);

	cTau	= new double[grid.cCount];

	ro		= new double[grid.cCount];
	ru		= new double[grid.cCount];
	rv		= new double[grid.cCount];
	re		= new double[grid.cCount];

	tmpArr = new double[grid.cCount];
	tmpArrInt = new int[grid.cCount];

	for (int i = 0; i < grid.cCount; i++)
	{
		Region & reg = getRegion(i);
		convertParToCons(i, reg.par);
	}

	calcTimeStep();
	save(0);
}

void FVM_TVD_IMPLICIT::calcGrad() 
{
	static bool isCall;
	if (!isCall) {		//hardcore memory optimization.
		gradR		= new Vector[grid.cCount];
		gradP		= new Vector[grid.cCount];
		gradU		= new Vector[grid.cCount];
		gradV		= new Vector[grid.cCount];
		isCall = true;
	}

	int nc = grid.cCount;
	int ne = grid.eCount;

	memset(gradR, 0, grid.cCount*sizeof(Vector));
	memset(gradP, 0, grid.cCount*sizeof(Vector));
	memset(gradU, 0, grid.cCount*sizeof(Vector));
	memset(gradV, 0, grid.cCount*sizeof(Vector));
	
	for (int iEdge = 0; iEdge < ne; iEdge++)
	{
			
		int c1	= grid.edges[iEdge].c1;
		int c2	= grid.edges[iEdge].c2;
			
		Param pL, pR;
		reconstruct(iEdge, pL, pR);

		Vector n	= grid.edges[iEdge].n;
		double l	= grid.edges[iEdge].l;
			
		gradR[c1].x += (pL.r+pR.r)/2*n.x*l;
		gradR[c1].y += (pL.r+pR.r)/2*n.y*l;
		gradP[c1].x += (pL.p+pR.p)/2*n.x*l;
		gradP[c1].y += (pL.p+pR.p)/2*n.y*l;
		gradU[c1].x += (pL.u+pR.u)/2*n.x*l;
		gradU[c1].y += (pL.u+pR.u)/2*n.y*l;
		gradV[c1].x += (pL.v+pR.v)/2*n.x*l;
		gradV[c1].y += (pL.v+pR.v)/2*n.y*l;
		if (c2 > -1) 
		{
			gradR[c2].x -= (pL.r+pR.r)/2*n.x*l;
			gradR[c2].y -= (pL.r+pR.r)/2*n.y*l;
			gradP[c2].x -= (pL.p+pR.p)/2*n.x*l;
			gradP[c2].y -= (pL.p+pR.p)/2*n.y*l;
			gradU[c2].x -= (pL.u+pR.u)/2*n.x*l;
			gradU[c2].y -= (pL.u+pR.u)/2*n.y*l;
			gradV[c2].x -= (pL.v+pR.v)/2*n.x*l;
			gradV[c2].y -= (pL.v+pR.v)/2*n.y*l;
		}

	}
	for (int iCell = 0; iCell < nc; iCell++)
	{
		register double si = grid.cells[iCell].S;
		gradR[iCell].x /= si;
		gradR[iCell].y /= si;
		gradP[iCell].x /= si;
		gradP[iCell].y /= si;
		gradU[iCell].x /= si;
		gradU[iCell].y /= si;
		gradV[iCell].x /= si;
		gradV[iCell].y /= si;
	}
}

double **FVM_TVD_IMPLICIT::allocMtx4()
{
	double		**tempMtx4 = new double* [4];
	for (int i = 0; i < 4; ++i) tempMtx4[i] = new double [4];
	return tempMtx4;
}
void   FVM_TVD_IMPLICIT::freeMtx4(double **mtx4)
{
	for (int i = 0; i < 4; ++i)
		delete [] mtx4[i];
	delete [] mtx4;
}
void FVM_TVD_IMPLICIT::multMtx4(double **dst4, double **srcA4, double **srcB4)
{
	double sum;
	for (int i = 0; i < 4; ++i)
	{
		for (int j = 0; j < 4; ++j)
		{
			sum = 0;
			for (int k = 0; k < 4; ++k)
				sum += srcA4[i][k] * srcB4[k][j];
			dst4[i][j] = sum;
		}
	}
}
void FVM_TVD_IMPLICIT::clearMtx4(double **mtx4)
{
	for (int i = 0; i < 4; ++i)
		for (int j = 0; j < 4; ++j)
			mtx4[i][j] = 0;
}
void FVM_TVD_IMPLICIT::printMtx4(double **mtx4, char *msg)
{
	if (msg != 0)	printf("%s\n", msg);
	else			printf("\n");
	for (int i = 0; i < 4; ++i)
	{
		for (int j = 0; j < 4; ++j)
			printf("%16.8e ", mtx4[i][j]);
		printf("\n");
	}
}

void FVM_TVD_IMPLICIT::calcTimeStep()
{
/*
	double tau = 1.0e+20;
	for (int iCell = 0; iCell < grid.cCount; iCell++)
	{
		Param p;
		convertConsToPar(iCell, p);
		double tmp = grid.cells[iCell].S/_max_(abs(p.u)+p.cz, abs(p.v)+p.cz);
		if (tmp < tau) tau = tmp;
	}
	tau  *= CFL;
	TAU = _min_(TAU, tau);
	for (int iCell = 0; iCell < grid.cCount; iCell++)
		cTau[iCell] = TAU;
*/
	//double tauMin = 1.0e+20;
	//for (int iCell = 0; iCell < grid.cCount; iCell++)
	//{
	//	if (STEADY) {
	//		Param p;
	//		convertConsToPar(iCell, p);
	//		double tmp = grid.cells[iCell].S/(sqrt(p.u*p.u+p.v*p.v)+p.cz+FLT_EPSILON);
	//		cTau[iCell] = _min_(CFL*tmp, TAU);
	//		if (cTau[iCell] < tauMin) tauMin = cTau[iCell];
	//	} else {
	//		cTau[iCell] = TAU;
	//	}
	//}
	//if (STEADY)	TAU_MIN = tauMin;

	if (!STEADY) {
		for (int iCell = 0; iCell < grid.cCount; iCell++)
		{
			cTau[iCell] = TAU;
		}
	}
}

void FVM_TVD_IMPLICIT::eigenValues(double** dst4, double c, double u, double nx, double v, double ny)
{
	clearMtx4(dst4);
	double qn = u*nx + v*ny;
	dst4[0][0] = qn - c;
	dst4[1][1] = qn;
	dst4[2][2] = qn + c;
	dst4[3][3] = qn;
}
void FVM_TVD_IMPLICIT::rightEigenVector(double **dst4, double c, double u, double nx, double v, double ny, double H)
{
	double qn = u*nx + v*ny;
	double q2 = u*u + v*v;
	double lx = -ny; 
	double ly = nx;
	double ql = u*lx + v*ly;
	dst4[0][0] = 1.0;			dst4[0][1] = 1.0;		dst4[0][2] = 1.0;			dst4[0][3] = 0.0;
	dst4[1][0] = u-c*nx;	dst4[1][1] = u;		dst4[1][2] = u+c*nx;	dst4[1][3] = lx;
	dst4[2][0] = v-c*ny;	dst4[2][1] = v;		dst4[2][2] = v+c*ny;	dst4[2][3] = ly;
	dst4[3][0] = H-qn*c;	dst4[3][1] = q2/2;	dst4[3][2] = H+qn*c;	dst4[3][3] = ql;
}
void FVM_TVD_IMPLICIT::leftEigenVector(double **dst4, double c, double GAM, double u, double nx, double v, double ny)
{
	double qn = u*nx + v*ny;
	double q2 = u*u + v*v;
	double lx = -ny; 
	double ly = nx;
	double ql = u*lx + v*ly;
	double g1 = GAM-1.0;
	double c2 = c*c;
	dst4[0][0] = 0.5*(0.5*q2*g1/c2 + qn/c);	dst4[0][1] = -0.5*(g1*u/c2 + nx/c);	dst4[0][2] = -0.5*(g1*v/c2 + ny/c);	dst4[0][3] = 0.5*g1/c2;
	dst4[1][0] = 1-0.5*q2*g1/c2;			dst4[1][1] = g1*u/c2;				dst4[1][2] = g1*v/c2;				dst4[1][3] = -g1/c2;
	dst4[2][0] = 0.5*(0.5*q2*g1/c2 - qn/c);	dst4[2][1] = -0.5*(g1*u/c2 - nx/c);	dst4[2][2] = -0.5*(g1*v/c2 - ny/c);	dst4[2][3] = 0.5*g1/c2;
	dst4[3][0] = -ql;						dst4[3][1] = lx;					dst4[3][2] = ly;					dst4[3][3] = 0;	
}

void FVM_TVD_IMPLICIT::calcAP(double **dst4, double **rightEgnVecl4, double **egnVal4, double **leftEgnVecl4)
{
	double		**tempMtx4 = allocMtx4();
	
	//egnVal4 +.
	for (int i = 0; i < 4; ++i)
		for (int j = 0; j < 4; ++j)
			dst4[i][j] = egnVal4[i][j] > 0 ? egnVal4[i][j] : 0;
	
	multMtx4(tempMtx4, rightEgnVecl4, dst4);
	multMtx4(dst4, tempMtx4, leftEgnVecl4);
	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j) {
			dst4[i][j] *= 0.5;
		}
	}
	freeMtx4(tempMtx4);
}
void FVM_TVD_IMPLICIT::calcAM(double **dst4, double **rightEgnVecl4, double **egnVal4, double **leftEgnVecl4)
{
	double		**tempMtx4 = allocMtx4();
	
	//egnVal4 -.
	for (int i = 0; i < 4; ++i)
		for (int j = 0; j < 4; ++j)
			dst4[i][j] = egnVal4[i][j] < 0 ? egnVal4[i][j] : 0;
	
	multMtx4(tempMtx4, rightEgnVecl4, dst4);
	multMtx4(dst4, tempMtx4, leftEgnVecl4);
	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j) {
			dst4[i][j] *= 0.5;
		}
	}
	freeMtx4(tempMtx4);
}

void FVM_TVD_IMPLICIT::calcRoeAverage(Param& average, Param pL, Param pR, double GAM, Vector n)
{
	double WI, UN, UT;
	//double unl = pL.u*n.x+pL.v*n.y;
	//double unr = pR.u*n.x+pR.v*n.y;
	//double utl = pL.u*n.y-pL.v*n.x;
	//double utr = pR.u*n.y-pR.v*n.x;
	//rim_orig(average.r, average.e, average.p, UN, UT, WI, pL.r, pL.p, unl, utl, 0, pR.r, pR.p, unr, utr, 0, GAM);
	//	
	//double UI = UN*n.x+UT*n.y;
	//double VI = UN*n.y-UT*n.x;

	//average.u = UI;
	//average.v = VI;

	roe_orig(average.r, average.e, average.p, average.u, average.v, WI,
		pL.r, pL.p, pL.u, pL.v, 0.0,
		pR.r, pR.p, pR.u, pR.v, 0.0, GAM);

	average.cz = sqrt(GAM*average.p/average.r);
	average.E = average.e + 0.5*(average.u*average.u + average.v*average.v);
}

void FVM_TVD_IMPLICIT::reconstruct(int iFace, Param& pL, Param& pR)
{
	int		c1 = grid.edges[iFace].c1;
	int		c2 = grid.edges[iFace].c2;
	convertConsToPar(c1, pL);
	if (grid.edges[iFace].type == Edge::TYPE_INNER) {
		convertConsToPar(c2, pR);
	}
	else {
		boundaryCond(iFace, pL, pR);
	}
}

void FVM_TVD_IMPLICIT::reconstruct(int iFace, Param& pL, Param& pR, Point p)
{
	if (grid.edges[iFace].type == Edge::TYPE_INNER) 
	{
		int c1	= grid.edges[iFace].c1;
		int c2	= grid.edges[iFace].c2;
		convertConsToPar(c1, pL);
		convertConsToPar(c2, pR);
		Point &PE = p;
		Point P1 = grid.cells[c1].c;
		Point P2 = grid.cells[c2].c;
		Vector DL1;
		Vector DL2;
		DL1.x=PE.x-P1.x;
		DL1.y=PE.y-P1.y;
		DL2.x=PE.x-P2.x;
		DL2.y=PE.y-P2.y;
		pL.r+=gradR[c1].x*DL1.x+gradR[c1].y*DL1.y;
		pL.p+=gradP[c1].x*DL1.x+gradP[c1].y*DL1.y;
		pL.u+=gradU[c1].x*DL1.x+gradU[c1].y*DL1.y;
		pL.v+=gradV[c1].x*DL1.x+gradV[c1].y*DL1.y;
		pR.r+=gradR[c2].x*DL2.x+gradR[c2].y*DL2.y;
		pR.p+=gradP[c2].x*DL2.x+gradP[c2].y*DL2.y;
		pR.u+=gradU[c2].x*DL2.x+gradU[c2].y*DL2.y;
		pR.v+=gradV[c2].x*DL2.x+gradV[c2].y*DL2.y;


	} else {
		int c1	= grid.edges[iFace].c1;
		convertConsToPar(c1, pL);
		boundaryCond(iFace, pL, pR);
	}
}

void FVM_TVD_IMPLICIT::run() 
{
	int						nc = grid.cCount; // количество €чеек.
	int						ne = grid.eCount; // количество ребер.
	double					t = 0.0;
	unsigned int			step = 0;

	MatrixSolver			*solverMtx = MatrixSolver::create("HYPRE_GMRES");
	double					**eigenMtx4, **rEigenVector4, **lEigenVector4;
	double					**Amtx4P, **Amtx4M;
	double					**right4, **mtx4;

	int solverErr = 0;

	right4 = new double*[nc];
	for (int i = 0; i < nc; i++) {
		right4[i] = new double[4];
	}
	mtx4 = allocMtx4();
	eigenMtx4 = allocMtx4();
	rEigenVector4 = allocMtx4();
	lEigenVector4 = allocMtx4();
	Amtx4P = allocMtx4();
	Amtx4M = allocMtx4();

	solverMtx->init(nc, 4);

	if (STEADY)		log("Steady-State Flow\n");
	else			log("Unsteady-State Flow\n");
	log("TMAX = %e STEP_MAX = %d\n", TMAX, STEP_MAX);	
	log("Flux calculation method: %s\n", FLUX_NAMES[FLUX]);
	
	// инициализируем портрет матрицы
	log("Matrix structure initialization:\n");
	CSRMatrix::DELTA = 65536;
	for (int iEdge = 0; iEdge < ne; iEdge++) {
		int		c1 = grid.edges[iEdge].c1;
		int		c2 = grid.edges[iEdge].c2;
		solverMtx->createMatrElement(c1, c1);
		if (c2 >= 0){
			solverMtx->createMatrElement(c1, c2);
			solverMtx->createMatrElement(c2, c2);
			solverMtx->createMatrElement(c2, c1);
		}
		if (iEdge % 100 == 0) {
			log("\tfor edge: %d\n", iEdge);
		}
	}
	log("\tcomplete...\n");
	while (t < TMAX && step < STEP_MAX)
	{
		long timeStart, timeEnd;
		timeStart = clock();
		if (!STEADY) t += TAU;
		if (!solverErr) step++;

		//заполнение матрицы.
		Param		average, cellL, cellR;
		double		__GAM = 1.4;			// TODO: сделать правильное вычисление показател€ адиабаты

		solverMtx->zero();
		for (int iCell = 0; iCell < nc; iCell++){
			memset(right4[iCell], 0, 4 * sizeof(double));
		}

		memset(tmpArr, 0, grid.cCount*sizeof(double)); 
		
		calcGrad();
		for (int iEdge = 0; iEdge < ne; iEdge++) {
			int		c1 = grid.edges[iEdge].c1;
			int		c2 = grid.edges[iEdge].c2;
			double	l = grid.edges[iEdge].l;

			//сделаем нормаль внешней.
			Vector	n = grid.edges[iEdge].n;
			reconstruct(iEdge, cellL, cellR);
			calcRoeAverage(average, cellL, cellR, __GAM, n);

			// вычисл€ем спектральный радиус дл€ вычислени€ шага по времени
			if (STEADY) {
				double lambda = sqrt(average.u*average.u + average.v*average.v) + average.cz;
				lambda *= l;

				tmpArr[c1] += lambda;
				if (c2 >= 0) tmpArr[c2] += lambda;
			}

			double H = average.E + average.p / average.r;

			eigenValues(eigenMtx4, average.cz, average.u, n.x, average.v, n.y);
			rightEigenVector(rEigenVector4, average.cz, average.u, n.x, average.v, n.y, H);
			leftEigenVector(lEigenVector4, average.cz, __GAM, average.u, n.x, average.v, n.y);
			calcAP(Amtx4P, rEigenVector4, eigenMtx4, lEigenVector4);
			calcAM(Amtx4M, rEigenVector4, eigenMtx4, lEigenVector4);
			
			for (int i = 0; i < 4; ++i)
			{
				for (int j = 0; j < 4; ++j)
				{
					Amtx4M[i][j] *= l;
					Amtx4P[i][j] *= l;
				}
			}

			solverMtx->addMatrElement(c1, c1, Amtx4P);

			if (c2 >= 0) {
				solverMtx->addMatrElement(c1, c2, Amtx4M);	//du~(i,j,n+1)

				eigenValues(eigenMtx4, average.cz, average.u, -n.x, average.v, -n.y);
				rightEigenVector(rEigenVector4, average.cz, average.u, -n.x, average.v, -n.y, H);
				leftEigenVector(lEigenVector4, average.cz, __GAM, average.u, -n.x, average.v, -n.y);
				calcAP(Amtx4P, rEigenVector4, eigenMtx4, lEigenVector4);
				calcAM(Amtx4M, rEigenVector4, eigenMtx4, lEigenVector4);

				for (int i = 0; i < 4; ++i)
				{
					for (int j = 0; j < 4; ++j)
					{
						Amtx4M[i][j] *= l;
						Amtx4P[i][j] *= l;
					}
				}

				solverMtx->addMatrElement(c2, c2, Amtx4P);
				solverMtx->addMatrElement(c2, c1, Amtx4M);

			}

			double	fr2, fu2, fv2, fe2;
			calcFlux(fr2, fu2, fv2, fe2, cellL, cellR, n, __GAM);
			
			double	fr, fu, fv, fe;
			fr = fu = fv = fe = 0.0;
			for (int iGP = 1; iGP < grid.edges[iEdge].cCount; ++iGP) 
			{
				double fr1, fu1, fv1, fe1;
				reconstruct(iEdge, cellL, cellR, grid.edges[iEdge].c[iGP]);
				calcFlux(fr1, fu1, fv1, fe1, cellL, cellR, n, __GAM);
				fr += fr1;
				fu += fu1;
				fv += fv1;
				fe += fe1;
			}

			right4[c1][0] -= 0.5*l*fr;
			right4[c1][1] -= 0.5*l*fu;
			right4[c1][2] -= 0.5*l*fv;
			right4[c1][3] -= 0.5*l*fe;
			
			if (c2 >= 0) {
				right4[c2][0] += 0.5*l*fr;
				right4[c2][1] += 0.5*l*fu;
				right4[c2][2] += 0.5*l*fv;
				right4[c2][3] += 0.5*l*fe;
			}

			//log("edge = %d of %d\n", iEdge, ne);
		}

		// вычисл€ем шаги по временив €чейках по насчитанным ранее значени€м спектра
		if (STEADY) {
			for (int iCell = 0; iCell < grid.cCount; iCell++)
			{
				Param p;
				convertConsToPar(iCell, p);
				cTau[iCell] = 2.0*CFL*grid.cells[iCell].S / (tmpArr[iCell] + 1.0e-100);
			}
		}

		for (int iCell = 0; iCell < nc; ++iCell)
		{
			for (int i = 0; i < 4; ++i)
			{
				for (int j = 0; j < 4; ++j)
				{
					if (i == j) {
						mtx4[i][j] = grid.cells[iCell].S / cTau[iCell];
					}
					else {
						mtx4[i][j] = 0.0;
					}
				}
			}
			solverMtx->addMatrElement(iCell, iCell, mtx4);		//du(i, n+1)
			solverMtx->setRightElement(iCell, right4[iCell]);	//F(i,j,n)			
		}

		int maxIter = 20;
		const double eps = 1.0e-7;
		
		solverErr = solverMtx->solve(eps, maxIter);

		if (solverErr == MatrixSolver::RESULT_OK) {

			if (SMOOTHING) smoothingDelta(solverMtx->x);

			for (int cellIndex = 0, ind = 0; cellIndex < nc; cellIndex++, ind += 4)
			{
				if (cellIsLim(cellIndex))	continue;
				ro[cellIndex] += solverMtx->x[ind + 0];
				ru[cellIndex] += solverMtx->x[ind + 1];
				rv[cellIndex] += solverMtx->x[ind + 2];
				re[cellIndex] += solverMtx->x[ind + 3];

				Param par;
				convertConsToPar(cellIndex, par);
				if (par.r > limitRmax)			{ par.r = limitRmax;	setCellFlagLim(cellIndex); }
				if (par.r < limitRmin)			{ par.r = limitRmin;	setCellFlagLim(cellIndex); }
				if (fabs(par.u) > limitUmax)		{ par.u = limitUmax*par.u / fabs(par.u);	setCellFlagLim(cellIndex); }
				if (fabs(par.v) > limitUmax)		{ par.v = limitUmax*par.v / fabs(par.v);	setCellFlagLim(cellIndex); }
				if (par.p > limitPmax)			{ par.p = limitPmax;	setCellFlagLim(cellIndex); }
				if (par.p < limitPmin)			{ par.p = limitPmin;	setCellFlagLim(cellIndex); }
				if (cellIsLim(cellIndex)) 		{ par.e = par.p / ((__GAM - 1)*par.r); convertParToCons(cellIndex, par); }
			}
			remediateLimCells();
			int limCells = getLimitedCellsCount();
			if (STEADY && (limCells >= maxLimCells)) decCFL();
			
			timeEnd = clock(); 
			
			if (step % FILE_SAVE_STEP == 0)
			{
				save(step);
			}
			if (step % PRINT_STEP == 0)
			{
				calcLiftForce();
				if (!STEADY) {

					log("step: %d\ttime step: %.16f\tmax iter: %d\tlim: %d\ttime: %d ms\n", step, t, maxIter, limCells, timeEnd-timeStart);
				}
				else {
					log("step: %d\tmax iter: %d\tlim: %d\ttime: %d ms\n", step, maxIter, limCells, timeEnd - timeStart);
				}
			}

			if (STEADY && (step % stepCFL == 0)) incCFL();
		}
		else {
			if (STEADY) {
				decCFL();
			}
			else {
				solverErr = 0;
			}
		}
	}

	freeMtx4(eigenMtx4);
	freeMtx4(rEigenVector4);
	freeMtx4(lEigenVector4);
	freeMtx4(Amtx4P);
	freeMtx4(Amtx4M);
	freeMtx4(mtx4);
	for (int i = 0; i < nc; i++) {
		delete[] right4[i];
	}
	delete[] right4;
	delete solverMtx;
}

void FVM_TVD_IMPLICIT::remediateLimCells()
{
	for (int iCell = 0; iCell < grid.cCount; iCell++) 
	{
		if (cellIsLim(iCell)) 
		{
			// пересчитываем по сосед€м			
			double sRO = 0.0;
			double sRU = 0.0;
			double sRV = 0.0;
			double sRE = 0.0;
			double S   = 0.0;
			for (int i = 0; i < grid.cells[iCell].eCount; i++)
			{
				int		iEdge = grid.cells[iCell].edgesInd[i];
				int		j = grid.edges[iEdge].c2;
				if (j == iCell)	{
					//std::swap(j, grid.edges[iEdge].c1); // так нужно еще нормаль поворчивать тогда
					j = grid.edges[iEdge].c1;
				}
				if (j >= 0) {
					double  s = grid.cells[j].S;
					S += s;
					sRO += ro[j]*s;
					sRU += ru[j]*s;
					sRV += rv[j]*s;
					sRE += re[j]*s;
				} 
			}
			if (S >= TAU*TAU) {
				ro[iCell] = sRO/S;
				ru[iCell] = sRU/S;
				rv[iCell] = sRV/S;
				re[iCell] = sRE/S;
			}
			// после 0x20 итераций пробуем вернуть €чейку в счет
			grid.cells[iCell].flag += 0x010000;
			if (grid.cells[iCell].flag & 0x200000) grid.cells[iCell].flag &= 0x001110;
		}
	}
}

int FVM_TVD_IMPLICIT::getLimitedCellsCount() {
	int n = 0;
	for (int iCell = 0; iCell < grid.cCount; iCell++) {
		if ((grid.cells[iCell].flag & CELL_FLAG_LIM) > 0) n++;
	}
	return n;
}

void FVM_TVD_IMPLICIT::smoothingDelta(double* p)
{
	for (register int shift = 0; shift < 4; shift++) {
		memset(tmpArr, 0, grid.cCount*sizeof(double));
		memset(tmpArrInt, 0, grid.cCount*sizeof(int));
		for (register int iEdge = 0; iEdge < grid.eCount; iEdge++) {
			if (grid.edges[iEdge].c2 >= 0) {
				register int c1 = grid.edges[iEdge].c1;
				register int c2 = grid.edges[iEdge].c2;
				if (!cellIsLim(c2)) { tmpArr[c1] += p[c2*4 + shift]; tmpArrInt[c1]++; }
				if (!cellIsLim(c1)) { tmpArr[c2] += p[c1*4 + shift]; tmpArrInt[c2]++; }
			}
		}
		for (register int i = 0; i < grid.cCount; i++) {
			register int iCell = i*4 + shift;
			p[iCell] = (p[iCell] + SMOOTHING_PAR * tmpArr[iCell]) / (1.0 + SMOOTHING_PAR*tmpArrInt[iCell]);
		}
	}
}

void FVM_TVD_IMPLICIT::save(int step)
{
	char fName[50];

	sprintf(fName, "res_%010d.vtk", step);
	FILE * fp = fopen(fName, "w");
	fprintf(fp, "# vtk DataFile Version 2.0\n");
	fprintf(fp, "GASDIN data file\n");
	fprintf(fp, "ASCII\n");
	fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");
	
	fprintf(fp, "POINTS %d float\n", grid.nCount);
	for (int i = 0; i < grid.nCount; i++)
	{
		fprintf(fp, "%f %f %f  ", grid.nodes[i].x,  grid.nodes[i].y, 0.0);
		if (i+1 % 8 == 0) fprintf(fp, "\n");
	}
	fprintf(fp, "\n");
	fprintf(fp, "CELLS %d %d\n", grid.cCount, 4*grid.cCount);
	for (int i = 0; i < grid.cCount; i++)
	{
		fprintf(fp, "3 %d %d %d\n", grid.cells[i].nodesInd[0], grid.cells[i].nodesInd[1], grid.cells[i].nodesInd[2]);
	}
	fprintf(fp, "\n");

	fprintf(fp, "CELL_TYPES %d\n", grid.cCount);
	for (int i = 0; i < grid.cCount; i++) fprintf(fp, "5\n");
	fprintf(fp, "\n");

	fprintf(fp, "CELL_DATA %d\nSCALARS Density float 1\nLOOKUP_TABLE default\n", grid.cCount);
	for (int i = 0; i < grid.cCount; i++)
	{
		Param p;
		convertConsToPar(i, p);
		fprintf(fp, "%25.16f ", p.r);
		if (i+1 % 8 == 0 || i+1 == grid.cCount) fprintf(fp, "\n");
	}

	fprintf(fp, "SCALARS Pressure float 1\nLOOKUP_TABLE default\n", grid.cCount);
	for (int i = 0; i < grid.cCount; i++)
	{
		Param p;
		convertConsToPar(i, p);
		fprintf(fp, "%25.16f ", p.p);
		if (i+1 % 8 == 0 || i+1 == grid.cCount) fprintf(fp, "\n");
	}

	fprintf(fp, "SCALARS Temperature float 1\nLOOKUP_TABLE default\n", grid.cCount);
	for (int i = 0; i < grid.cCount; i++)
	{
		Param p;
		convertConsToPar(i, p);
		fprintf(fp, "%25.16f ", p.T);
		if (i+1 % 8 == 0 || i+1 == grid.cCount) fprintf(fp, "\n");
	}

	fprintf(fp, "SCALARS MachNumber float 1\nLOOKUP_TABLE default\n", grid.cCount);
	for (int i = 0; i < grid.cCount; i++)
	{
		Param p;
		convertConsToPar(i, p);
		fprintf(fp, "%25.16f ", sqrt(p.u*p.u+p.v*p.v)/p.cz);
		if (i+1 % 8 == 0 || i+1 == grid.cCount) fprintf(fp, "\n");
	}

	fprintf(fp, "VECTORS Velosity float\n");
	for (int i = 0; i < grid.cCount; i++)
	{
		Param p;
		convertConsToPar(i, p);
		fprintf(fp, "%25.16f %25.16f %25.16f ", p.u, p.v, 0.0);
		if (i + 1 % 8 == 0 || i + 1 == grid.cCount) fprintf(fp, "\n");
	}
	
	fprintf(fp, "SCALARS Total_energy float 1\nLOOKUP_TABLE default\n");
	for (int i = 0; i < grid.cCount; i++)
	{
		Param p;
		convertConsToPar(i, p);
		fprintf(fp, "%25.16f ", p.E);
		if (i + 1 % 8 == 0 || i + 1 == grid.cCount) fprintf(fp, "\n");	
	}

	fprintf(fp, "SCALARS Total_pressure float 1\nLOOKUP_TABLE default\n");
	for (int i = 0; i < grid.cCount; i++)
	{
		Material &mat = getMaterial(i);
		double gam = mat.getGamma();
		double agam = gam - 1.0;
		Param p;
		convertConsToPar(i, p);
		double M2 = (p.u*p.u + p.v*p.v) / (p.cz*p.cz);
		fprintf(fp, "%25.16e ", p.p*::pow(1.0+0.5*M2*agam, gam/agam) );
		if ((i+1) % 8 == 0  ||  i+1 == grid.cCount) fprintf(fp, "\n");
	}

	fprintf(fp, "SCALARS TAU float 1\nLOOKUP_TABLE default\n", grid.cCount);
	for (int i = 0; i < grid.cCount; i++)
	{
		fprintf(fp, "%25.16f ", cTau[i]);
		if (i+1 % 8 == 0 || i+1 == grid.cCount) fprintf(fp, "\n");
	}

	fclose(fp);
	printf("File '%s' saved...\n", fName);
}

void FVM_TVD_IMPLICIT::calcFlux(double& fr, double& fu, double& fv, double& fe, Param pL, Param pR, Vector n, double GAM)
{
	if (FLUX == FLUX_GODUNOV) {	// GODUNOV FLUX
		double RI, EI, PI, UI, VI, WI, UN, UT;
		//double unl = pL.u*n.x+pL.v*n.y;
		//double unr = pR.u*n.x+pR.v*n.y;
		//rim(RI, EI, PI, UN, UI, VI,  pL.r, pL.p, unl, pL.u, pL.v,  pR.r, pR.p, unr, pR.u, pR.v, n, GAM);
			
		double unl = pL.u*n.x+pL.v*n.y;
		double unr = pR.u*n.x+pR.v*n.y;
		double utl = pL.u*n.y-pL.v*n.x;
		double utr = pR.u*n.y-pR.v*n.x;
		rim_orig(RI, EI, PI, UN, UT, WI,  pL.r, pL.p, unl, utl, 0,  pR.r, pR.p, unr, utr, 0, GAM);
		
		UI = UN*n.x+UT*n.y;
		VI = UN*n.y-UT*n.x;

		fr = RI*UN;
		fu = fr*UI+PI*n.x;
		fv = fr*VI+PI*n.y;
		fe = (RI*(EI+0.5*(UI*UI+VI*VI))+PI)*UN;
	}
	if (FLUX == FLUX_LAX) {	// LAX-FRIEDRIX FLUX
		double unl = pL.u*n.x+pL.v*n.y;
		double unr = pR.u*n.x+pR.v*n.y;
		double rol, rul, rvl, rel,  ror, rur, rvr, rer;
		double alpha = _max_(fabs(unl)+sqrt(GAM*pL.p/pL.r), fabs(unr)+sqrt(GAM*pR.p/pR.r));
		//pL.getToCons(rol, rul, rvl, rel);
		//pR.getToCons(ror, rur, rvr, rer);
		rol = pL.r;
		rul = pL.r*pL.u;
		rvl = pL.r*pL.v;
		rel = pL.p / (GAM - 1.0) + 0.5*pL.r*(pL.u*pL.u + pL.v*pL.v);
		ror = pR.r;
		rur = pR.r*pR.u;
		rvr = pR.r*pR.v;
		rer = pR.p /(GAM - 1.0) + 0.5*pR.r*(pR.u*pR.u + pR.v*pR.v);
		double frl = rol*unl;
		double frr = ror*unr;
		fr = 0.5*(frr+frl								- alpha*(ror-rol));
		fu = 0.5*(frr*pR.u+frl*pL.u + (pR.p+pL.p)*n.x	- alpha*(rur-rul));
		fv = 0.5*(frr*pR.v+frl*pL.v + (pR.p+pL.p)*n.y	- alpha*(rvr-rvl));
		fe = 0.5*((rer+pR.p)*unr + (rel+pL.p)*unl		- alpha*(rer-rel));
	}
}

void FVM_TVD_IMPLICIT::boundaryCond(int iEdge, Param& pL, Param& pR)
{
	int iBound = -1;
	for (int i = 0; i < bCount; i++)
	{
		if (grid.edges[iEdge].type == boundaries[i].edgeType) 
		{
			iBound = i;
			break;
		}
	}
	if (iBound < 0)
	{
		log("ERROR (boundary condition): unknown edge type of edge %d...\n", iEdge);
		EXIT(1);
	}
	Boundary& b = boundaries[iBound];
	int c1	= grid.edges[iEdge].c1;
	Material& m = getMaterial(c1);
	switch (b.type)
	{
	case Boundary::BOUND_INLET:
		pR.T  = b.par[0];		//!< температура
		pR.p  = b.par[1];		//!< давление
		pR.u  = b.par[2];		//!< перва€ компонента вектора скорости
		pR.v  = b.par[3];		//!< втора€ компонента вектора скорости
		
		m.URS(pR, 2);
		m.URS(pR, 1);
		break;
	
	case Boundary::BOUND_OUTLET:
		pR = pL;
		break;
	
	case Boundary::BOUND_WALL:
		pR = pL;
		double Un = pL.u*grid.edges[iEdge].n.x+pL.v*grid.edges[iEdge].n.y;
		Vector V;
		V.x = grid.edges[iEdge].n.x*Un*2.0; 
		V.y = grid.edges[iEdge].n.y*Un*2.0;
		pR.u = pL.u-V.x;
		pR.v = pL.v-V.y;
		break;
	}
}

void FVM_TVD_IMPLICIT::done()
{
	delete[] ro;
	delete[] ru;
	delete[] rv;
	delete[] re;
	
	delete[] gradR;
	delete[] gradP;
	delete[] gradU;
	delete[] gradV;

	delete[] cTau;
	delete[] tmpArr;
	delete[] tmpArrInt;
}

Region & FVM_TVD_IMPLICIT::getRegionByCellType(int type)
{
	for (int i = 0; i < regCount; i++)
	{
		if (regions[i].cellType == type) return regions[i];
	}
	log("ERROR: unknown cell type %d...\n", type);
	EXIT(1);
}


Region   &	FVM_TVD_IMPLICIT::getRegion	(int iCell)
{
	return getRegionByCellType( grid.cells[iCell].type );
}

Material &	FVM_TVD_IMPLICIT::getMaterial	(int iCell)
{
	Region & reg = getRegion(iCell);
	return materials[reg.matId];
}


void FVM_TVD_IMPLICIT::convertParToCons(int iCell, Param & par)
{
	ro[iCell] = par.r;
	ru[iCell] = par.r*par.u;
	rv[iCell] = par.r*par.v;
	re[iCell] = par.r*(par.e+0.5*(par.u*par.u+par.v*par.v));
}

void FVM_TVD_IMPLICIT::convertConsToPar(int iCell, Param & par)
{
	par.r = ro[iCell];
	par.u = ru[iCell]/ro[iCell];
	par.v = rv[iCell]/ro[iCell];
	par.E = re[iCell]/ro[iCell];
	par.e = par.E-0.5*(par.u*par.u+par.v*par.v);
	Material& mat = getMaterial(iCell);
	mat.URS(par, 0);
	mat.URS(par, 1);
}

void FVM_TVD_IMPLICIT::decCFL() 
{
	if (CFL > 0.1) {
		CFL *= 0.75;
		if (CFL < 0.1) CFL = 0.1;
		log(" CFL Number has been decreased : %f \n", CFL);
	}
}

void FVM_TVD_IMPLICIT::incCFL() 
{
	if (CFL < maxCFL) {
		CFL *= scaleCFL;
		if (CFL > maxCFL) CFL = maxCFL;
		log(" CFL Number has been increased : %f \n", CFL);
	}
}

void FVM_TVD_IMPLICIT::calcLiftForce()
{
	const double width = 1.0; // предполагаема€ ширина профил€ по z.
	Param		 par;
	Fx = Fy = 0.0;
	for (int iEdge = 0; iEdge < grid.eCount; ++iEdge)
	{
		if (grid.edges[iEdge].type == Edge::TYPE_WALL)
		{
			int			cellIndex = grid.edges[iEdge].c1;
			double		nx = grid.edges[iEdge].n.x;
			double		ny = grid.edges[iEdge].n.y;
			if (cellIndex < 0) {
				cellIndex = grid.edges[iEdge].c2;
				nx *= -1.0;
				ny *= -1.0;
			}
			convertConsToPar(cellIndex, par);
			Fx += par.p*width*grid.edges[iEdge].l*nx;
			Fy += par.p*width*grid.edges[iEdge].l*ny;
		}
	}
}


