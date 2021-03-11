#include "fem_rkdg.h"
#include "tinyxml.h"
#include <string>
#include "global.h"
#include "MeshReader.h"

#define FLUX_RIM

#define _EPS_ 1.0e-11
#define _NORM_(X) ((fabs(X) <= _EPS_) ? 0.0 : (X))
#define POW_2(X) ((X)*(X))

VECTOR normVector(VECTOR v)
{
	VECTOR w(v.n);
	for (int i = 0; i < v.n; i++) {
		w[i] = _NORM_(v[i]);
	}
	return w;
}


void FEM_RKDG::init(char * xmlFileName) // TODO: освободить память для всех массивов
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


	node0 = task->FirstChild("control");
	node0->FirstChild("STEADY")->ToElement()->Attribute("value", &STEADY);
	node0->FirstChild("TAU")->ToElement()->Attribute("value", &TAU);
	node0->FirstChild("TMAX")->ToElement()->Attribute("value", &TMAX);
	node0->FirstChild("CFL")->ToElement()->Attribute("value", &CFL);
	node0->FirstChild("FILE_OUTPUT_STEP")->ToElement()->Attribute("value", &FILE_SAVE_STEP);
	node0->FirstChild("LOG_OUTPUT_STEP")->ToElement()->Attribute("value", &PRINT_STEP);
	node0->FirstChild("BASE_FUNC_COUNT")->ToElement()->Attribute("value", &FUNC_COUNT);
	node0->FirstChild("GP_CELL_COUNT")->ToElement()->Attribute("value", &GP_CELL_COUNT);
	node0->FirstChild("GP_EDGE_COUNT")->ToElement()->Attribute("value", &GP_EDGE_COUNT);
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

	// чтение параметров о МАТЕРИАЛАХ
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

	// чтение параметров о РЕГИОНАХ
	double maxP = 0.0;
	double maxR = 0.0;
	double maxT = 0.0;
	double maxU = 0.0;

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

		reg.name = regNode->FirstChild("name")->ToElement()->GetText();
		
		node1 = regNode->FirstChild("parameters");
		node1->FirstChild( "Vx" )->ToElement()->Attribute( "value", &reg.par.u );
		node1->FirstChild( "Vy" )->ToElement()->Attribute( "value", &reg.par.v );
		node1->FirstChild( "T"  )->ToElement()->Attribute( "value", &reg.par.T );
		node1->FirstChild( "P"  )->ToElement()->Attribute( "value", &reg.par.p );
		
		Material& mat = materials[reg.matId];
		mat.URS(reg.par, 2);	// r=r(p,T)
		mat.URS(reg.par, 1);	// e=e(p,r)
		
		regNode = regNode->NextSibling("region");

		if (reg.par.p > maxP) maxP = reg.par.p;
		if (reg.par.r > maxR) maxR = reg.par.r;
		if (reg.par.T > maxT) maxT = reg.par.T;
		if (reg.par.U2() > maxU) maxU = reg.par.U2();
	}

	maxU = sqrt(maxU);

	//// параметры обезразмеривания
	//L_ = 1.0;
	//R_ = maxR;					// характерная плотность = начальная плотность
	//P_ = maxP;					// характерное давление = начальное давление 
	//T_ = maxT;					// характерная температура = начальная температура
	//U_ = sqrt(P_ / R_);		// характерная скорость = sqrt( P_ / R_ )
	//E_ = POW_2(U_);			// характерная энергия  = U_**2
	//CV_ = POW_2(U_) / T_;	// характерная теплоёмкость  = U_**2 / T_
	//TIME_ = L_ / U_;			// характерное время
	//MU_ = R_ * U_ * L_;		// характерная вязкость = R_ * U_ * L_
	//KP_ = R_ * POW_2(U_) * U_ * L_ / T_;	// коэффициент теплопроводности = R_ * U_**3 * L_ / T_
	//CV_ = POW_2(U_) / T_;	// характерная теплоёмкость  = U_**2 / T_


	//// Обезразмеривание всех параметров
	//TAU /= TIME_;
	//TMAX /= TIME_;

	// параметры обезразмеривания
	L_ = 1.0;
	R_ = 1.0;					// характерная плотность = начальная плотность
	P_ = 1.0;					// характерное давление = начальное давление 
	T_ = 1.0;					// характерная температура = начальная температура
	U_ = 1.0;		// характерная скорость = sqrt( P_ / R_ )
	E_ = 1.0;			// характерная энергия  = U_**2
	CV_ = 1.0;	// характерная теплоёмкость  = U_**2 / T_
	TIME_ = 1.0;			// характерное время
	MU_ = 1.0;		// характерная вязкость = R_ * U_ * L_
	KP_ = 1.0;	// коэффициент теплопроводности = R_ * U_**3 * L_ / T_
	CV_ = 1.0;	// характерная теплоёмкость  = U_**2 / T_


	// Обезразмеривание всех параметров
	TAU /= TIME_;
	TMAX /= TIME_;

	for (int i = 0; i < regCount; i++) {
		Region &reg = regions[i];
		Param &par = reg.par;

		Material& mat = materials[reg.matId];
		mat.URS(reg.par, 2);	// r=r(p,T)
		mat.URS(reg.par, 1);	// e=e(p,r)

		par.p /= P_;
		par.u /= U_;
		par.v /= U_;
		par.T /= T_;

		par.r /= R_;
		par.e /= E_;
		par.cz /= U_;
		par.ML /= MU_;
		par.E = par.e + par.U2()*0.5;

		int zhrv = 0;
	}

	Material::gR *= R_ * T_ / P_;	// Газовая постоянная
	for (int i = 0; i < matCount; i++)
	{
		Material & mat = materials[i];
		mat.Cp /= CV_;
		mat.K /= KP_;
		mat.ML /= MU_;
	}



	/* Чтение параметров ГУ */
	node0 = task->FirstChild("boundaries");
	TiXmlNode* bNode = node0->FirstChild("boundCond");
	while (bNode != NULL)
	{
		int edgeType;
		bNode->ToElement()->Attribute("edgeType", &edgeType);

		CFDBoundary * b;

		try {
			b = CFDBoundary::create(bNode, &grid);
		}
		catch (Exception e) {
			log("ERROR: %s\n", e.getMessage());
			exit(e.getType());
		}

		// TODO: !!!! костыль для обезразмеривания
		if (b->parCount == 4) {
			b->par[0] /= U_;
			b->par[1] /= U_;
			b->par[2] /= T_;
			b->par[3] /= P_;
		}

		boundaries.push_back(b);

		bNode = bNode->NextSibling("boundCond");
	}

	bCount = boundaries.size();

	
	/* Чтение данных сетки. */
	//node0 = task->FirstChild("mesh");
	//const char* fName = node0->FirstChild("name")->ToElement()->Attribute("value");
	//const char* tName = node0->FirstChild("filesType")->ToElement()->Attribute("value");
	//MeshReader* mr = MeshReader::create(MeshReader::getType((char*)tName), (char*)fName);
	//mr->read(&grid);

	grid.readMeshFiles();
	//grid.saveMeshInfo();

	/* Определение ГУ для каждого ребра. */
	for (int iEdge = 0; iEdge < grid.eCount; iEdge++) {
		Edge & e = grid.edges[iEdge];
		if (e.type == Edge::TYPE_INNER) {
			e.bnd = NULL;
			continue;
		}
		if (e.type == Edge::TYPE_NAMED) {
			int iBound = -1;
			for (int i = 0; i < bCount; i++)
			{
				if (strcmp(e.typeName, boundaries[i]->name) == 0)
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

			e.bnd = boundaries[iBound];
		}
		else {
			int iBound = -1;
			for (int i = 0; i < bCount; i++)
			{
				if (e.type == boundaries[i]->edgeType)
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

			e.bnd = boundaries[iBound];
		}
	}



	Parallel::buf = new double[10*grid.cCountEx];

	cTau = new double[grid.cCountEx];

	ro = new VECTOR[grid.cCountEx];
	ru = new VECTOR[grid.cCountEx];
	rv = new VECTOR[grid.cCountEx];
	re = new VECTOR[grid.cCountEx];

	ro_old = new VECTOR[grid.cCountEx];
	ru_old = new VECTOR[grid.cCountEx];
	rv_old = new VECTOR[grid.cCountEx];
	re_old = new VECTOR[grid.cCountEx];

	ro_int = new VECTOR[grid.cCountEx];
	ru_int = new VECTOR[grid.cCountEx];
	rv_int = new VECTOR[grid.cCountEx];
	re_int = new VECTOR[grid.cCountEx];

	for (int i = 0; i < grid.cCountEx; i++)
	{
		ro[i].init(FUNC_COUNT);
		ru[i].init(FUNC_COUNT);
		rv[i].init(FUNC_COUNT);
		re[i].init(FUNC_COUNT);

		ro_old[i].init(FUNC_COUNT);
		ru_old[i].init(FUNC_COUNT);
		rv_old[i].init(FUNC_COUNT);
		re_old[i].init(FUNC_COUNT);

		ro_int[i].init(FUNC_COUNT);
		ru_int[i].init(FUNC_COUNT);
		rv_int[i].init(FUNC_COUNT);
		re_int[i].init(FUNC_COUNT);
	}

	
	cellJ = new double[grid.cCountEx];
	cellGW = new double*[grid.cCountEx];
	cellGP = new Point*[grid.cCountEx];
	for (int i = 0; i < grid.cCountEx; i++) {
		cellGW[i] = new double[GP_CELL_COUNT];
		cellGP[i] = new Point[GP_CELL_COUNT];
	}

	edgeJ = new double[grid.eCountEx];
	edgeGW = new double*[grid.eCountEx];
	edgeGP = new Point*[grid.eCountEx];
	for (int i = 0; i < grid.eCountEx; i++) {
		edgeGW[i] = new double[GP_EDGE_COUNT];
		edgeGP[i] = new Point[GP_EDGE_COUNT];
	}

	matrA = new double**[grid.cCount];
	matrInvA = new double**[grid.cCount];

	for (int i = 0; i < grid.cCount; i++) {
		matrA[i] = new double*[FUNC_COUNT];
		matrInvA[i] = new double*[FUNC_COUNT];
		for (int j = 0; j < FUNC_COUNT; j++) {
			matrA[i][j] = new double[FUNC_COUNT];
			matrInvA[i][j] = new double[FUNC_COUNT];
		}
	}

	dBuf = new double[grid.cCountEx];
	iBuf = new int[grid.cCountEx];
	vBuf = new VECTOR[grid.cCountEx];


	calcGaussPar();

	calcMassMatr();
	for (int i = 0; i < grid.cCount; i++)
	{
		Cell & c = grid.cells[i];
		Region & reg = getRegion(c.typeName);
		convertParToCons(i, reg.par);
	}

	//int* test = new int[grid.cCountEx];
	//for (int i = 0; i < grid.cCountEx; i++) {
	//	test[i] = i;
	//}
	//exchange(test);

	procExchangeData();

	memcpy(ro_old, ro, grid.cCountEx*sizeof(double));
	memcpy(ru_old, ru, grid.cCountEx*sizeof(double));
	memcpy(rv_old, rv, grid.cCountEx*sizeof(double));
	memcpy(re_old, re, grid.cCountEx*sizeof(double));

	
	
	calcTimeStep();
	save(0);

}


void FEM_RKDG::calcTimeStep()
{
	if (!STEADY) {
		double tau = 1.0e+20;
		for (int iCell = 0; iCell < grid.cCount; iCell++)
		{
			Param p;
			double U, maxU = 0.0;
			for (int i = 0; i < GP_CELL_COUNT; i++) {
				convertConsToPar(iCell, cellGP[iCell][i], p);
				U = p.magU() + p.cz; //_max_(abs(p.u) + p.cz, abs(p.v) + p.cz);
				if (maxU < U) maxU = U;
			}
			for (int i = 0; i < 3; i++) {
				int iEdge = grid.cells[iCell].edgesInd[i];
				for (int iG = 0; iG < GP_EDGE_COUNT; iG++) {
					convertConsToPar(iCell, edgeGP[iEdge][i], p);
					U = p.magU() + p.cz; //_max_(abs(p.u) + p.cz, abs(p.v) + p.cz);
					if (maxU < U) maxU = U;
				}
			}

			convertConsToPar(iCell, grid.cells[iCell].c, p);
			U = p.magU() + p.cz; //_max_(abs(p.u) + p.cz, abs(p.v) + p.cz);
			if (maxU < U) maxU = U;

			double tmp = grid.cells[iCell].S / maxU;
			if (tmp < tau) tau = tmp;
		}
		tau *= CFL;
		TAU = _min_(TAU, tau);
		MPI_Allreduce(&TAU, &tau, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
		TAU = tau;
		for (int i = 0; i < grid.cCount; i++) {
			cTau[i] = TAU;
		}
		printf("\n\nUNSTEADY. Proc = %4d. Time step TAU = %e.\n\n", Parallel::procId, TAU);
		fflush(stdout);
	}
	else {
		for (int iCell = 0; iCell < grid.cCount; iCell++)
		{
			Param p;
			double U, maxU = 0.0;
			for (int i = 0; i < GP_CELL_COUNT; i++) {
				convertConsToPar(iCell, cellGP[iCell][i], p);
				U = p.magU() + p.cz; //_max_(abs(p.u) + p.cz, abs(p.v) + p.cz);
				if (maxU < U) maxU = U;
			}
			for (int i = 0; i < 3; i++) {
				int iEdge = grid.cells[iCell].edgesInd[i];
				for (int iG = 0; iG < GP_EDGE_COUNT; iG++) {
					convertConsToPar(iCell, edgeGP[iEdge][i], p);
					U = p.magU() + p.cz; //_max_(abs(p.u) + p.cz, abs(p.v) + p.cz);
					if (maxU < U) maxU = U;
				}
			}

			convertConsToPar(iCell, grid.cells[iCell].c, p);
			U = p.magU() + p.cz; //_max_(abs(p.u) + p.cz, abs(p.v) + p.cz);
			if (maxU < U) maxU = U;

			cTau[iCell] = CFL * grid.cells[iCell].S / maxU;
		}
	}
}

void FEM_RKDG::calcMassMatr()
{
	for (int iCell = 0; iCell < grid.cCount; iCell++) {
		double **	A		= matrA[iCell];
		double **	invA	= matrInvA[iCell];
		for (int i = 0; i < FUNC_COUNT; i++) {
			for (int j = 0; j < FUNC_COUNT; j++) {
				A[i][j] = 0.0;
				for (int iGP = 0; iGP < GP_CELL_COUNT; iGP++) {
					A[i][j] += cellGW[iCell][iGP] * getF(i, iCell, cellGP[iCell][iGP])
						* getF(j, iCell, cellGP[iCell][iGP]);
				}
				A[i][j] *= cellJ[iCell];
				//A[i][j] = _NORM_(A[i][j]);
			}
		}
		inverseMatr(A, invA, FUNC_COUNT);
	}
}


void FEM_RKDG::calcGaussPar()
{

	// для ячеек
	if (GP_CELL_COUNT == 4) {
		for (int i = 0; i < grid.cCount; i++) {

			double a = 1.0 / 3.0;
			double b = 1.0 / 5.0;
			double c = 3.0 / 5.0;
			double x1 = grid.nodes[grid.cells[i].nodesInd[0]].x;
			double y1 = grid.nodes[grid.cells[i].nodesInd[0]].y;
			double x2 = grid.nodes[grid.cells[i].nodesInd[1]].x;
			double y2 = grid.nodes[grid.cells[i].nodesInd[1]].y;
			double x3 = grid.nodes[grid.cells[i].nodesInd[2]].x;
			double y3 = grid.nodes[grid.cells[i].nodesInd[2]].y;
			double a1 = x1 - x3;
			double a2 = y1 - y3;
			double b1 = x2 - x3;
			double b2 = y2 - y3;
			double c1 = x3;
			double c2 = y3;

			cellGW[i][0] = -27.0 / 48.0;
			cellGP[i][0].x = a1*a + b1*a + c1;
			cellGP[i][0].y = a2*a + b2*a + c2;

			cellGW[i][1] = 25.0 / 48.0;
			cellGP[i][1].x = a1*c + b1*b + c1;
			cellGP[i][1].y = a2*c + b2*b + c2;

			cellGW[i][2] = 25.0 / 48.0;
			cellGP[i][2].x = a1*b + b1*c + c1;
			cellGP[i][2].y = a2*b + b2*c + c2;

			cellGW[i][3] = 25.0 / 48.0;
			cellGP[i][3].x = a1*b + b1*b + c1;
			cellGP[i][3].y = a2*b + b2*b + c2;

			cellJ[i] = 0.5*fabs(a1*b2 - a2*b1);
			//cellJ[i] = 0.5*(a1*b2 - a2*b1);
			if (cellJ[i] <= 0) {
				int zhrv = 0;
			}
		}
	}
	else if (GP_CELL_COUNT == 3) {
		for (int i = 0; i < grid.cCount; i++) {
			double a = 1.0 / 6.0;
			double b = 2.0 / 3.0;
			double x1 = grid.nodes[grid.cells[i].nodesInd[0]].x;
			double y1 = grid.nodes[grid.cells[i].nodesInd[0]].y;
			double x2 = grid.nodes[grid.cells[i].nodesInd[1]].x;
			double y2 = grid.nodes[grid.cells[i].nodesInd[1]].y;
			double x3 = grid.nodes[grid.cells[i].nodesInd[2]].x;
			double y3 = grid.nodes[grid.cells[i].nodesInd[2]].y;
			double a1 = x1 - x3;
			double a2 = y1 - y3;
			double b1 = x2 - x3;
			double b2 = y2 - y3;
			double c1 = x3;
			double c2 = y3;

			cellGW[i][0] = 1.0 / 6.0;
			cellGP[i][0].x = a1*a + b1*a + c1;
			cellGP[i][0].y = a2*a + b2*a + c2;

			cellGW[i][1] = 1.0 / 6.0;
			cellGP[i][1].x = a1*a + b1*b + c1;
			cellGP[i][1].y = a2*a + b2*b + c2;

			cellGW[i][2] = 1.0 / 6.0;
			cellGP[i][2].x = a1*b + b1*a + c1;
			cellGP[i][2].y = a2*b + b2*a + c2;

			cellJ[i] = fabs(a1*b2 - a2*b1);

		}
	}

	// для ребер
	if (GP_EDGE_COUNT == 3) {
		for (int i = 0; i < grid.eCount; i++) {
			double gp1 = -3.0 / 5.0;
			double gp2 = 0.0;
			double gp3 = 3.0 / 5.0;
			double x1 = grid.nodes[grid.edges[i].n1].x;
			double y1 = grid.nodes[grid.edges[i].n1].y;
			double x2 = grid.nodes[grid.edges[i].n2].x;
			double y2 = grid.nodes[grid.edges[i].n2].y;

			edgeGW[i][0] = 5.0 / 9.0;
			edgeGP[i][0].x = ((x1 + x2) + gp1*(x2 - x1)) / 2.0;
			edgeGP[i][0].y = ((y1 + y2) + gp1*(y2 - y1)) / 2.0;

			edgeGW[i][1] = 8.0 / 9.0;
			edgeGP[i][1].x = (x1 + x2) / 2.0;
			edgeGP[i][1].y = (y1 + y2) / 2.0;

			edgeGW[i][2] = 5.0 / 9.0;
			edgeGP[i][2].x = ((x1 + x2) + gp3*(x2 - x1)) / 2.0;
			edgeGP[i][2].y = ((y1 + y2) + gp3*(y2 - y1)) / 2.0;

			edgeJ[i] = sqrt(POW_2(x2 - x1) + POW_2(y2 - y1))*0.5;

		}
	}
	else if (GP_EDGE_COUNT == 2) {
		for (int i = 0; i < grid.eCount; i++) {
			double _sqrt3 = 1.0 / sqrt(3.0);
			double x1 = grid.nodes[grid.edges[i].n1].x;
			double y1 = grid.nodes[grid.edges[i].n1].y;
			double x2 = grid.nodes[grid.edges[i].n2].x;
			double y2 = grid.nodes[grid.edges[i].n2].y;

			edgeGW[i][0] = 1.0;
			edgeGP[i][0].x = ((x1 + x2) - _sqrt3*(x2 - x1))*0.5;
			edgeGP[i][0].y = ((y1 + y2) - _sqrt3*(y2 - y1))*0.5;

			edgeGW[i][1] = 1.0;
			edgeGP[i][1].x = ((x1 + x2) + _sqrt3*(x2 - x1))*0.5;
			edgeGP[i][1].y = ((y1 + y2) + _sqrt3*(y2 - y1))*0.5;

			edgeJ[i] = sqrt(POW_2(x2 - x1) + POW_2(y2 - y1))*0.5;

		}
	}
}


void FEM_RKDG::run() 
{


	int nc = grid.cCount;
	int ne = grid.eCount;

	double			t		= 0.0;
	unsigned int	step	= 0;
	while (t < TMAX) 
	{
		if (!STEADY) {
			t += TAU;
		}
		else {
			calcTimeStep();
		}
		step++;
		copyToOld();

		{	// RK-step 1
			zeroIntegrals();
			calcLimiters();
			procExchangeData();
			calcConvectionVol();
			calcConvectionSurf();
			calcNewFields();
			procExchangeData();
		}

		{	// RK-step 2
			zeroIntegrals();
			calcLimiters();
			procExchangeData();
			calcConvectionVol();
			calcConvectionSurf();
			calcNewFields2();
			procExchangeData();
		}

		calcResiduals();

		Parallel::barrier();

		if (step % FILE_SAVE_STEP == 0)
		{
			save(step);
		}
		if (step % PRINT_STEP == 0)
		{
			if (!STEADY) {
				log("step: %d\t\ttime step: %.16f\n", step, t);
			}
			else {
				log("step: %6d\tres: ro:%10.5e ru:%10.5e rv:%10.5e re:%10.5e\n", step, ro_res, ru_res, rv_res, re_res); 
			}
		}
		
	}

}

void FEM_RKDG::save(int step)
{
	char fName[50];
	FILE * fp;
	if (Parallel::isRoot()) {
		sprintf(fName, "res_%010d.pvtu", step);
		fp = fopen(fName, "w");
		fprintf(fp, "<?xml version=\"1.0\"?>\n");
		fprintf(fp, "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
		fprintf(fp, "  <PUnstructuredGrid GhostLevel=\"0\">\n");
		fprintf(fp, "    <PPoints><PDataArray type=\"Float32\" NumberOfComponents = \"3\" /></PPoints>\n");
		fprintf(fp, "    <PCellData Scalars=\"Density, Pressure, Proc\" Vectors=\"Velocity\">\n");
		fprintf(fp, "      <PDataArray type=\"Float32\" Name=\"Density\" format=\"ascii\" />\n");
		fprintf(fp, "      <PDataArray type=\"Float32\" Name=\"Pressure\" format=\"ascii\" />\n");
		fprintf(fp, "      <PDataArray type=\"Float32\" Name=\"Velocity\" format=\"ascii\" NumberOfComponents=\"3\"/>\n");
		fprintf(fp, "      <PDataArray type=\"Int32\" Name=\"Proc\" format=\"ascii\" />\n");
		fprintf(fp, "      <PDataArray type=\"Int32\" Name=\"Loc_ID\" format=\"ascii\" />\n");
		fprintf(fp, "    </PCellData>\n");
		
		
		for (int pid = 0; pid < Parallel::procCount; pid++) {
			fprintf(fp, "    <Piece Source=\"vtk_data/res_%010d.%04d.vtu\"/>\n", step, pid);
		}
		fprintf(fp, "  </PUnstructuredGrid>\n");

		fprintf(fp, "</VTKFile>\n");
		fclose(fp);
		log("File '%s' saved...\n", fName);
	}

	sprintf(fName, "vtk_data/res_%010d.%04d.vtu", step, Parallel::procId);
	log("File '%s' saved...\n", fName);
	fp = fopen(fName, "w");
	fprintf(fp, "<?xml version=\"1.0\"?>\n");
	fprintf(fp, "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
	fprintf(fp, "  <UnstructuredGrid GhostLevel=\"0\">\n");
	fprintf(fp, "    <Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", grid.nCount, grid.cCount);
	fprintf(fp, "      <Points>\n");
	fprintf(fp, "        <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n");
	fprintf(fp, "          ");
	for (int i = 0; i < grid.nCount; i++) {
		fprintf(fp, "%25.15f %25.15f %25.15f ", grid.nodes[i].x, grid.nodes[i].y, 0.0);
	}
	fprintf(fp, "\n");
	fprintf(fp, "        </DataArray>\n");
	fprintf(fp, "      </Points>\n");
	fprintf(fp, "      <Cells>\n");
	fprintf(fp, "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n");
	fprintf(fp, "          ");
	for (int i = 0; i < grid.cCount; i++) {
		fprintf(fp, "%d ", (i+1)*3);
	}
	fprintf(fp, "\n");
	fprintf(fp, "        </DataArray>\n");
	fprintf(fp, "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n");
	fprintf(fp, "          ");
	for (int i = 0; i < grid.cCount; i++) {
		fprintf(fp, "%d %d %d ", grid.cells[i].nodesInd[0], grid.cells[i].nodesInd[1], grid.cells[i].nodesInd[2]);
	}
	fprintf(fp, "\n");
	fprintf(fp, "        </DataArray>\n");
	fprintf(fp, "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n");
	fprintf(fp, "          ");
	for (int i = 0; i < grid.cCount; i++) {
		fprintf(fp, "%d ", 5);
	}
	fprintf(fp, "\n");
	fprintf(fp, "        </DataArray>\n");
	fprintf(fp, "      </Cells>\n");
	
	
	fprintf(fp, "      <CellData Scalars=\"Density, Pressure, Proc\" Vectors=\"Velocity\">\n");
	
	// плотность
	fprintf(fp, "        <DataArray type=\"Float32\" Name=\"Density\" format=\"ascii\">\n");
	fprintf(fp, "          ");
	for (int i = 0; i < grid.cCount; i++) {
		Param p;
		convertConsToPar(i, p);
		fprintf(fp, "%25.16f ", p.r*R_);
		if (i + 1 % 8 == 0 || i + 1 == grid.cCount) fprintf(fp, "\n");
	}
	fprintf(fp, "        </DataArray>\n");

	// давление
	fprintf(fp, "        <DataArray type=\"Float32\" Name=\"Pressure\" format=\"ascii\">\n");
	fprintf(fp, "          ");
	for (int i = 0; i < grid.cCount; i++) {
		Param p;
		convertConsToPar(i, p);
		fprintf(fp, "%25.16f ", p.p*P_);
		if (i + 1 % 8 == 0 || i + 1 == grid.cCount) fprintf(fp, "\n");
	}
	fprintf(fp, "        </DataArray>\n");
	
	// скорость
	fprintf(fp, "        <DataArray type=\"Float32\" Name=\"Velocity\" format=\"ascii\" NumberOfComponents=\"3\">\n");
	fprintf(fp, "          ");
	for (int i = 0; i < grid.cCount; i++) {
		Param p;
		convertConsToPar(i, p);
		fprintf(fp, "%25.16f %25.16f %25.16f ", p.u*U_, p.v*U_, 0.0);
		if (i + 1 % 8 == 0 || i + 1 == grid.cCount) fprintf(fp, "\n");
	}
	fprintf(fp, "        </DataArray>\n");
	
	// номер процессора
	fprintf(fp, "        <DataArray type=\"Int32\" Name=\"Proc\" format=\"ascii\">\n");
	fprintf(fp, "          ");
	for (int i = 0; i < grid.cCount; i++) {
		fprintf(fp, "%d ", Parallel::procId);
		if (i + 1 % 8 == 0 || i + 1 == grid.cCount) fprintf(fp, "\n");
	}
	fprintf(fp, "        </DataArray>\n");

	// локальный номер ячейки
	fprintf(fp, "        <DataArray type=\"Int32\" Name=\"Loc_ID\" format=\"ascii\">\n");
	fprintf(fp, "          ");
	for (int i = 0; i < grid.cCount; i++) {
		fprintf(fp, "%d ", i);
		if (i + 1 % 8 == 0 || i + 1 == grid.cCount) fprintf(fp, "\n");
	}
	fprintf(fp, "        </DataArray>\n");

	fprintf(fp, "      </CellData>\n");
	
	
	fprintf(fp, "    </Piece>\n");


	fprintf(fp, "  </UnstructuredGrid>\n");

	fprintf(fp, "</VTKFile>\n");
	fclose(fp);
	return;
	
	sprintf(fName, "vtk_debug/res_%010d.%04d.vtk", step, Parallel::procId);
	fp = fopen(fName, "w");
	fprintf(fp, "# vtk DataFile Version 2.0\n");
	fprintf(fp, "GASDIN data file\n");
	fprintf(fp, "ASCII\n");
	fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");

	fprintf(fp, "POINTS %d float\n", grid.nCountEx);
	for (int i = 0; i < grid.nCountEx; i++)
	{
		fprintf(fp, "%f %f %f  ", grid.nodes[i].x, grid.nodes[i].y, 0.0);
		if (i + 1 % 8 == 0) fprintf(fp, "\n");
	}
	fprintf(fp, "\n");
	fprintf(fp, "CELLS %d %d\n", grid.cCountEx, 4 * grid.cCountEx);
	for (int i = 0; i < grid.cCountEx; i++)
	{
		fprintf(fp, "3 %d %d %d\n", grid.cells[i].nodesInd[0], grid.cells[i].nodesInd[1], grid.cells[i].nodesInd[2]);
	}
	fprintf(fp, "\n");

	fprintf(fp, "CELL_TYPES %d\n", grid.cCountEx);
	for (int i = 0; i < grid.cCountEx; i++) fprintf(fp, "5\n");
	fprintf(fp, "\n");

	fprintf(fp, "CELL_DATA %d\nSCALARS Density float 1\nLOOKUP_TABLE default\n", grid.cCountEx);
	for (int i = 0; i < grid.cCountEx; i++)
	{
		Param p;
		convertConsToPar(i, p);
		fprintf(fp, "%25.16f ", p.r*R_);
		if (i + 1 % 8 == 0 || i + 1 == grid.cCountEx) fprintf(fp, "\n");
	}

	fprintf(fp, "SCALARS Pressure float 1\nLOOKUP_TABLE default\n", grid.cCountEx);
	for (int i = 0; i < grid.cCountEx; i++)
	{
		Param p;
		convertConsToPar(i, p);
		fprintf(fp, "%25.16f ", p.p*P_);
		if (i + 1 % 8 == 0 || i + 1 == grid.cCountEx) fprintf(fp, "\n");
	}

	fprintf(fp, "SCALARS Temperature float 1\nLOOKUP_TABLE default\n", grid.cCountEx);
	for (int i = 0; i < grid.cCountEx; i++)
	{
		Param p;
		convertConsToPar(i, p);
		fprintf(fp, "%25.16f ", p.T*T_);
		if (i + 1 % 8 == 0 || i + 1 == grid.cCountEx) fprintf(fp, "\n");
	}

	fprintf(fp, "SCALARS MachNumber float 1\nLOOKUP_TABLE default\n", grid.cCountEx);
	for (int i = 0; i < grid.cCountEx; i++)
	{
		Param p;
		convertConsToPar(i, p);
		fprintf(fp, "%25.16f ", sqrt(p.u*p.u + p.v*p.v) / p.cz);
		if (i + 1 % 8 == 0 || i + 1 == grid.cCountEx) fprintf(fp, "\n");
	}

	fprintf(fp, "VECTORS Velosity float\n");
	for (int i = 0; i < grid.cCountEx; i++)
	{
		Param p;
		convertConsToPar(i, p);
		fprintf(fp, "%25.16f %25.16f %25.16f ", p.u*U_, p.v*U_, 0.0);
		if (i + 1 % 8 == 0 || i + 1 == grid.cCountEx) fprintf(fp, "\n");
	}

	fprintf(fp, "SCALARS Total_energy float 1\nLOOKUP_TABLE default\n");
	for (int i = 0; i < grid.cCountEx; i++)
	{
		Param p;
		convertConsToPar(i, p);
		fprintf(fp, "%25.16f ", p.E*E_);
		if (i + 1 % 8 == 0 || i + 1 == grid.cCountEx) fprintf(fp, "\n");
	}

	fprintf(fp, "SCALARS Total_pressure float 1\nLOOKUP_TABLE default\n");
	for (int i = 0; i < grid.cCountEx; i++)
	{
		Material &mat = getMaterial(i);
		double gam = mat.getGamma();
		double agam = gam - 1.0;
		Param p;
		convertConsToPar(i, p);
		double M2 = (p.u*p.u + p.v*p.v) / (p.cz*p.cz);
		fprintf(fp, "%25.16e ", P_*p.p*::pow(1.0 + 0.5*M2*agam, gam / agam));
		if ((i + 1) % 8 == 0 || i + 1 == grid.cCountEx) fprintf(fp, "\n");
	}

	fprintf(fp, "SCALARS TAU float 1\nLOOKUP_TABLE default\n", grid.cCountEx);
	for (int i = 0; i < grid.cCountEx; i++)
	{
		fprintf(fp, "%25.16f ", cTau[i]*TIME_);
		if (i + 1 % 8 == 0 || i + 1 == grid.cCountEx) fprintf(fp, "\n");
	}

	fprintf(fp, "SCALARS SQUARE float 1\nLOOKUP_TABLE default\n", grid.cCount);
	for (int i = 0; i < grid.cCountEx; i++)
	{
		fprintf(fp, "%25.16f ", grid.cells[i].S*L_*L_);
		if (i + 1 % 8 == 0 || i + 1 == grid.cCountEx) fprintf(fp, "\n");
	}

	fclose(fp);
	
	
	return;
}

void FEM_RKDG::procExchangeData()
{
	exchange(ro);
	exchange(ru);
	exchange(rv);
	exchange(re);
}

void FEM_RKDG::calcFlux(double& fr, double& fu, double& fv, double& fe, Param pL, Param pR, Vector n, double GAM)
{
	if (FLUX == FLUX_GODUNOV) {	// GODUNOV FLUX
		double RI, EI, PI, UI, VI, WI, UN, UT;
		double unl = pL.u*n.x + pL.v*n.y;
		double unr = pR.u*n.x + pR.v*n.y;
		double utl = pL.u*n.y - pL.v*n.x;
		double utr = pR.u*n.y - pR.v*n.x;
		rim_orig(RI, EI, PI, UN, UT, WI, pL.r, pL.p, unl, utl, 0, pR.r, pR.p, unr, utr, 0, GAM);

		UI = UN*n.x + UT*n.y;
		VI = UN*n.y - UT*n.x;

		fr = RI*UN;
		fu = fr*UI + PI*n.x;
		fv = fr*VI + PI*n.y;
		fe = (RI*(EI + 0.5*(UI*UI + VI*VI)) + PI)*UN;
	}
	if (FLUX == FLUX_LAX) {	// LAX-FRIEDRIX FLUX
		double unl = pL.u*n.x + pL.v*n.y;
		double unr = pR.u*n.x + pR.v*n.y;
		double rol, rul, rvl, rel, ror, rur, rvr, rer;
		double alpha = _max_(fabs(unl) + sqrt(GAM*pL.p / pL.r), fabs(unr) + sqrt(GAM*pR.p / pR.r));
		//pL.getToCons(rol, rul, rvl, rel);
		//pR.getToCons(ror, rur, rvr, rer);
		rol = pL.r;
		rul = pL.r*pL.u;
		rvl = pL.r*pL.v;
		rel = pL.p / (GAM - 1.0) + 0.5*pL.r*(pL.U2());
		ror = pR.r;
		rur = pR.r*pR.u;
		rvr = pR.r*pR.v;
		rer = pR.p / (GAM - 1.0) + 0.5*pR.r*(pR.U2());
		double frl = rol*unl;
		double frr = ror*unr;
		fr = 0.5*(frr + frl - alpha*(ror - rol));
		fu = 0.5*(frr*pR.u + frl*pL.u + (pR.p + pL.p)*n.x - alpha*(rur - rul));
		fv = 0.5*(frr*pR.v + frl*pL.v + (pR.p + pL.p)*n.y - alpha*(rvr - rvl));
		fe = 0.5*((rer + pR.p)*unr + (rel + pL.p)*unl - alpha*(rer - rel));
	}
}



void FEM_RKDG::zeroIntegrals()
{
	for (int i = 0; i < grid.cCount; i++) {
		for (int iF = 0; iF < FUNC_COUNT; iF++) {
			ro_int[i][iF] = 0.0;
			ru_int[i][iF] = 0.0;
			rv_int[i][iF] = 0.0;
			re_int[i][iF] = 0.0;

			ro_old[i][iF] = ro[i][iF];
			ru_old[i][iF] = ru[i][iF];
			rv_old[i][iF] = rv[i][iF];
			re_old[i][iF] = re[i][iF];
		}

	}
}


void FEM_RKDG::copyToOld()
{
	for (int i = 0; i < grid.cCount; i++) {
		for (int iF = 0; iF < FUNC_COUNT; iF++) {
			ro_old[i][iF] = ro[i][iF];
			ru_old[i][iF] = ru[i][iF];
			rv_old[i][iF] = rv[i][iF];
			re_old[i][iF] = re[i][iF];
		}

	}
}


void FEM_RKDG::calcConvectionSurf()
{
	double fr, fu, fv, fe;
	for (int i = 0; i < grid.eCount; i++)
	{
		Edge& edge = grid.edges[i];
		Vector n = edge.n;
		if (edge.type == Edge::TYPE_INNER)
		{

			int c1 = edge.c1;
			int c2 = edge.c2;
			VECTOR intRO1(FUNC_COUNT), intRU1(FUNC_COUNT), intRV1(FUNC_COUNT), intRE1(FUNC_COUNT);
			VECTOR intRO2(FUNC_COUNT), intRU2(FUNC_COUNT), intRV2(FUNC_COUNT), intRE2(FUNC_COUNT);
			for (int iGP = 0; iGP < GP_EDGE_COUNT; iGP++)
			{
				Point& pt = edgeGP[i][iGP];
				Param pL, pR;
				convertConsToPar(c1, pt, pL);
				convertConsToPar(c2, pt, pR);
				Material &mat = getMaterial(c1);
				double __GAM = mat.getGamma();   // TODO: сделать правильное вычисление показателя адиабаты
				calcFlux(fr, fu, fv, fe, pL, pR, n, __GAM);
				
				VECTOR tmp(FUNC_COUNT), tmp1(FUNC_COUNT), tmp2(FUNC_COUNT);
				tmp1     = getF(c1, pt);
				tmp1	*= edgeGW[i][iGP];
				tmp2     = getF(c2, pt);
				tmp2	*= edgeGW[i][iGP];
				
				tmp     = tmp1;
				tmp    *= fr;
				intRO1 += tmp;
				tmp     = tmp1;
				tmp    *= fu;
				intRU1 += tmp;
				tmp     = tmp1;
				tmp    *= fv;
				intRV1 += tmp;
				tmp     = tmp1;
				tmp    *= fe;
				intRE1 += tmp;



				tmp     = tmp2;
				tmp    *= fr;
				intRO2 += tmp;
				tmp     = tmp2;
				tmp    *= fu;
				intRU2 += tmp;
				tmp     = tmp2;
				tmp    *= fv;
				intRV2 += tmp;
				tmp     = tmp2;
				tmp    *= fe;
				intRE2 += tmp;
			}
			
			intRO1 *= edgeJ[i];
			intRU1 *= edgeJ[i];
			intRV1 *= edgeJ[i];
			intRE1 *= edgeJ[i];
			
			intRO2 *= edgeJ[i];
			intRU2 *= edgeJ[i];
			intRV2 *= edgeJ[i];
			intRE2 *= edgeJ[i];

			ro_int[c1] -= intRO1;
			ru_int[c1] -= intRU1;
			rv_int[c1] -= intRV1;
			re_int[c1] -= intRE1;

			ro_int[c2] += intRO2;
			ru_int[c2] += intRU2;
			rv_int[c2] += intRV2;
			re_int[c2] += intRE2;

		} else {

			int c1 = edge.c1;
			VECTOR intRO1(FUNC_COUNT), intRU1(FUNC_COUNT), intRV1(FUNC_COUNT), intRE1(FUNC_COUNT);
			VECTOR intRO2(FUNC_COUNT), intRU2(FUNC_COUNT), intRV2(FUNC_COUNT), intRE2(FUNC_COUNT);
			for (int iGP = 0; iGP < GP_EDGE_COUNT; iGP++)
			{
				Point& pt = edgeGP[i][iGP];
				Param pL, pR;
				convertConsToPar(c1, pt, pL);
				boundaryCond(i, pt, pL, pR);
				Material &mat = getMaterial(c1);
				double __GAM = mat.getGamma();
				calcFlux(fr, fu, fv, fe, pL, pR, n, __GAM);

				VECTOR tmp(FUNC_COUNT), tmp1(FUNC_COUNT);
				tmp1 = getF(c1, pt);
				tmp1 *= edgeGW[i][iGP];

				tmp = tmp1;
				tmp *= fr;
				intRO1 += tmp;

				tmp = tmp1;
				tmp *= fu;
				intRU1 += tmp;

				tmp = tmp1;
				tmp *= fv;
				intRV1 += tmp;

				tmp = tmp1;
				tmp *= fe;
				intRE1 += tmp;
			}
			intRO1 *= edgeJ[i];
			intRU1 *= edgeJ[i];
			intRV1 *= edgeJ[i];
			intRE1 *= edgeJ[i];
			
			ro_int[c1] -= intRO1;
			ru_int[c1] -= intRU1;
			rv_int[c1] -= intRV1;
			re_int[c1] -= intRE1;

		}
	}
}


void FEM_RKDG::calcConvectionVol()
{
	for (int i = 0; i < grid.cCount; i++)
	{
		VECTOR intRO(FUNC_COUNT), intRU(FUNC_COUNT), intRV(FUNC_COUNT), intRE(FUNC_COUNT);
		for (int iGP = 0; iGP < GP_CELL_COUNT; iGP++)
		{
			Point& pt = cellGP[i][iGP];
			Param par;
			convertConsToPar(i, pt, par);
			
			double fr = par.r*par.u;
			double fu = fr*par.u+par.p;
			double fv = fr*par.v;
			double fe = par.u*(par.r*par.E+par.p);
			
			double gr = par.r*par.v;
			double gu = gr*par.u;
			double gv = gr*par.v+par.p;
			double ge = par.v*(par.r*par.E+par.p);
			
			VECTOR tmpInt(FUNC_COUNT), tmpInt1(FUNC_COUNT);
			
			
			// RO
			tmpInt = 0.0;
			tmpInt1  = getDFDX(i, pt); 
			tmpInt1 *= fr;
			tmpInt += tmpInt1;
			
			tmpInt1  = getDFDY(i, pt); 
			tmpInt1 *= gr;
			tmpInt += tmpInt1;
			
			tmpInt *= cellGW[i][iGP];
			intRO  += tmpInt;

			// RU
			tmpInt = 0.0;
			tmpInt1  = getDFDX(i, pt); 
			tmpInt1 *= fu;
			tmpInt += tmpInt1;

			tmpInt1  = getDFDY(i, pt); 
			tmpInt1 *= gu;
			tmpInt += tmpInt1;

			tmpInt *= cellGW[i][iGP];
			intRU  += tmpInt;

			// RV
			tmpInt = 0.0;
			tmpInt1  = getDFDX(i, pt); 
			tmpInt1 *= fv;
			tmpInt += tmpInt1;

			tmpInt1  = getDFDY(i, pt); 
			tmpInt1 *= gv;
			tmpInt += tmpInt1;

			tmpInt *= cellGW[i][iGP];
			intRV  += tmpInt;

			// RE
			tmpInt = 0.0;
			tmpInt1  = getDFDX(i, pt); 
			tmpInt1 *= fe;
			tmpInt += tmpInt1;

			tmpInt1  = getDFDY(i, pt); 
			tmpInt1 *= ge;
			tmpInt += tmpInt1;

			tmpInt *= cellGW[i][iGP];
			intRE  += tmpInt;
		}
		intRO *= cellJ[i];
		intRU *= cellJ[i];
		intRV *= cellJ[i];
		intRE *= cellJ[i];

		ro_int[i] += intRO;
		ru_int[i] += intRU;
		rv_int[i] += intRV;
		re_int[i] += intRE;
	}
}

void FEM_RKDG::calcDiffusionVol() {} //TODO:

void FEM_RKDG::calcDiffusionSurf() {} //TODO:

void FEM_RKDG::calcNewFields()
{
	VECTOR tmp(FUNC_COUNT);

	for (int i = 0; i < grid.cCount; i++)
	{
		tmp = ro_int[i];
		tmp *= matrInvA[i];
		tmp *= cTau[i];
		ro[i] += tmp;

		tmp = ru_int[i];
		tmp *= matrInvA[i];
		tmp *= cTau[i];
		ru[i] += tmp;

		tmp = rv_int[i];
		tmp *= matrInvA[i];
		tmp *= cTau[i];
		rv[i] += tmp;

		tmp = re_int[i];
		tmp *= matrInvA[i];
		tmp *= cTau[i];
		re[i] += tmp;
	}
}

void FEM_RKDG::calcNewFields2()
{
	VECTOR tmp(FUNC_COUNT);

	for (int i = 0; i < grid.cCount; i++)
	{
		tmp = ro_int[i];
		tmp *= matrInvA[i];
		tmp *= cTau[i];
		ro[i] += tmp;

		tmp = ru_int[i];
		tmp *= matrInvA[i];
		tmp *= cTau[i];
		ru[i] += tmp;

		tmp = rv_int[i];
		tmp *= matrInvA[i];
		tmp *= cTau[i];
		rv[i] += tmp;

		tmp = re_int[i];
		tmp *= matrInvA[i];
		tmp *= cTau[i];
		re[i] += tmp;
	}

	for (int i = 0; i < grid.cCount; i++) {
		ro[i] += ro_old[i];
		ro[i] *= 0.5;

		ru[i] += ru_old[i];
		ru[i] *= 0.5;

		rv[i] += rv_old[i];
		rv[i] *= 0.5;

		re[i] += re_old[i];
		re[i] *= 0.5;
	}
}

void FEM_RKDG::calcResiduals()
{
	VECTOR v;
	double res;
	ro_res = 0.0;
	ru_res = 0.0;
	rv_res = 0.0;
	re_res = 0.0;
	for (int i = 0; i < grid.cCount; i++) {
		v = ro[i];
		v -= ro_old[i];
		v.abs();
		res = v.norm();
		if (res > ro_res) ro_res = res;

		v = ru[i];
		v -= ru_old[i];
		v.abs();
		res = v.norm();
		if (res > ru_res) ru_res = res;

		v = rv[i];
		v -= rv_old[i];
		v.abs();
		res = v.norm();
		if (res > rv_res) rv_res = res;

		v = re[i];
		v -= re_old[i];
		v.abs();
		res = v.norm();
		if (res > re_res) re_res = res;
	}
}

void FEM_RKDG::calcLimiters() {}


void FEM_RKDG::boundaryCond(int iEdge, Point pt, Param& pL, Param& pR)
{
	Edge &edge = grid.edges[iEdge];
	int c1 = edge.c1;
	Material& m = getMaterial(c1);
	if (edge.bnd) {
		edge.bnd->run(iEdge, pL, pR);
		m.URS(pR, 2);
		m.URS(pR, 1);
		pR.E = pR.e + 0.5*(pR.U2());
		return;
	}
	else {
		char msg[128];
		sprintf(msg, "Not defined boundary condition for edge %d\n", iEdge);
		throw Exception(msg, Exception::TYPE_BOUND_UNKNOWN);
	}

}


void FEM_RKDG::done()
{
	delete[] ro;
	delete[] ru;
	delete[] rv;
	delete[] re;

	delete[] ro_old;
	delete[] ru_old;
	delete[] rv_old;
	delete[] re_old;

	delete[] ro_int;
	delete[] ru_int;
	delete[] rv_int;
	delete[] re_int;
}




Region & FEM_RKDG::getRegionByCellType(int type)
{
	for (int i = 0; i < regCount; i++)
	{
		if (regions[i].cellType == type) return regions[i];
	}
	log("ERROR: unknown cell type %d...\n", type);
	EXIT(1);
}


Region   &	FEM_RKDG::getRegion	(int iCell)
{
	return getRegionByCellType( grid.cells[iCell].type );
}

Region & FEM_RKDG::getRegionByName(char* name)
{
	for (int i = 0; i < regCount; i++)
	{
		if (strcmp(regions[i].name.c_str(), name) == 0) return regions[i];
	}
	log("ERROR: unknown cell name '%s'...\n", name);
	EXIT(1);
}

Region & FEM_RKDG::getRegion(char * name)
{
	return getRegionByName(name);
}

Material &	FEM_RKDG::getMaterial(int iCell)
{
	Region & reg = getRegion(grid.cells[iCell].typeName);
	return materials[reg.matId];
}


void FEM_RKDG::convertParToCons(int iCell, Param & par)
{
	ro[iCell][0] = par.r;
	ru[iCell][0] = par.r*par.u;
	rv[iCell][0] = par.r*par.v;
	re[iCell][0] = par.r*(par.e+0.5*(par.u*par.u+par.v*par.v));

	for (int i = 1; i < FUNC_COUNT; i++)
	{
		ro[iCell][i] = 0.0;
		ru[iCell][i] = 0.0;
		rv[iCell][i] = 0.0;
		re[iCell][i] = 0.0;
	}
}

void FEM_RKDG::convertConsToPar(int iCell, Point pt, Param & par)
{
	par.r = getRO(iCell, pt);
	par.u = getRU(iCell, pt)/par.r;
	par.v = getRV(iCell, pt)/par.r;
	par.E = getRE(iCell, pt)/par.r;
	par.e = par.E-0.5*(par.u*par.u+par.v*par.v);
	Material& mat = getMaterial(iCell);
	mat.URS(par, 0);
	mat.URS(par, 1);
}

void FEM_RKDG::convertConsToPar(int iCell, Param & par)
{
	convertConsToPar(iCell, grid.cells[iCell].c, par);
}





VECTOR& FEM_RKDG::getF(int iCell, Point pt)
{
	if (FUNC_COUNT == 3 )
	{
		VECTOR v(3);
		v[0] = 1.0;
		v[1] = (pt.x - grid.cells[iCell].c.x)/grid.cells[iCell].HX;
		v[2] = (pt.y - grid.cells[iCell].c.y)/grid.cells[iCell].HY;
		return v;
	}
	else if (FUNC_COUNT == 6)
	{
		VECTOR v(6);
		v[0] = 1.0;
		v[1] = (pt.x - grid.cells[iCell].c.x)/grid.cells[iCell].HX;
		v[2] = (pt.y - grid.cells[iCell].c.y)/grid.cells[iCell].HY;
		v[3] = (pt.x - grid.cells[iCell].c.x)*(pt.x - grid.cells[iCell].c.x)/(grid.cells[iCell].HX*grid.cells[iCell].HX);
		v[4] = (pt.y - grid.cells[iCell].c.y)*(pt.y - grid.cells[iCell].c.y)/(grid.cells[iCell].HY*grid.cells[iCell].HY);
		v[5] = (pt.x - grid.cells[iCell].c.x)*(pt.y - grid.cells[iCell].c.y)/(grid.cells[iCell].HX*grid.cells[iCell].HY);
		return v;
	} else {
		log("ERROR: wrong basic functions count!\n");
		EXIT(1);
	}
}

VECTOR FEM_RKDG::getDFDX(int iCell, Point pt)
{
	if (FUNC_COUNT == 3 )
	{
		VECTOR v(3);
		v[0] = 0.0;
		v[1] = 1.0/grid.cells[iCell].HX;
		v[2] = 0.0;
		return v;
	}
	else if (FUNC_COUNT == 6)
	{
		VECTOR v(6);
		v[0] = 0.0;
		v[1] = 1.0/grid.cells[iCell].HX;
		v[2] = 0.0;
		v[3] = 2.0*(pt.x - grid.cells[iCell].c.x)/(grid.cells[iCell].HX*grid.cells[iCell].HX);
		v[4] = 0.0;
		v[5] = (pt.y - grid.cells[iCell].c.y)/(grid.cells[iCell].HX*grid.cells[iCell].HY);
		return v;
	} else {
		log("ERROR: wrong basic functions count!\n");
		EXIT(1);
	}
}

VECTOR FEM_RKDG::getDFDY(int iCell, Point pt)
{
	if (FUNC_COUNT == 3 )
	{
		VECTOR v(3);
		v[0] = 0.0;
		v[1] = 0.0;
		v[2] = 1.0/grid.cells[iCell].HY;
		return v;
	}
	else if (FUNC_COUNT == 6)
	{
		VECTOR v(6);
		v[0] = 0.0;
		v[1] = 0.0;
		v[2] = 1.0/grid.cells[iCell].HY;
		v[3] = 0.0;
		v[4] = 2.0*(pt.y - grid.cells[iCell].c.y)/(grid.cells[iCell].HY*grid.cells[iCell].HY);
		v[5] = (pt.x - grid.cells[iCell].c.x)/(grid.cells[iCell].HX*grid.cells[iCell].HY);
		return v;
	} else {
		log("ERROR: wrong basic functions count!\n");
		EXIT(1);
	}
}
