#include "fem_dg_implicit.h"
#include "tinyxml.h"
#include "global.h"
#include <ctime>
#include <cfloat>
#include "Limiter.h"
#include "MatrixSolver.h"
#include "MeshReader.h"
#include "ANN/ANN.h"

#define POW_2(x) ((x)*(x))
#define F_HLLC_U(UK, FK, SK, SS, PK, RK, VK) (((SS)*((SK)*(UK)-(FK)) + (SK)*( (PK)+(RK)*((SK)-(VK))*((SS)-(VK)) )) / ((SK)-(SS)))
#define F_HLLC_V(UK, FK, SK, SS, PK, RK, VK) (((SS)*((SK)*(UK)-(FK))) / ((SK)-(SS)))
#define F_HLLC_E(UK, FK, SK, SS, PK, RK, VK) (((SS)*((SK)*(UK)-(FK)) + (SK)*( (PK)+(RK)*((SK)-(VK))*((SS)-(VK)) )*(SS)) / ((SK)-(SS)))

void FEM_DG_IMPLICIT::init(char * xmlFileName)
{
	TiXmlDocument doc(xmlFileName);
	bool loadOkay = doc.LoadFile(TIXML_ENCODING_UTF8);
	if (!loadOkay)
	{
		log("ERROR: %s\n", doc.ErrorDesc());
		exit(doc.ErrorId());
	}

	TiXmlNode* task = 0;
	TiXmlElement* el = 0;
	TiXmlNode* node0 = 0;
	TiXmlNode* node1 = 0;
	task = doc.FirstChild("task");

	int steadyVal = 1;
	int viscosityVal = 1;
	int turbulenceVal = 1;
	node0 = task->FirstChild("control");
	node0->FirstChild("STEADY")->ToElement()->Attribute("value", &steadyVal);
    node0->FirstChild("VISCOSITY")->ToElement()->Attribute("value", &viscosityVal);
	node0->FirstChild("SIGMA")->ToElement()->Attribute("value", &SIGMA);
	node0->FirstChild("TAU")->ToElement()->Attribute("value", &TAU);
	node0->FirstChild("TMAX")->ToElement()->Attribute("value", &TMAX);
	node0->FirstChild("STEP_MAX")->ToElement()->Attribute("value", &STEP_MAX);
	node0->FirstChild("FILE_OUTPUT_STEP")->ToElement()->Attribute("value", &FILE_SAVE_STEP);
	node0->FirstChild("LOG_OUTPUT_STEP")->ToElement()->Attribute("value", &PRINT_STEP);
	node0->FirstChild("TURBULENCE")->ToElement()->Attribute("value", &turbulenceVal);

	const char * limiterName = node0->FirstChild("LIMITER")->ToElement()->Attribute("value");

	const char * flxStr = node0->FirstChild("FLUX")->ToElement()->Attribute("value");
	if (strcmp(flxStr, "GODUNOV") == 0) {
		FLUX = FLUX_GODUNOV;
		log("Flux: Godunov\n");
	}
	else if (strcmp(flxStr, "LAX") == 0) {
		FLUX = FLUX_LAX;
		log("Flux: Lax-Friedrichs\n");
	}
	else if (strcmp(flxStr, "CD") == 0) {
		FLUX = FLUX_CD;
		log("Flux: CD\n");
	}
	else if (strcmp(flxStr, "KIR") == 0) {
		FLUX = FLUX_KIR;
		log("Flux: KIR\n");
	}
	else if (strcmp(flxStr, "AUSM+") == 0) {
		FLUX = FLUX_AUSMP;
		log("Flux: AUSM+\n");
	}
	else if (strcmp(flxStr, "AUSMPW+") == 0) {
		FLUX = FLUX_AUSMPW;
		log("Flux: AUSMPW+\n");
	}
	else if (strcmp(flxStr, "HLLC") == 0) {
		FLUX = FLUX_HLLC;
		log("Flux: HLLC\n");
	}
	else {
		FLUX = FLUX_GODUNOV;
		log("Flux: Godunov\n");
	}

	if (steadyVal == 0) {
		STEADY = false;
	}
	else {
		STEADY = true;
		node1 = node0->FirstChild("CFL");
		node1->FirstChild("start")->ToElement()->Attribute("value", &CFL);
		node1->FirstChild("scale")->ToElement()->Attribute("value", &scaleCFL);
		node1->FirstChild("max")->ToElement()->Attribute("value", &maxCFL);
		node1->FirstChild("step")->ToElement()->Attribute("value", &stepCFL);
		node1->FirstChild("max_limited_cells")->ToElement()->Attribute("value", &maxLimCells);
	}

	if (viscosityVal == 1){
	    VISCOSITY = true;
		fieldCount = FIELD_COUNT_EXT;
		matrDim = FIELD_COUNT_EXT * BASE_FUNC_COUNT;
		log("Viscosity: true\n");
	}
	else{
	    VISCOSITY = false;
		fieldCount = FIELD_COUNT;
		matrDim = FIELD_COUNT * BASE_FUNC_COUNT;
		log("Viscosity: false\n");
	}

	if (turbulenceVal == 0) {
		TURBULENCE = false;
	}
	else {
		TURBULENCE = true;
		const char* turbModelStr = node0->FirstChild("TURBULENCE_MODEL")->ToElement()->Attribute("value");
		if (strcmp(turbModelStr, "SA") == 0) {
			TURBULENCE_MODEL = MODEL_SA;
			log("Turbulence model: Spalart–Allmaras\n");
		}
	}

	// сглаживание невязок
	int smUsing = 1;
	node0 = task->FirstChild("smoothing");
	node0->FirstChild("using")->ToElement()->Attribute("value", &smUsing);
	node0->FirstChild("coefficient")->ToElement()->Attribute("value", &SMOOTHING_PAR);
	SMOOTHING = (smUsing == 1);

	// чтение параметров о ПРЕДЕЛЬНЫХ ЗНАЧЕНИЯХ
	node0 = task->FirstChild("limits");
	node0->FirstChild("ro")->ToElement()->Attribute("min", &limitRmin);
	node0->FirstChild("ro")->ToElement()->Attribute("max", &limitRmax);
	node0->FirstChild("p")->ToElement()->Attribute("min", &limitPmin);
	node0->FirstChild("p")->ToElement()->Attribute("max", &limitPmax);
	node0->FirstChild("u")->ToElement()->Attribute("max", &limitUmax);

	// чтение параметров о СТАНДАРТНЫХ ЗНАЧЕНИЯХ
	node0 = task->FirstChild("references");
	node0->FirstChild("ro")->ToElement()->Attribute("value", &RRef);
	node0->FirstChild("p")->ToElement()->Attribute("value", &PRef);

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
		node1->FirstChild("M")->ToElement()->Attribute("value", &mat.M);
		node1->FirstChild("Cp")->ToElement()->Attribute("value", &mat.Cp);
		node1->FirstChild("K")->ToElement()->Attribute("value", &mat.K);
		node1->FirstChild("ML")->ToElement()->Attribute("value", &mat.ML);
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

	for (int i = 0; i < regCount; i++) {
		Region& reg = regions[i];
		regNode->ToElement()->Attribute("id", &reg.id);
		regNode->FirstChild("material")->ToElement()->Attribute("id", &reg.matId);
		regNode->FirstChild("cell")->ToElement()->Attribute("type", &reg.cellType);

		reg.name = regNode->FirstChild("name")->ToElement()->GetText();

		node1 = regNode->FirstChild("parameters");
		node1->FirstChild("Vx")->ToElement()->Attribute("value", &reg.par.u);
		node1->FirstChild("Vy")->ToElement()->Attribute("value", &reg.par.v);
		node1->FirstChild("Vz")->ToElement()->Attribute("value", &reg.par.w);
		node1->FirstChild("T")->ToElement()->Attribute("value", &reg.par.T);
		node1->FirstChild("P")->ToElement()->Attribute("value", &reg.par.p);

		Material& mat = materials[reg.matId];
		mat.URS(reg.par, 2);    // r=r(p,T)
		mat.URS(reg.par, 1);    // e=e(p,r)

		regNode = regNode->NextSibling("region");

		if (reg.par.p > maxP) maxP = reg.par.p;
		if (reg.par.r > maxR) maxR = reg.par.r;
		if (reg.par.T > maxT) maxT = reg.par.T;
		if (reg.par.U2() > maxU) maxU = reg.par.U2();
	}

	maxU = sqrt(maxU);

	// параметры обезразмеривания
	L_ = 1.0;
	R_ = maxR;					// характерная плотность = начальная плотность
	P_ = maxP;					// характерное давление = начальное давление
	T_ = maxT;					// характерная температура = начальная температура
	U_ = sqrt(P_ / R_);		// характерная скорость = sqrt( P_ / R_ )
	E_ = POW_2(U_);			// характерная энергия  = U_**2
	CV_ = POW_2(U_) / T_;	// характерная теплоёмкость  = U_**2 / T_
	TIME_ = L_ / U_;			// характерное время
	MU_ = R_ * U_ * L_;		// характерная вязкость = R_ * U_ * L_
	KP_ = R_ * POW_2(U_) * U_ * L_ / T_;	// коэффициент теплопроводности = R_ * U_**3 * L_ / T_
	CV_ = POW_2(U_) / T_;	// характерная теплоёмкость  = U_**2 / T_
    TAU_ = R_ * U_ / L_;    // характерный компонент тензора вязких напряжений

	// заполняем параметры для лимитера
	limParam.Pmin = limitPmin / P_;
	limParam.PRef = PRef / P_;
	limParam.RRef = RRef / R_;

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
		par.w /= U_;
		par.T /= T_;

		par.r /= R_;
		par.e /= E_;
		par.cz /= U_;
		par.ML /= MU_;
		par.E = par.e + par.U2()*0.5;
	}

	Material::gR *= R_ * T_ / P_;	// Газовая постоянная
	for (int i = 0; i < matCount; i++)
	{
		Material & mat = materials[i];
		mat.Cp /= CV_;
		mat.K /= KP_;
		mat.ML /= MU_;
	}



	// чтение параметров о ГРАНИЧНЫХ УСЛОВИЯХ
    node0 = task->FirstChild("boundaries");
    TiXmlNode* bNode = node0->FirstChild("boundCond");
    while (bNode != NULL)
    {
        int faceType;
        bNode->ToElement()->Attribute("faceType", &faceType);

        CFDBoundary * b;

        try {
            b = CFDBoundary::create(bNode, &grid);
        }
        catch (Exception e) {
            log("ERROR: %s\n", e.getMessage());
            exit(e.getType());
        }

		if (b->faceType == CFDBoundary::TYPE_ID_INLET_SUB) {
			b->par[0] /= R_;
			b->par[1] /= P_;
			b->par[2] /= T_;
		}
		if (b->faceType == CFDBoundary::TYPE_ID_OUTLET_SUB) {
			b->par[0] /= P_;
		}
		if (b->faceType == CFDBoundary::TYPE_ID_PRESSURE) {
			b->par[0] /= P_;
			b->par[1] /= T_;
		}
		if (b->faceType == CFDBoundary::TYPE_ID_WALL_NO_SLIP) {
			b->par[0] /= T_;
		}
		if (b->faceType == CFDBoundary::TYPE_ID_MASSFLOW) {
			b->par[0] /= P_;
			b->par[2] /= T_;
		}
		if (b->faceType == CFDBoundary::TYPE_ID_FREE_STREAM) {
			b->par[0] /= P_;
			b->par[1] /= T_;
		}
		if (b->faceType == CFDBoundary::TYPE_ID_STAGNATION) {
			b->par[0] /= P_;
			b->par[1] /= P_;
			b->par[2] /= T_;
		}

        boundaries.push_back(b);

        bNode = bNode->NextSibling("boundCond");
    }

    bCount = boundaries.size();
	
	node0 = task->FirstChild("mesh");
	const char* fName = node0->FirstChild("name")->ToElement()->Attribute("value");
	const char* tName = node0->FirstChild("filesType")->ToElement()->Attribute("value");
	MeshReader* mr = MeshReader::create(MeshReader::getType((char*)tName), (char*)fName);
	mr->read(&grid);

	// задаем фиксированный вектор и определяем ориентацию на гранях ячеек
	//defineOrientationOnFaces();

	memAlloc();

	/* определение ГУ для каждой грани. */
	for (int iFace = 0; iFace < grid.fCount; iFace++) {
		Face & f = grid.faces[iFace];
		if (f.type == Face::TYPE_INNER) {
			f.bnd = NULL;
			continue;
		}
		if (f.type == Face::TYPE_NAMED) {
			int iBound = -1;
			for (int i = 0; i < bCount; i++)
			{
				if (strcmp(f.typeName, boundaries[i]->name) == 0)
				{
					iBound = i;
					break;
				}
			}
			if (iBound < 0)
			{
				log("ERROR (boundary condition): unknown face type of face %d...\n", iFace);
				EXIT(1);
			}

			f.bnd = boundaries[iBound];
		}
		else {
			int iBound = -1;
			for (int i = 0; i < bCount; i++)
			{
				if (strcmp(f.typeName, boundaries[i]->name) == 0)
				{
					iBound = i;
					break;
				}
			}
			if (iBound < 0)
			{
				log("ERROR (boundary condition): unknown face type of face %d...\n", iFace);
				EXIT(1);
			}
			f.bnd = boundaries[iBound];
			strcpy(f.typeName, boundaries[iBound]->name);
		}
	}

	// инициализация лимитера
	limiter = Limiter::create(limiterName, this);
	if (limiter != NULL) {
		log("Limiter: %s\n", limiter->getName());
	}
	else {
		log("Without limiter\n");
	}

	// инициализация решателя
	node0 = task->FirstChild("solver");
	const char * solverName = node0->ToElement()->Attribute("type");
	solverMtx = MatrixSolver::create(solverName);
	//solverMtx->init(grid.cCount, MATR_DIM);
	solverMtx->init(&(grid), matrDim, 1);
	log("Solver type: %s.\n", solverMtx->getName());

	node0->FirstChild("iterations")->ToElement()->Attribute("value", &SOLVER_ITER);
	node0->FirstChild("epsilon")->ToElement()->Attribute("value", &SOLVER_EPS);
	if (strstr(solverName, "HYPRE") != NULL) {
		node1 = node0->FirstChild("hypre");

		int tmp = 0;
		node1->FirstChild("printLevel")->ToElement()->Attribute("value", &tmp);
		solverMtx->setParameter("PRINT_LEVEL", tmp);
		node1->FirstChild("krylovDim")->ToElement()->Attribute("value", &tmp);
		solverMtx->setParameter("KRYLOV_DIM", tmp);
	}

	calcGaussPar();

	calcMassMatr();

	for (int i = 0; i < grid.cCount; i++)
	{
		Cell& c = grid.cells[i];
		Region& reg = getRegion(c.typeName);
		convertParToCons(i, reg.par);
		if (VISCOSITY) {
			memset(tau_xx[i], 0, sizeof(double) * BASE_FUNC_COUNT);
			memset(tau_yy[i], 0, sizeof(double) * BASE_FUNC_COUNT);
			memset(tau_zz[i], 0, sizeof(double) * BASE_FUNC_COUNT);
			memset(tau_xy[i], 0, sizeof(double) * BASE_FUNC_COUNT);
			memset(tau_xz[i], 0, sizeof(double) * BASE_FUNC_COUNT);
			memset(tau_yz[i], 0, sizeof(double) * BASE_FUNC_COUNT);
		}
	}

	calcTimeStep();

	save(0);
}

void FEM_DG_IMPLICIT::calcMassMatr()
{
	for (int iCell = 0; iCell < grid.cCount; iCell++) {
		double **	A = matrA[iCell];
		for (int i = 0; i < BASE_FUNC_COUNT; i++) {
			for (int j = 0; j < BASE_FUNC_COUNT; j++) {
				A[i][j] = 0.0;
				for (int iGP = 0; iGP < GP_CELL_COUNT; iGP++) {
					A[i][j] += cellGW[iCell][iGP] * getF(i, iCell, cellGP[iCell][iGP])
						                          * getF(j, iCell, cellGP[iCell][iGP]);
				}
				A[i][j] *= cellJ[iCell];
			}
		}
	}
}

void FEM_DG_IMPLICIT::calcGaussPar()
{
	
	// для ячеек
	if (GP_CELL_COUNT == 5) {
		for (int i = 0; i < grid.cCount; i++) {

			double a = 1.0 / 4.0;
			double b = 1.0 / 2.0;
			double c = 1.0 / 6.0;

			double x1 = grid.nodes[grid.cells[i].nodesInd[0]].x;
			double y1 = grid.nodes[grid.cells[i].nodesInd[0]].y;
			double z1 = grid.nodes[grid.cells[i].nodesInd[0]].z;
			double x2 = grid.nodes[grid.cells[i].nodesInd[1]].x;
			double y2 = grid.nodes[grid.cells[i].nodesInd[1]].y;
			double z2 = grid.nodes[grid.cells[i].nodesInd[1]].z;
			double x3 = grid.nodes[grid.cells[i].nodesInd[2]].x;
			double y3 = grid.nodes[grid.cells[i].nodesInd[2]].y;
			double z3 = grid.nodes[grid.cells[i].nodesInd[2]].z;
			double x4 = grid.nodes[grid.cells[i].nodesInd[3]].x;
			double y4 = grid.nodes[grid.cells[i].nodesInd[3]].y;
			double z4 = grid.nodes[grid.cells[i].nodesInd[3]].z;
			double a1 = x1 - x4; double a2 = x2 - x4; double a3 = x3 - x4; double a4 = x4;
			double b1 = y1 - y4; double b2 = y2 - y4; double b3 = y3 - y4; double b4 = y4;
			double c1 = z1 - z4; double c2 = z2 - z4; double c3 = z3 - z4; double c4 = z4;

			cellGW[i][0] = -4.0 / 30.0;
			cellGP[i][0].x = a * a1 + a * a2 + a * a3 + a4;
			cellGP[i][0].y = a * b1 + a * b2 + a * b3 + b4;
			cellGP[i][0].z = a * c1 + a * c2 + a * c3 + c4;

			cellGW[i][1] = 9.0 / 120.0;
			cellGP[i][1].x = b * a1 + c * a2 + c * a3 + a4;
			cellGP[i][1].y = b * b1 + c * b2 + c * b3 + b4;
			cellGP[i][1].z = b * c1 + c * c2 + c * c3 + c4;

			cellGW[i][2] = 9.0 / 120.0;
			cellGP[i][2].x = c * a1 + b * a2 + c * a3 + a4;
			cellGP[i][2].y = c * b1 + b * b2 + c * b3 + b4;
			cellGP[i][2].z = c * c1 + b * c2 + c * c3 + c4;

			cellGW[i][3] = 9.0 / 120.0;
			cellGP[i][3].x = c * a1 + c * a2 + b * a3 + a4;
			cellGP[i][3].y = c * b1 + c * b2 + b * b3 + b4;
			cellGP[i][3].z = c * c1 + c * c2 + b * c3 + c4;

			cellGW[i][4] = 9.0 / 120.0;
			cellGP[i][4].x = c * a1 + c * a2 + c * a3 + a4;
			cellGP[i][4].y = c * b1 + c * b2 + c * b3 + b4;
			cellGP[i][4].z = c * c1 + c * c2 + c * c3 + c4;

			cellJ[i] = fabs(a1 * (b2 * c3 - b3 * c2) - a2 * (b1 * c3 - b3 * c1) + a3 * (b1 * c2 - b2 * c1));
		}
	}
	else if (GP_CELL_COUNT == 4) {
		for (int i = 0; i < grid.cCount; i++) {

			double a = 0.5854101966249685;
			double b = 0.1381966011250105;

			double x1 = grid.nodes[grid.cells[i].nodesInd[0]].x;
			double y1 = grid.nodes[grid.cells[i].nodesInd[0]].y;
			double z1 = grid.nodes[grid.cells[i].nodesInd[0]].z;

			double x2 = grid.nodes[grid.cells[i].nodesInd[1]].x;
			double y2 = grid.nodes[grid.cells[i].nodesInd[1]].y;
			double z2 = grid.nodes[grid.cells[i].nodesInd[1]].z;

			double x3 = grid.nodes[grid.cells[i].nodesInd[2]].x;
			double y3 = grid.nodes[grid.cells[i].nodesInd[2]].y;
			double z3 = grid.nodes[grid.cells[i].nodesInd[2]].z;

			double x4 = grid.nodes[grid.cells[i].nodesInd[3]].x;
			double y4 = grid.nodes[grid.cells[i].nodesInd[3]].y;
			double z4 = grid.nodes[grid.cells[i].nodesInd[3]].z;

			double a1 = x1 - x4; double a2 = x2 - x4; double a3 = x3 - x4; double a4 = x4;
			double b1 = y1 - y4; double b2 = y2 - y4; double b3 = y3 - y4; double b4 = y4;
			double c1 = z1 - z4; double c2 = z2 - z4; double c3 = z3 - z4; double c4 = z4;

			cellGW[i][0] = 1.0 / 24.0;
			cellGP[i][0].x = a * a1 + b * a2 + b * a3 + a4;
			cellGP[i][0].y = a * b1 + b * b2 + b * b3 + b4;
			cellGP[i][0].z = a * c1 + b * c2 + b * c3 + c4;

			cellGW[i][1] = 1.0 / 24.0;
			cellGP[i][1].x = b * a1 + a * a2 + b * a3 + a4;
			cellGP[i][1].y = b * b1 + a * b2 + b * b3 + b4;
			cellGP[i][1].z = b * c1 + a * c2 + b * c3 + c4;

			cellGW[i][2] = 1.0 / 24.0;
			cellGP[i][2].x = b * a1 + b * a2 + a * a3 + a4;
			cellGP[i][2].y = b * b1 + b * b2 + a * b3 + b4;
			cellGP[i][2].z = b * c1 + b * c2 + a * c3 + c4;

			cellGW[i][3] = 1.0 / 24.0;
			cellGP[i][3].x = b * a1 + b * a2 + b * a3 + a4;
			cellGP[i][3].y = b * b1 + b * b2 + b * b3 + b4;
			cellGP[i][3].z = b * c1 + b * c2 + b * c3 + c4;

			cellJ[i] = fabs(a1 * (b2 * c3 - b3 * c2) - a2 * (b1 * c3 - b3 * c1) + a3 * (b1 * c2 - b2 * c1));
		}
	}

	// для граней
	if (GP_FACE_COUNT == 4) {
		for (int i = 0; i < grid.fCount; i++) {

			double a = 1.0 / 3.0;
			double b = 1.0 / 5.0;
			double c = 3.0 / 5.0;
			double x1 = grid.nodes[grid.faces[i].nodesInd[0]].x;
			double y1 = grid.nodes[grid.faces[i].nodesInd[0]].y;
			double z1 = grid.nodes[grid.faces[i].nodesInd[0]].z;

			double x2 = grid.nodes[grid.faces[i].nodesInd[1]].x;
			double y2 = grid.nodes[grid.faces[i].nodesInd[1]].y;
			double z2 = grid.nodes[grid.faces[i].nodesInd[1]].z;

			double x3 = grid.nodes[grid.faces[i].nodesInd[2]].x;
			double y3 = grid.nodes[grid.faces[i].nodesInd[2]].y;
			double z3 = grid.nodes[grid.faces[i].nodesInd[2]].z;

			double a1 = x1 - x3; double a2 = x2 - x3; double a3 = x3;
			double b1 = y1 - y3; double b2 = y2 - y3; double b3 = y3;
			double c1 = z1 - z3; double c2 = z2 - z3; double c3 = z3;

			faceGW[i][0] = -27.0 / 48.0;
			faceGP[i][0].x = a1 * a + a2 * a + a3;
			faceGP[i][0].y = b1 * a + b2 * a + b3;
			faceGP[i][0].z = c1 * a + c2 * a + c3;

			faceGW[i][1] = 25.0 / 48.0;
			faceGP[i][1].x = a1 * c + a2 * b + a3;
			faceGP[i][1].y = b1 * c + b2 * b + b3;
			faceGP[i][1].z = c1 * c + c2 * b + c3;

			faceGW[i][2] = 25.0 / 48.0;
			faceGP[i][2].x = a1 * b + a2 * c + a3;
			faceGP[i][2].y = b1 * b + b2 * c + b3;
			faceGP[i][2].z = c1 * b + c2 * c + c3;

			faceGW[i][3] = 25.0 / 48.0;
			faceGP[i][3].x = a1 * b + a2 * b + a3;
			faceGP[i][3].y = b1 * b + b2 * b + b3;
			faceGP[i][3].z = c1 * b + c2 * b + c3;

			faceJ[i] = sqrt((b1 * c2 - c1 * b2) * (b1 * c2 - c1 * b2) + (a1 * c2 - c1 * a2) * (a1 * c2 - c1 * a2) + (a1 * b2 - b1 * a2) * (a1 * b2 - b1 * a2));
		}
	}
	else if (GP_FACE_COUNT == 3) {
		for (int i = 0; i < grid.fCount; i++) {
			double a = 1.0 / 6.0;
			double b = 2.0 / 3.0;
			double x1 = grid.nodes[grid.faces[i].nodesInd[0]].x;
			double y1 = grid.nodes[grid.faces[i].nodesInd[0]].y;
			double z1 = grid.nodes[grid.faces[i].nodesInd[0]].z;

			double x2 = grid.nodes[grid.faces[i].nodesInd[1]].x;
			double y2 = grid.nodes[grid.faces[i].nodesInd[1]].y;
			double z2 = grid.nodes[grid.faces[i].nodesInd[1]].z;

			double x3 = grid.nodes[grid.faces[i].nodesInd[2]].x;
			double y3 = grid.nodes[grid.faces[i].nodesInd[2]].y;
			double z3 = grid.nodes[grid.faces[i].nodesInd[2]].z;

			double a1 = x1 - x3; double a2 = x2 - x3; double a3 = x3;
			double b1 = y1 - y3; double b2 = y2 - y3; double b3 = y3;
			double c1 = z1 - z3; double c2 = z2 - z3; double c3 = z3;

			faceGW[i][0] = 1.0 / 6.0;
			faceGP[i][0].x = a1 * b + a2 * a + a3;
			faceGP[i][0].y = b1 * b + b2 * a + b3;
			faceGP[i][0].z = c1 * b + c2 * a + c3;

			faceGW[i][1] = 1.0 / 6.0;
			faceGP[i][1].x = a1 * a + a2 * b + a3;
			faceGP[i][1].y = b1 * a + b2 * b + b3;
			faceGP[i][1].z = c1 * a + c2 * b + c3;

			faceGW[i][2] = 1.0 / 6.0;
			faceGP[i][2].x = a1 * a + a2 * a + a3;
			faceGP[i][2].y = b1 * a + b2 * a + b3;
			faceGP[i][2].z = c1 * a + c2 * a + c3;

			faceJ[i] = sqrt((b1 * c2 - c1 * b2) * (b1 * c2 - c1 * b2) + (a1 * c2 - c1 * a2) * (a1 * c2 - c1 * a2) + (a1 * b2 - b1 * a2) * (a1 * b2 - b1 * a2));

		}
	}
}

double** allocMtx(int N)
{
	double **m;
	m = new double*[N];
	for (int i = 0; i < N; i++)
		m[i] = new double[N];

	return m;
}

void freeMtx(double** m, int N)
{
	for (int i = 0; i < N; i++)
		delete[] m[i];
	delete[] m;
}

void FEM_DG_IMPLICIT::memAlloc()
{
	int n = grid.cCount;
	cTau = new double[n];

	ro = new double* [n];
	ru = new double* [n];
	rv = new double* [n];
	rw = new double* [n];
	re = new double* [n];
	if (VISCOSITY) {
		tau_xx = new double* [n];
		tau_yy = new double* [n];
		tau_zz = new double* [n];
		tau_xy = new double* [n];
		tau_xz = new double* [n];
		tau_yz = new double* [n];
	}

	cellGP = new Point*[n];
	cellGW = new double*[n];
	cellJ = new double[n];

	faceGW = new double*[grid.fCount];
	faceJ  = new double[grid.fCount];
	faceGP = new Point*[grid.fCount];

	matrA		= new double**[n];
	matrInvA	= new double**[n];

	for (int i = 0; i < n; i++) {
		ro[i] = new double[BASE_FUNC_COUNT];
		ru[i] = new double[BASE_FUNC_COUNT];
		rv[i] = new double[BASE_FUNC_COUNT];
		rw[i] = new double[BASE_FUNC_COUNT];
		re[i] = new double[BASE_FUNC_COUNT];

		if (VISCOSITY) {
			tau_xx[i] = new double[BASE_FUNC_COUNT];
			tau_yy[i] = new double[BASE_FUNC_COUNT];
			tau_zz[i] = new double[BASE_FUNC_COUNT];
			tau_xy[i] = new double[BASE_FUNC_COUNT];
			tau_xz[i] = new double[BASE_FUNC_COUNT];
			tau_yz[i] = new double[BASE_FUNC_COUNT];
		}

		cellGP[i] = new Point[GP_CELL_COUNT];
		cellGW[i] = new double[GP_CELL_COUNT];

		matrA[i] = allocMtx(BASE_FUNC_COUNT);
		matrInvA[i] = allocMtx(BASE_FUNC_COUNT);
	}

	for (int i = 0; i < grid.fCount; i++) {
		faceGP[i] = new Point[GP_FACE_COUNT];
		faceGW[i] = new double[GP_FACE_COUNT];
	}

	tmpArr = new double[n];
	tmpArr1 = new double[n];
	tmpArr2 = new double[n];
	tmpCFL = new double[n];
	tmpArrInt = new int[n];

	fields = new double**[fieldCount];
	fields[FIELD_RO] = ro;
	fields[FIELD_RU] = ru;
	fields[FIELD_RV] = rv;
	fields[FIELD_RW] = rw;
	fields[FIELD_RE] = re;
	if (VISCOSITY) {
		fields[FIELD_TAU_XX] = tau_xx;
		fields[FIELD_TAU_YY] = tau_yy;
		fields[FIELD_TAU_ZZ] = tau_zz;
		fields[FIELD_TAU_XY] = tau_xy;
		fields[FIELD_TAU_XZ] = tau_xz;
		fields[FIELD_TAU_YZ] = tau_yz;
	}

	matrSmall = allocMtx(BASE_FUNC_COUNT);
	matrSmall2 = allocMtx(BASE_FUNC_COUNT);
	matrBig = allocMtx(matrDim);
	matrBig2 = allocMtx(matrDim);

}

void FEM_DG_IMPLICIT::memFree() 
{
	int n = grid.cCount; 
	
	for (int i = 0; i < n; i++) {
		delete[] ro[i];
		delete[] ru[i];
		delete[] rv[i];
		delete[] rw[i];
		delete[] re[i];

		if (VISCOSITY) {
			delete[] tau_xx[i];
			delete[] tau_yy[i];
			delete[] tau_zz[i];
			delete[] tau_xy[i];
			delete[] tau_xz[i];
			delete[] tau_yz[i];
		}

		delete[] cellGP[i];
		delete[] cellGW[i];

		freeMtx(matrA[i], BASE_FUNC_COUNT);
		freeMtx(matrInvA[i], BASE_FUNC_COUNT);
	}
	delete[] cTau;

	delete[] ro;
	delete[] ru;
	delete[] rv;
	delete[] rw;
	delete[] re;
	if (VISCOSITY) {
		delete[] tau_xx;
		delete[] tau_yy;
		delete[] tau_zz;
		delete[] tau_xy;
		delete[] tau_xz;
		delete[] tau_yz;
	}

	delete[] cellGP;
	delete[] cellGW;
	delete[] cellJ;

	delete[] matrA;
	delete[] matrInvA;

	for (int i = 0; i < grid.fCount; i++) {
		delete[] faceGP[i];
		delete[] faceGW[i];
	}
	delete[] faceGW;
	delete[] faceJ;
	delete[] faceGP;

	delete[] tmpArr;
	delete[] tmpArr1;
	delete[] tmpArr2;
	delete[] tmpCFL;

	for (int i = 0; i < fieldCount; i++) {
		delete[] fields[i];
	}
	delete[] fields;

	freeMtx(matrSmall, BASE_FUNC_COUNT);
	freeMtx(matrSmall2, BASE_FUNC_COUNT);
	freeMtx(matrBig, matrDim);
	freeMtx(matrBig2, matrDim);

	delete[] tmpArr;
	delete[] tmpArr1;
	delete[] tmpArr2;
	delete[] tmpArrInt;
}

Region & FEM_DG_IMPLICIT::getRegionByCellType(int type)
{
	for (int i = 0; i < regCount; i++)
	{
		if (regions[i].cellType == type) return regions[i];
	}
	log("ERROR: unknown cell type %d...\n", type);
	EXIT(1);
}


Region   &	FEM_DG_IMPLICIT::getRegion(int iCell)
{
	return getRegionByCellType(grid.cells[iCell].type);
}

Region& FEM_DG_IMPLICIT::getRegionByName(char* name)
{
	for (int i = 0; i < regCount; i++)
	{
		if (strcmp(regions[i].name.c_str(), name) == 0) return regions[i];
	}
	log("ERROR: unknown cell name '%s'...\n", name);
	EXIT(1);
}

Region& FEM_DG_IMPLICIT::getRegion(char* name)
{
	return getRegionByName(name);
}

Material &	FEM_DG_IMPLICIT::getMaterial(int iCell)
{
	Region& reg = getRegion(grid.cells[iCell].typeName);
	return materials[reg.matId];
}


void FEM_DG_IMPLICIT::convertParToCons(int iCell, Param & par)
{
	memset(ro[iCell], 0, sizeof(double) * BASE_FUNC_COUNT);
	memset(ru[iCell], 0, sizeof(double) * BASE_FUNC_COUNT);
	memset(rv[iCell], 0, sizeof(double) * BASE_FUNC_COUNT);
	memset(rw[iCell], 0, sizeof(double) * BASE_FUNC_COUNT);
	memset(re[iCell], 0, sizeof(double) * BASE_FUNC_COUNT);
	ro[iCell][0] = par.r;
	ru[iCell][0] = par.r * par.u;
	rv[iCell][0] = par.r * par.v;
	rw[iCell][0] = par.r * par.w;
	re[iCell][0] = par.r * (par.e + 0.5 * (par.u * par.u + par.v * par.v + par.w * par.w));
}

void FEM_DG_IMPLICIT::convertConsToPar(int iCell, Param & par)
{
	double fRO = getField(FIELD_RO, iCell, grid.cells[iCell].c);
	double fRU = getField(FIELD_RU, iCell, grid.cells[iCell].c);
	double fRV = getField(FIELD_RV, iCell, grid.cells[iCell].c);
	double fRW = getField(FIELD_RW, iCell, grid.cells[iCell].c);
	double fRE = getField(FIELD_RE, iCell, grid.cells[iCell].c);

	par.r = fRO;
	par.u = fRU / fRO;
	par.v = fRV / fRO;
	par.w = fRW / fRO;
	par.E = fRE / fRO;
	par.e = par.E - 0.5 * par.U2();
	Material& mat = getMaterial(iCell);
	mat.URS(par, 0);
}

double FEM_DG_IMPLICIT::getField(int fldId, int iCell, double x, double y, double z)
{
	double *fld = fields[fldId][iCell];
	double result = 0.0;
	for (int i = 0; i < BASE_FUNC_COUNT; i++) {
		result += fld[i] * getF(i, iCell, x, y, z);
	}
	return result;
}

double FEM_DG_IMPLICIT::getField(int fldId, int iCell, Point p)
{
	return getField(fldId, iCell, p.x, p.y, p.z);
}

double FEM_DG_IMPLICIT::getFieldLinear(int fldId, int iCell, double x, double y, double z)
{
	double* fld = fields[fldId][iCell];
	double result = 0.0;
	for (int i = 0; i < 4; i++) { // разложение по линейным функциям
		result += fld[i] * getF(i, iCell, x, y, z);
	}
	return result;
}

double FEM_DG_IMPLICIT::getFieldLinear(int fldId, int iCell, Point p)
{
	return getFieldLinear(fldId, iCell, p.x, p.y, p.z);
}

void FEM_DG_IMPLICIT::getFields(double &fRO, double &fRU, double &fRV, double &fRW, double &fRE, int iCell, double x, double y, double z)
{
	fRO = getField(FIELD_RO, iCell, x, y, z);
	fRU = getField(FIELD_RU, iCell, x, y, z);
	fRV = getField(FIELD_RV, iCell, x, y, z);
	fRW = getField(FIELD_RW, iCell, x, y, z);
	fRE = getField(FIELD_RE, iCell, x, y, z);
}

void FEM_DG_IMPLICIT::getFields(double &fRO, double &fRU, double &fRV, double &fRW, double &fRE, int iCell, Point p)
{
	getFields(fRO, fRU, fRV, fRW, fRE, iCell, p.x, p.y, p.z);
}

void FEM_DG_IMPLICIT::getFieldsLinear(double& fRO, double& fRU, double& fRV, double& fRW, double& fRE, int iCell, double x, double y, double z)
{
	fRO = getFieldLinear(FIELD_RO, iCell, x, y, z);
	fRU = getFieldLinear(FIELD_RU, iCell, x, y, z);
	fRV = getFieldLinear(FIELD_RV, iCell, x, y, z);
	fRW = getFieldLinear(FIELD_RW, iCell, x, y, z);
	fRE = getFieldLinear(FIELD_RE, iCell, x, y, z);
}

void FEM_DG_IMPLICIT::getFieldsLinear(double& fRO, double& fRU, double& fRV, double& fRW, double& fRE, int iCell, Point p)
{
	getFieldsLinear(fRO, fRU, fRV, fRW, fRE, iCell, p.x, p.y, p.z);
}

double FEM_DG_IMPLICIT::getF(int id, int iCell, Point p)
{
	return getF(id, iCell, p.x, p.y, p.z);
}

double FEM_DG_IMPLICIT::getF(int id, int iCell, double x, double y, double z)
{
	Point	&c = grid.cells[iCell].c;
	double	&hx = grid.cells[iCell].HX;
	double	&hy = grid.cells[iCell].HY;
	double	&hz = grid.cells[iCell].HZ;

	switch (id) {
	case 0:
		return 1.0;
		break;
	case 1:
		return (x - c.x) / hx;
		break;
	case 2:
		return (y - c.y) / hy;
		break;
	case 3:
		return (z - c.z) / hz;
		break;
	case 4:
		return (x - c.x) * (x - c.x) / hx / hx;
		break;
	case 5:
		return (y - c.y) * (y - c.y) / hy / hy;
		break;
	case 6:
		return (z - c.z) * (z - c.z) / hz / hz;
		break;
	case 7:
		return (x - c.x) * (y - c.y) / hx / hy;
		break;
	case 8:
		return (x - c.x) * (z - c.z) / hx / hz;
		break;
	case 9:
		return (y - c.y) * (z - c.z) / hy / hz;
		break;
	default:
		return 0.0;
		break;
	}
}

double FEM_DG_IMPLICIT::getDfDx(int id, int iCell, Point p)
{
	return getDfDx(id, iCell, p.x, p.y, p.z);
}

double FEM_DG_IMPLICIT::getDfDx(int id, int iCell, double x, double y, double z)
{
	Point &c = grid.cells[iCell].c;
	double &hx = grid.cells[iCell].HX;
	double &hy = grid.cells[iCell].HY;
	double &hz = grid.cells[iCell].HZ;

	switch (id) {
	case 0:
		return 0.0;
		break;
	case 1:
		return 1.0 / hx;
		break;
	case 2:
		return 0.0;
		break;
	case 3:
		return 0.0;
		break;
	case 4:
		return 2.0*(x - c.x) / hx / hx;
		break;
	case 5:
		return 0.0;
		break;
	case 6:
		return 0.0;
		break;
	case 7:
		return (y - c.y) / hx / hy;
		break;
	case 8:
		return (z - c.z) / hx / hz;
		break;
	case 9:
		return 0.0;
		break;
	default:
		return 0.0;
		break;
	}
}

double FEM_DG_IMPLICIT::getDfDy(int id, int iCell, Point p)
{
	return getDfDy(id, iCell, p.x, p.y, p.z);
}

double FEM_DG_IMPLICIT::getDfDy(int id, int iCell, double x, double y, double z)
{
	Point &c = grid.cells[iCell].c;
	double &hx = grid.cells[iCell].HX;
	double &hy = grid.cells[iCell].HY;
	double &hz = grid.cells[iCell].HZ;

	switch (id) {
	case 0:
		return 0.0;
		break;
	case 1:
		return 0.0;
		break;
	case 2:
		return 1.0 / hy;
		break;
	case 3:
		return 0.0;
		break;
	case 4:
		return 0.0;
		break;
	case 5:
		return 2.0*(y - c.y) / hy / hy;
		break;
	case 6:
		return 0.0;
		break;
	case 7:
		return (x - c.x) / hx / hy;
		break;
	case 8:
		return 0.0;
		break;
	case 9:
		return (z - c.z) / hy / hz;
		break;
	default:
		return 0.0;
		break;
	}
}

double FEM_DG_IMPLICIT::getDfDz(int id, int iCell, Point p)
{
	return getDfDz(id, iCell, p.x, p.y, p.z);
}

double FEM_DG_IMPLICIT::getDfDz(int id, int iCell, double x, double y, double z)
{
	Point& c = grid.cells[iCell].c;
	double& hx = grid.cells[iCell].HX;
	double& hy = grid.cells[iCell].HY;
	double& hz = grid.cells[iCell].HZ;

	switch (id) {
	case 0:
		return 0.0;
		break;
	case 1:
		return 0.0;
		break;
	case 2:
		return 0.0;
		break;
	case 3:
		return 1.0 / hz;
		break;
	case 4:
		return 0.0;
		break;
	case 5:
		return 0.0;
		break;
	case 6:
		return 2.0 * (z - c.z) / hz / hz;
		break;
	case 7:
		return 0.0;
		break;
	case 8:
		return (x - c.x) / hx / hz;
		break;
	case 9:
		return (y - c.y) / hy / hz;
		break;
	default:
		return 0.0;
		break;
	}
}

void FEM_DG_IMPLICIT::calcTimeStep()
{
	double tauMin = 1.0e+20;
	for (int iCell = 0; iCell < grid.cCount; iCell++)
	{
		if (STEADY) {
			Param p;
			convertConsToPar(iCell, p);
			double tmp = 0.5 * grid.cells[iCell].V / (sqrt(p.u * p.u + p.v * p.v + p.w * p.w) + p.cz + FLT_EPSILON);
			cTau[iCell] = _min_(CFL*tmp, TAU);
			if (cTau[iCell] < tauMin) tauMin = cTau[iCell];
		} else {
			cTau[iCell] = TAU;
		}
	}
	if (STEADY)	TAU_MIN = tauMin;
}

void FEM_DG_IMPLICIT::save(int step)
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
		fprintf(fp, "%f %f %f  ", grid.nodes[i].x * L_, grid.nodes[i].y * L_, grid.nodes[i].z * L_);
		if (i + 1 % 8 == 0) fprintf(fp, "\n");
	}
	fprintf(fp, "\n");
	fprintf(fp, "CELLS %d %d\n", grid.cCount, 5 * grid.cCount);
	for (int i = 0; i < grid.cCount; i++)
	{
		fprintf(fp, "4 %d %d %d %d\n", grid.cells[i].nodesInd[0], grid.cells[i].nodesInd[1], grid.cells[i].nodesInd[2], grid.cells[i].nodesInd[3]);
	}
	fprintf(fp, "\n");

	fprintf(fp, "CELL_TYPES %d\n", grid.cCount);
	for (int i = 0; i < grid.cCount; i++) fprintf(fp, "10\n");
	fprintf(fp, "\n");

	fprintf(fp, "CELL_DATA %d\nSCALARS Density float 1\nLOOKUP_TABLE default\n", grid.cCount);
	for (int i = 0; i < grid.cCount; i++)
	{
		Param p;
		convertConsToPar(i, p);
		fprintf(fp, "%25.16f ", p.r * R_);
		if (i + 1 % 8 == 0 || i + 1 == grid.cCount) fprintf(fp, "\n");
	}

	fprintf(fp, "SCALARS Pressure float 1\nLOOKUP_TABLE default\n", grid.cCount);
	for (int i = 0; i < grid.cCount; i++)
	{
		Param p;
		convertConsToPar(i, p);
		fprintf(fp, "%25.16f ", p.p * P_);
		if (i + 1 % 8 == 0 || i + 1 == grid.cCount) fprintf(fp, "\n");
	}

	fprintf(fp, "SCALARS Temperature float 1\nLOOKUP_TABLE default\n", grid.cCount);
	for (int i = 0; i < grid.cCount; i++)
	{
		Param p;
		convertConsToPar(i, p);
		Material &mat = getMaterial(i);
		mat.URS(p, 1);
		fprintf(fp, "%25.16f ", p.T * T_);
		if (i + 1 % 8 == 0 || i + 1 == grid.cCount) fprintf(fp, "\n");
	}

	fprintf(fp, "SCALARS MachNumber float 1\nLOOKUP_TABLE default\n", grid.cCount);
	for (int i = 0; i < grid.cCount; i++)
	{
		Param p;
		convertConsToPar(i, p);
		fprintf(fp, "%25.16f ", sqrt(p.u * p.u + p.v * p.v + p.w * p.w) / p.cz);
		if (i + 1 % 8 == 0 || i + 1 == grid.cCount) fprintf(fp, "\n");
	}

	fprintf(fp, "VECTORS Velosity float\n");
	for (int i = 0; i < grid.cCount; i++)
	{
		Param p;
		convertConsToPar(i, p);
		fprintf(fp, "%25.16f %25.16f %25.16f ", p.u * U_, p.v * U_, p.w * U_);
		if (i + 1 % 8 == 0 || i + 1 == grid.cCount) fprintf(fp, "\n");
	}

	if (VISCOSITY) {
		fprintf(fp, "VECTORS TauX float\n");
		for (int i = 0; i < grid.cCount; i++)
		{
			fprintf(fp, "%25.16f %25.16f %25.16f ", tau_xx[i][0] * TAU_, tau_xy[i][0] * TAU_, tau_xz[i][0] * TAU_);
			if (i + 1 % 8 == 0 || i + 1 == grid.cCount) fprintf(fp, "\n");
		}

		fprintf(fp, "VECTORS TauY float\n");
		for (int i = 0; i < grid.cCount; i++)
		{
			fprintf(fp, "%25.16f %25.16f %25.16f ", tau_xy[i][0] * TAU_, tau_yy[i][0] * TAU_, tau_yz[i][0] * TAU_);
			if (i + 1 % 8 == 0 || i + 1 == grid.cCount) fprintf(fp, "\n");
		}

		fprintf(fp, "VECTORS TauZ float\n");
		for (int i = 0; i < grid.cCount; i++)
		{
			fprintf(fp, "%25.16f %25.16f %25.16f ", tau_xz[i][0] * TAU_, tau_yz[i][0] * TAU_, tau_zz[i][0] * TAU_);
			if (i + 1 % 8 == 0 || i + 1 == grid.cCount) fprintf(fp, "\n");
		}
	}

	fprintf(fp, "SCALARS Total_energy float 1\nLOOKUP_TABLE default\n");
	for (int i = 0; i < grid.cCount; i++)
	{
		Param p;
		convertConsToPar(i, p);
		fprintf(fp, "%25.16f ", p.E * E_);
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
		double M2 = (p.u * p.u + p.v * p.v + p.w * p.w) / (p.cz * p.cz);
		fprintf(fp, "%25.16e ", p.p*::pow(1.0 + 0.5*M2*agam, gam / agam)*P_);
		if ((i + 1) % 8 == 0 || i + 1 == grid.cCount) fprintf(fp, "\n");
	}

	fprintf(fp, "SCALARS TAU float 1\nLOOKUP_TABLE default\n", grid.cCount);
	for (int i = 0; i < grid.cCount; i++)
	{
		fprintf(fp, "%25.16f ", cTau[i]*TIME_);
		if (i + 1 % 8 == 0 || i + 1 == grid.cCount) fprintf(fp, "\n");
	}

	fclose(fp);
	printf("File '%s' saved...\n", fName);

}

int FEM_DG_IMPLICIT::getLimitedCellsCount() {
	int n = 0;
	for (int iCell = 0; iCell < grid.cCount; iCell++) {
		if ((grid.cells[iCell].flag & CELL_FLAG_LIM) > 0) n++;
	}
	return n;
}

void FEM_DG_IMPLICIT::remediateLimCells()
{
	for (int iCell = 0; iCell < grid.cCount; iCell++)
	{
		if (cellIsLim(iCell))
		{
			// пересчитываем по соседям			
			double sRO = 0.0;
			double sRU = 0.0;
			double sRV = 0.0;
			double sRW = 0.0;
			double sRE = 0.0;
			double V = 0.0;
			for (int i = 0; i < grid.cells[iCell].fCount; i++)
			{
				int		iFace = grid.cells[iCell].facesInd[i];
				int		j = grid.faces[iFace].c2;
				if (j == iCell)	{
					//std::swap(j, grid.edges[iEdge].c1); // так нужно еще нормаль поворчивать тогда
					j = grid.faces[iFace].c1;
				}
				if (j >= 0) {
					double  v = grid.cells[j].V;
					V += v;
					sRO += getField(FIELD_RO, j, grid.cells[j].c) * v;
					sRU += getField(FIELD_RU, j, grid.cells[j].c) * v;
					sRV += getField(FIELD_RV, j, grid.cells[j].c) * v;
					sRW += getField(FIELD_RW, j, grid.cells[j].c) * v;
					sRE += getField(FIELD_RE, j, grid.cells[j].c) * v;
				}
			}
			if (V >= TAU*TAU) {
				memset(ro[iCell], 0, sizeof(double) * BASE_FUNC_COUNT);
				memset(ru[iCell], 0, sizeof(double) * BASE_FUNC_COUNT);
				memset(rv[iCell], 0, sizeof(double) * BASE_FUNC_COUNT);
				memset(rw[iCell], 0, sizeof(double) * BASE_FUNC_COUNT);
				memset(re[iCell], 0, sizeof(double) * BASE_FUNC_COUNT);

				ro[iCell][0] = sRO / V;
				ru[iCell][0] = sRU / V;
				rv[iCell][0] = sRV / V;
				rw[iCell][0] = sRW / V;
				re[iCell][0] = sRE / V;
			}
			// после 0x20 итераций пробуем вернуть ячейку в счет
			grid.cells[iCell].flag += 0x010000;
			if (grid.cells[iCell].flag & 0x200000) grid.cells[iCell].flag &= 0x001110;
		}
	}
}


void FEM_DG_IMPLICIT::decCFL()
{
	if (CFL > 0.01) {
		CFL *= 0.75;
		if (CFL < 0.01) CFL = 0.01;
		log(" CFL Number has been decreased : %25.25e \n", CFL);
	}
	//log(" <<< CFL Number has been decreased : %25.25e \n", CFL);
}

void FEM_DG_IMPLICIT::incCFL()
{
	if (CFL < maxCFL) {
		CFL *= scaleCFL;
		if (CFL > maxCFL) CFL = maxCFL;
		log(" CFL Number has been increased : %25.25e \n", CFL);
	}
}

double **FEM_DG_IMPLICIT::allocMtx5()
{
	double		**tempMtx5 = new double*[5];
	for (int i = 0; i < 5; ++i) tempMtx5[i] = new double[5];
	return tempMtx5;
}
void   FEM_DG_IMPLICIT::freeMtx5(double **mtx5)
{
	for (int i = 0; i < 5; ++i)
		delete[] mtx5[i];
	delete[] mtx5;
}
void FEM_DG_IMPLICIT::multMtx5(double **dst5, double **srcA5, double **srcB5)
{
	double sum;
	for (int i = 0; i < 5; ++i)
	{
		for (int j = 0; j < 5; ++j)
		{
			sum = 0;
			for (int k = 0; k < 5; ++k)
				sum += srcA5[i][k] * srcB5[k][j];
			dst5[i][j] = sum;
		}
	}
}
void FEM_DG_IMPLICIT::clearMtx5(double **mtx5)
{
	for (int i = 0; i < 5; ++i)
	for (int j = 0; j < 5; ++j)
		mtx5[i][j] = 0;
}

void FEM_DG_IMPLICIT::multMtxToVal(double **dst, double x, int N)
{
	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < N; ++j)
		{
			dst[i][j] *= x;
		}
	}
}

void FEM_DG_IMPLICIT::fillMtx(double** dst, double x, int N)
{
	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < N; ++j)
		{
			dst[i][j] = x;
		}
	}
}

void FEM_DG_IMPLICIT::eigenValues(double** dst5, double c, double u, double nx, double v, double ny, double w, double nz)
{
	clearMtx5(dst5);
	double un = u * nx + v * ny + w * nz;
	dst5[0][0] = un;
	dst5[1][1] = un;
	dst5[2][2] = un;
	dst5[3][3] = un + c;
	dst5[4][4] = un - c;
}

void FEM_DG_IMPLICIT::eigenValuesAbs(double** dst5, double c, double u, double nx, double v, double ny, double w, double nz)
{
	clearMtx5(dst5);
	double un = u * nx + v * ny + w * nz;
	dst5[0][0] = fabs(un);
	dst5[1][1] = fabs(un);
	dst5[2][2] = fabs(un);
	dst5[3][3] = fabs(un + c);
	dst5[4][4] = fabs(un - c);
}
void FEM_DG_IMPLICIT::rightEigenVector(double **dst5, double c, double u, double nx, double v, double ny, double w, double nz, double H)
{
	double un = u * nx + v * ny + w * nz;
	double q2 = 0.5 * (u * u + v * v + w * w);
	
	dst5[0][0] = nx;
	dst5[0][1] = ny;
	dst5[0][2] = nz;
	dst5[0][3] = 1.0;
	dst5[0][4] = 1.0;

	dst5[1][0] = u * nx;
	dst5[1][1] = u * ny - c * nz;
	dst5[1][2] = u * nz + c * ny;
	dst5[1][3] = u + c * nx;
	dst5[1][4] = u - c * nx;

	dst5[2][0] = v * nx + c * nz;
	dst5[2][1] = v * ny;
	dst5[2][2] = v * nz - c * nx;
	dst5[2][3] = v + c * ny;
	dst5[2][4] = v - c * ny;

	dst5[3][0] = w * nx - c * ny;
	dst5[3][1] = w * ny + c * nx;
	dst5[3][2] = w * nz;
	dst5[3][3] = w + c * nz;
	dst5[3][4] = w - c * nz;

	dst5[4][0] = q2 * nx + c * v * nz - c * w * ny;
	dst5[4][1] = q2 * ny + c * w * nx - c * u * nz;
	dst5[4][2] = q2 * nz + c * u * ny - c * v * nx;
	dst5[4][3] = H + c * un;
	dst5[4][4] = H - c * un;

}
void FEM_DG_IMPLICIT::leftEigenVector(double** dst5, double c, double GAM, double u, double nx, double v, double ny, double w, double nz)
{
	double un = u * nx + v * ny + w * nz;
	double q2 = u * u + v * v + w * w;
	double g1 = GAM - 1.0;
	double c2 = c * c;

	dst5[0][0] = (nx * (c2 - g1 * q2) + c * (w * ny - v * nz)) / c2;
	dst5[0][1] = (nx * g1 * u) / c2;
	dst5[0][2] = (nx * g1 * v + c * nz) / c2;
	dst5[0][3] = (nx * g1 * w - c * ny) / c2;
	dst5[0][4] = -nx * g1 / c2;

	dst5[1][0] = (ny * (c2 - g1 * q2) + c * (u * nz - w * nx)) / c2;
	dst5[1][1] = (ny * g1 * u - c * nz) / c2;
	dst5[1][2] = (ny * g1 * v) / c2;
	dst5[1][3] = (ny * g1 * w + c * nx) / c2;
	dst5[1][4] = -ny * g1 / c2;

	dst5[2][0] = (nz * (c2 - g1 * q2) + c * (v * nx - u * ny)) / c2;
	dst5[2][1] = (nz * g1 * u + c * ny) / c2;
	dst5[2][2] = (nz * g1 * v - c * nx) / c2;
	dst5[2][3] = (nz * g1 * w) / c2;
	dst5[2][4] = -nz * g1 / c2;

	dst5[3][0] = 0.5 * (g1 * q2 - c * un) / c2;
	dst5[3][1] = 0.5 * (c * nx - g1 * u) / c2;
	dst5[3][2] = 0.5 * (c * ny - g1 * v) / c2;
	dst5[3][3] = 0.5 * (c * nz - g1 * w) / c2;
	dst5[3][4] = 0.5 * g1 / c2;

	dst5[4][0] = 0.5 * (g1 * q2 + c * un) / c2;
	dst5[4][1] = 0.5 * (-c * nx - g1 * u) / c2;
	dst5[4][2] = 0.5 * (-c * ny - g1 * v) / c2;
	dst5[4][3] = 0.5 * (-c * nz - g1 * w) / c2;
	dst5[4][4] = 0.5 * g1 / c2;

}

void FEM_DG_IMPLICIT::calcAP(double **dst5, double **rightEgnVecl5, double **egnVal5, double **leftEgnVecl5)
{
	double		**tempMtx5 = allocMtx5();

	//egnVal5 +.
	for (int i = 0; i < 5; ++i)
	for (int j = 0; j < 5; ++j)
		dst5[i][j] = egnVal5[i][j] > 0 ? egnVal5[i][j] : 0;

	multMtx5(tempMtx5, rightEgnVecl5, dst5);
	multMtx5(dst5, tempMtx5, leftEgnVecl5);
	freeMtx5(tempMtx5);
}
void FEM_DG_IMPLICIT::calcAM(double **dst5, double **rightEgnVecl5, double **egnVal5, double **leftEgnVecl5)
{
	double		**tempMtx5 = allocMtx5();

	//egnVal5 -.
	for (int i = 0; i < 5; ++i)
	for (int j = 0; j < 5; ++j)
		dst5[i][j] = egnVal5[i][j] < 0 ? egnVal5[i][j] : 0;

	multMtx5(tempMtx5, rightEgnVecl5, dst5);
	multMtx5(dst5, tempMtx5, leftEgnVecl5);
	//multMtxToVal(dst5, 0.5, 5);
	freeMtx5(tempMtx5);
}

void FEM_DG_IMPLICIT::_calcA(double **dst5, double **rightEgnVecl5, double **egnVal5, double **leftEgnVecl5)
{
	double		**tempMtx5 = allocMtx5();

	multMtx5(tempMtx5, rightEgnVecl5, egnVal5);
	multMtx5(dst5, tempMtx5, leftEgnVecl5);
	freeMtx5(tempMtx5);
}

void FEM_DG_IMPLICIT::calcA(double **dst5, double c, double GAM, double u, double nx, double v, double ny, double w, double nz, double H)
{
	double **rightEgnVecl5 = allocMtx5();
	double **egnVal5 = allocMtx5();
	double **leftEgnVecl5 = allocMtx5();

	eigenValues(egnVal5, c, u, nx, v, ny, w, nz);
	rightEigenVector(rightEgnVecl5, c, u, nx, v, ny, w, nz, H);
	leftEigenVector(leftEgnVecl5, c, GAM, u, nx, v, ny, w, nz);
	_calcA(dst5, rightEgnVecl5, egnVal5, leftEgnVecl5);


	freeMtx5(rightEgnVecl5);
	freeMtx5(egnVal5);
	freeMtx5(leftEgnVecl5);
}

void FEM_DG_IMPLICIT::calcAbsA(double** dst5, double c, double GAM, double u, double nx, double v, double ny, double w, double nz, double H)
{
	double** rightEgnVecl5 = allocMtx5();
	double** egnVal5 = allocMtx5();
	double** leftEgnVecl5 = allocMtx5();

	eigenValuesAbs(egnVal5, c, u, nx, v, ny, w, nz);
	rightEigenVector(rightEgnVecl5, c, u, nx, v, ny, w, nz, H);
	leftEigenVector(leftEgnVecl5, c, GAM, u, nx, v, ny, w, nz);
	_calcA(dst5, rightEgnVecl5, egnVal5, leftEgnVecl5);


	freeMtx5(rightEgnVecl5);
	freeMtx5(egnVal5);
	freeMtx5(leftEgnVecl5);
}

void FEM_DG_IMPLICIT::calcAx(double **dst5, double c, double GAM, double u, double v, double w, double H)
{
	calcA(dst5, c, GAM, u, 1.0, v, 0.0, w, 0.0, H);
}

void FEM_DG_IMPLICIT::calcAy(double **dst5, double c, double GAM, double u, double v, double w, double H)
{
	calcA(dst5, c, GAM, u, 0.0, v, 1.0, w, 0.0, H);
}

void FEM_DG_IMPLICIT::calcAz(double** dst5, double c, double GAM, double u, double v, double w, double H)
{
	calcA(dst5, c, GAM, u, 0.0, v, 0.0, w, 1.0, H);
}

void FEM_DG_IMPLICIT::calcAx_(double **dst4, Param par, double GAM)
{
	double AGAM = GAM - 1.0;
	dst4[0][0] = 0.0;
	dst4[0][1] = 1.0;
	dst4[0][2] = 0.0;
	dst4[0][3] = 0.0;

	dst4[1][0] = -POW_2(par.u) + 0.5*par.U2()*AGAM;
	dst4[1][1] = 2.0*par.u - par.u*AGAM;
	dst4[1][2] = -par.v*AGAM;
	dst4[1][3] = AGAM;

	dst4[2][0] = -par.u*par.v;
	dst4[2][1] = par.v;
	dst4[2][2] = par.u;
	dst4[2][3] = 0.0;

	dst4[3][0] = par.U2()*par.u*AGAM - par.E*par.u*GAM;
	dst4[3][1] = -POW_2(par.u)*AGAM + par.E*GAM - 0.5*par.U2()*AGAM;
	dst4[3][2] = par.u*par.v*AGAM;
	dst4[3][3] = par.u*GAM;
}

void FEM_DG_IMPLICIT::calcAy_(double **dst4, Param par, double GAM)
{
	double AGAM = GAM - 1.0;
	dst4[0][0] = 0.0;
	dst4[0][1] = 0.0;
	dst4[0][2] = 1.0;
	dst4[0][3] = 0.0;

	dst4[1][0] = -par.u*par.v;
	dst4[1][1] = par.v;
	dst4[1][2] = par.u;
	dst4[1][3] = 0.0;

	dst4[2][0] = -POW_2(par.u) + (par.r*par.E + 0.5*par.U2())*AGAM;
	dst4[2][1] = 2.0*par.u + (par.r*par.E - par.u)*AGAM;
	dst4[2][2] = -par.v*AGAM;
	dst4[2][3] = AGAM;

	dst4[3][0] = par.U2()*par.u*AGAM - par.E*par.u*GAM;
	dst4[3][1] = -POW_2(par.u)*AGAM + par.E*GAM - 0.5*par.U2()*AGAM;
	dst4[3][2] = par.u*par.v*AGAM;
	dst4[3][3] = par.u*GAM;
}

void FEM_DG_IMPLICIT::calcRoeAverage(Param& average, Param pL, Param pR, double GAM, Vector n)
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

	roe_orig(average.r, average.e, average.p, average.u, average.v, average.w,
		pL.r, pL.p, pL.u, pL.v, pL.w,
		pR.r, pR.p, pR.u, pR.v, pL.w, GAM);

	average.cz = sqrt(GAM*average.p / average.r);
	average.E = average.e + 0.5 * average.U2();
}


void FEM_DG_IMPLICIT::consToPar(double fRO, double fRU, double fRV, double fRW, double fRE, Param& par)
{
	par.r = fRO;
	par.u = fRU / fRO;
	par.v = fRV / fRO;
	par.w = fRW / fRO;
	par.E = fRE / fRO;
	par.e = par.E - par.U2()*0.5;
}

void FEM_DG_IMPLICIT::calcGeomPar()
{
	double minDistance = 1.0e-7 / L_;

#define MIN_DISTANCE( _X_ ) { if ( _X_ < minDistance ) _X_ = minDistance; }

	//Patch& rPatch = Patch_GetAt(0);  // внутренниие грани

	for (int iFace = 0; iFace < grid.fCount; ++iFace)
	{
		int c1 = grid.faces[iFace].c1;
		int c2 = grid.faces[iFace].c2;

		if (c2 > -1) {
			Point ptP = grid.cells[c1].c;
			Point ptE = grid.cells[c2].c;
			Point ptG = grid.faces[iFace].c;

			Vector n = grid.faces[iFace].n;

			Vector vPG = ptG;  vPG -= ptP;		// vPG = ptG - ptP 
			double   fPG = _mag_(vPG);
			double   fLP_ = ::fabs(scalar_prod(n, vPG));  MIN_DISTANCE(fLP_);
			faceRP[iFace] = vPG;
			faceLP_[iFace] = fLP_;

			Vector vEG = ptG;  vEG -= ptE;		// vEG = ptG - ptE 
			double   fEG = _mag_(vEG);
			double   fLE_ = ::fabs(scalar_prod(n, vEG));  MIN_DISTANCE(fLE_);
			faceRE[iFace] = vEG;
			faceLE_[iFace] = fLE_;

			//m_pFaceW[iFace] = fEG / (fPG + fEG);	// вес интерполяции параметров на грань
		}
		else {
			Point ptP = grid.cells[c1].c;
			Point ptG = grid.faces[iFace].c;

			Vector n = grid.faces[iFace].n;	// получили нормаль к грани

			Vector vPG = ptG;  vPG -= ptP;		// vPG = ptG - ptP 
			faceRP[iFace] = vPG;
			faceLP_[iFace] = ::fabs(scalar_prod(n, vPG));	MIN_DISTANCE(faceLP_[iFace]);
		}

	} // for( iFace )

	//m_fSumVolume = 0.0;

	for (int iCell = 0; iCell < grid.cCount; iCell++)	// цикл по всем ячейкам
	{
		//m_fSumVolume += Cell_GetVolume(iCell);
		Y[iCell] = 1.0e10;
	} // for(iCell)

	int wallCount = 0;
	//static INDEX iFirstStep = 1;
	for (int iFace = 0; iFace < grid.fCount; ++iFace)
	{
		if (!grid.faces[iFace].bnd) continue;
		if (grid.faces[iFace].bnd->faceType == CFDBoundary::TYPE_ID_WALL_NO_SLIP || grid.faces[iFace].bnd->faceType == CFDBoundary::TYPE_ID_WALL_SLIP && !VISCOSITY)
		{
			wallCount++;
		}
	}

	// Расчёт расстояния до стенки
	if (TURBULENCE)
	{
		//Proc_AllReduce(&iWallCount, MPI_SUM); // теперь iWallCount - общее кол-во точек на стенке по всем процессорам

		//m_iWallCount = iWallCount;

		if (wallCount == 0) return;	// в задаче нет стенок с прилипанием

		log("Calculating Wall Distances .\n");	double start = MPI_Wtime();

		Point* pWall = new Point[wallCount];	 ::memset(pWall, 0L, wallCount * sizeof(Point));
		Vector* pWallN = new Vector[wallCount];	 ::memset(pWallN, 0L, wallCount * sizeof(Vector));


		// Определение угловых вершин  pNodeF[ iNode ] >= 100
		// в прилежащих гранях vN = 0  и расстояние до стенки вычисляется по KD дереву

		unsigned char* pNodeF = new unsigned char[grid.nCountEx];	::memset(pNodeF, 0L, grid.nCountEx * sizeof(unsigned char));
		Vector* pNodeN = new Vector[grid.nCountEx];	::memset(pNodeN, 0L, grid.nCountEx * sizeof(Vector));

		for (int iFace = 0; iFace < grid.fCount; iFace++)		// цикл по граням
		{
			int c1 = grid.faces[iFace].c1;
			if (!grid.faces[iFace].bnd) continue;
			if (grid.faces[iFace].bnd->faceType == CFDBoundary::TYPE_ID_WALL_NO_SLIP || grid.faces[iFace].bnd->faceType == CFDBoundary::TYPE_ID_WALL_SLIP && !VISCOSITY)
			{
					Y[c1] = faceLP_[iFace];

					Vector n = grid.faces[iFace].n;
					//Face_GetNodeIndex(iFace, &rNodeList);
					for (int i = 0; i < grid.faces[iFace].nCount; i++) {
						int iNode = grid.faces[iFace].nodesInd[i];
						if (pNodeF[iNode] && pNodeF[iNode] < 100) {
							double Cos = scalar_prod(pNodeN[iNode], n);
							double Cos2 = _sign_(Cos) * Cos * Cos / _mag_2_(pNodeN[iNode]);
							if (Cos2 < 0.64) pNodeF[iNode] = 100;		// если угол > 36 градусов, то устанавливаем признак 
						} // if
						pNodeN[iNode] += n;
						pNodeF[iNode] += 1;
					} // for( n )
				//} // for( iFace )
			} // if
		} // for( iFace )

		int iCount = 0;

		for (int iFace = 0; iFace < grid.fCount; iFace++)		// цикл по граням
		{
			int c1 = grid.faces[iFace].c1;
			if (!grid.faces[iFace].bnd) continue;
			if (grid.faces[iFace].bnd->faceType == CFDBoundary::TYPE_ID_WALL_NO_SLIP || grid.faces[iFace].bnd->faceType == CFDBoundary::TYPE_ID_WALL_SLIP && !VISCOSITY)
			{
				pWall[iCount] = grid.faces[iFace].c;
				pWallN[iCount] = grid.faces[iFace].n;
				for (int i = 0; i < grid.faces[iFace].nCount; i++) {
					int iNode = grid.faces[iFace].nodesInd[i];
					if (pNodeF[iNode] >= 100) { pWallN[iCount] = 0.0;  break; } // if
				} // for( n )
			} // if
		} // for( iFace )

		delete[] pNodeF;
		delete[] pNodeN;

		{ {
				int out = 0;

				int nPts = wallCount;
				ANNpointArray dataPts = annAllocPts(nPts, 3);
				for (int i = 0; i < wallCount; i++)
				{
					ANNpoint dp = dataPts[i]; Point wp = pWall[i];
					dp[0] = wp.x; dp[1] = wp.y; dp[2] = wp.z;
				} // for(i)

				ANNkd_tree* kdTree = new ANNkd_tree(dataPts, nPts, 3);

				ANNpoint      queryPt = annAllocPt(3);
				ANNidxArray   nnIdx = new ANNidx;
				ANNdistArray  dists = new ANNdist;

				//INDEX_LIST rFaceList;
				for (int iCell = 0; iCell < grid.cCount; iCell++)
				{
					if (!(++out & 0x3fff)) log(".");

					Point ptCell = grid.cells[iCell].c;

					queryPt[0] = ptCell.x;
					queryPt[1] = ptCell.y;
					queryPt[2] = ptCell.z;

					kdTree->annkSearch(queryPt, 1, nnIdx, dists);

					int index = nnIdx[0];

					//if (m_pWallIndex) m_pWallIndex[iCell] = index;

					Vector vN = pWallN[index];

					//Cell_GetFaceIndex(iCell, &rFaceList);
					//for (int nf = 0; nf < grid.cells[iCell].fCount; ++nf)
					//{
					//	int iFace = grid.cells[iCell].facesInd[nf];
					//	if (_pow_2_(scalar_prod(vN, grid.faces[iFace].n)) > 0.75)
					//		//if (Face_GetCell_2(iFace)) m_pFaceFlag[iFace] |= FLG_FACE_ORIENT;
					//} // for(nf) 

					double& fY = Y[iCell];
					if (fY < 1.0e10) continue; // для ячеек на стенке с прилипанием растояние уже расчитано

					double fD = sqrt(dists[0]);

					//			if ( fD > 1.0 || _mag_2_( vN ) < 0.1 ) 
					if (_mag_2_(vN) < 0.1)
					{
						fY = fD;
					}
					else {
						double fDN = fabs(scalar_prod(pWall[index] - ptCell, vN));
						fY = fDN;
						//				if ( fD > 2.0 * fDN ) fY = fD;
					} // if

					MIN_DISTANCE(fY);

				} // for(iCell)

				annDeallocPts(dataPts);
				delete kdTree;
				annDeallocPt(queryPt);
				delete nnIdx;
				delete dists;
				annClose();
			}}

		delete[] pWall;
		delete[] pWallN;


		double finish = MPI_Wtime();
		log(" \n\n cpu = %10.3e \n\n", finish - start);

	}
}

void FEM_DG_IMPLICIT::calcGrad()
{
	for (int iFace = 0; iFace < grid.fCount; iFace++)
	{
		int c1 = grid.faces[iFace].c1;
		int c2 = grid.faces[iFace].c2;

		Param pL, pR;
		if (grid.faces[iFace].type == Face::TYPE_INNER)
		{
			int c1 = grid.edges[iFace].c1;
			int c2 = grid.edges[iFace].c2;
			convertConsToPar(c1, pL);
			convertConsToPar(c2, pR);
		}
		else {
			int c1 = grid.edges[iFace].c1;
			convertConsToPar(c1, pL);
			boundaryCond(iFace, pL, pR);
		}

		Vector n = grid.faces[iFace].n;
		double S = grid.faces[iFace].S;

		gradU[c1].x += 0.5 * (pL.u + pR.u) * n.x * S;
		gradU[c1].y += 0.5 * (pL.u + pR.u) * n.y * S;
		gradU[c1].z += 0.5 * (pL.u + pR.u) * n.z * S;
		
		gradV[c1].x += 0.5 * (pL.v + pR.v) * n.x * S;
		gradV[c1].y += 0.5 * (pL.v + pR.v) * n.y * S;
		gradV[c1].z += 0.5 * (pL.v + pR.v) * n.z * S;

		gradW[c1].x += 0.5 * (pL.w + pR.w) * n.x * S;
		gradW[c1].y += 0.5 * (pL.w + pR.w) * n.y * S;
		gradW[c1].z += 0.5 * (pL.w + pR.w) * n.z * S;

		if (c2 > -1)
		{
			gradU[c1].x -= 0.5 * (pL.u + pR.u) * n.x * S;
			gradU[c1].y -= 0.5 * (pL.u + pR.u) * n.y * S;
			gradU[c1].z -= 0.5 * (pL.u + pR.u) * n.z * S;

			gradV[c1].x -= 0.5 * (pL.v + pR.v) * n.x * S;
			gradV[c1].y -= 0.5 * (pL.v + pR.v) * n.y * S;
			gradV[c1].z -= 0.5 * (pL.v + pR.v) * n.z * S;

			gradW[c1].x -= 0.5 * (pL.w + pR.w) * n.x * S;
			gradW[c1].y -= 0.5 * (pL.w + pR.w) * n.y * S;
			gradW[c1].z -= 0.5 * (pL.w + pR.w) * n.z * S;
		}

	}
	for (int iCell = 0; iCell < grid.cCount; iCell++)
	{
		double V = grid.cells[iCell].V;
		gradU[iCell] /= V;
		gradV[iCell] /= V;
		gradW[iCell] /= V;
	}
}
/*
void FEM_DG_IMPLICIT::calcTurbulent()
{
	//CONST_TURBULENT_SA& rConst = m_rConstTurbulentSA;	// ссылка на параметры модели турбулентности

	REAL* pNT = m_pTP[0];
	REAL* pNTg = m_pTPg[0];
	REAL* pBn = m_pBnTP[0];
	REAL* pNTn = m_pTPn[0];

	INDEX	ID_NT = ID_TP[0];

	::memset(pBn, 0, Cell_GetCount() * sizeof(REAL));

	// занулили матрицу
	//Solver_ZeroMatrix(); 		//*********************************************************

	// вычисляем правую часть и коэффициент на диагонали матрицы

	double Cw1 = Cb1 / _pow_2_(K) + (1.0 + Cb2) / SigmaNT;

	for (int iCell = 0; iCell < grid.cCount; iCell++)
	{
		double R = ro[iCell][0];			// плотность
		double y = Y[iCell];			// растояние до стенки
		double NT = pNTg[iCell];			// турбулентный параметр NuT

		double Volume = grid.cells[iCell].V;// Cell_GetVolume(iCell);

		double Kapa = R * NT / m_pML[iCell];
		double fKapa3 = _pow_3_(Kapa);		// рабочая переменная

		double fFv1 = fKapa3 / (fKapa3 + _pow_3_(rConst.fCv1));

		if (m_eWallTreatment == HIGH_Y_PLUS) fFv1 = 1.0;

		double Fv2 = 1.0 - fKapa / (1.0 + fKapa * fFv1);

		double Wxy2 = _pow_2_(gradU->y - gradV->x);
		double Wxz2 = _pow_2_(gradU->z - gradW->x);
		double Wyz2 = _pow_2_(gradV->z - gradW->y);

		double OmegaZ = ::sqrt(Wxy2 + Wxz2 + Wyz2);

		double Sxx2 = _pow_2_(gradU->x);
		double Syy2 = _pow_2_(gradV->y);
		double Szz2 = _pow_2_(gradW->z);

		double Sxy2 = _pow_2_(0.5 * (gradU->y + gradV->x));
		double Sxz2 = _pow_2_(0.5 * (gradU->z + gradW->x));
		double Syz2 = _pow_2_(0.5 * (gradV->z + gradW->y));

		double SZ = ::sqrt(2.0 * (Sxx2 + Syy2 + Szz2 + 2.0 * (Sxy2 + Sxz2 + Syz2)));

		OmegaZ += Cprod * _min_(0.0, SZ - OmegaZ);

		double KY2 = _pow_2_(K * y);	// рабочая переменная

		double S_ = OmegaZ + NT * Fv2 / KY2;



		//if ( fS_ < 1.0e-16 ) fS_ = 1.0e-16; // Закомментировано согласно рекомендациям СПбГТУ



		double Fw = 0.0;
		if (NT > 3.0 * fabs(S_) * KY2)
		{
			Fw = ::pow(1.0 + _pow_6_(Cw3), 1.0 / 6.0);
		}
		else {
			double Gama = NT / ((S_ + 1e-16) * KY2);

			double G = Gama + Cw2 * (_pow_6_(Gama) - Gama);
			double G6 = _pow_6_(G);			// рабочая переменная
			double C6 = _pow_6_(Cw3);	// рабочая переменная

			Fw = G * pow((1.0 + C6) / (G6 + C6), 1.0 / 6.0);
		} // if

		double HNT = 0.0;
		double BNT = 0.0;

		double Ft2 = Ct3 * exp(-Ct4 * _pow_2_(Kapa));

		//		double fA1 = fR * rConst.fCb1 * ( 1.0 - fFt2 ) * fS_;
		double A1 = R * Cb1 * S_;

		if (A1 >= 0.0)
		{
			HNT += A1 * NT;
			//		} else {
			//			fBNT += fA1;
		} // if

		double A2 = -R * (Cw1 * Fw - Cb1 * Ft2 / _pow_2_(K)) * NT / (y * y);

		if (A2 <= 0.0)
		{
			//			fHNT += fA2 * fNT;
			//		} else {
			BNT += A2;
		} // if

		if (grid.cells[iCell].flag & CELL_FLAG_BAD)
		{
			HNT = 0.0;
			BNT = 0.0;
		} // if

		double Tau_ = Volume / cTau[iCell];
		double TauPhys_ = Volume / TAU;

		double Eps0 = 1.0, Eps1 = 1.0,  Eps2 = 0.0; // коэффициенты аппроксимации по времени
		if (STEADY)
		{
			double K = Tau_ * Eps0 * ro[iCell][0] - BNT * Volume;

			Solver_SetValToDiag(iCell, K); 		//*********************************************************

			pBn[iCell] = Eps1 * Tau_ * ro[iCell][0] * pNT[iCell] + HNT * Volume;
			//if (m_iTOrder == 2)
			//{
			//	double* pNTn_ = m_pTPn_[0];
			//	pBn[iCell] -= m_fEps2 * fTau_ * m_pRn_[iCell] * pNTn_[iCell];
			//} // if m_iTOrder == 2

		}
		else { // m_iUnsteady != 0

			double K = (Tau_ + Eps0 * TauPhys_) * ro[iCell][0] - BNT * Volume;

			Solver_SetValToDiag(iCell, K); 		//*********************************************************

			pBn[iCell] = Tau_ * ro[iCell][0] * pNT[iCell] + TauPhys_ * ro[iCell][0] * pNTn[iCell] + HNT * Volume;

			//if (m_iTOrder == 1)
			//{
			//	pBn[iCell] = fTau_ * m_pR[iCell] * pNT[iCell] + fTauPhys_ * m_pRn[iCell] * pNTn[iCell] + fHNT * fVolume;
			//} // if

			//if (m_iTOrder == 2)
			//{
			//	REAL* pNTn_ = m_pTPn_[0];
			//	pBn[iCell] = fTau_ * m_pR[iCell] * pNT[iCell] +
			//		fTauPhys_ * (m_fEps1 * m_pRn[iCell] * pNTn[iCell] -
			//			m_fEps2 * m_pRn_[iCell] * pNTn_[iCell]) + fHNT * fVolume;
			//} // if
		} // if

		m_pFM[iCell] = fFv1;
	} // for( iCell )




	// Расчёт диффузионных и конвективных потоков
	for (int iFace = 0; iFace < grid.fCount; iFace++) {
		int c1 = grid.faces[iFace].c1;
		int c2 = grid.faces[iFace].c2;

		if (c2 > -1) {
			INDEX  iP = Face_GetCell_1(iFace);
			INDEX  iE = Face_GetCell_2(iFace);

			double fML = 0.5 * (m_pML[iP] + m_pML[iE]);		// молекулярная вязкость

			double fMT = 0.5 * (m_pRg[iP] * pNTg[iP] + m_pRg[iE] * pNTg[iE]);		// переносится ML + Ro * NTg а не ML + MT
			double fNT = 0.5 * (pNTg[iP] + pNTg[iE]);		// турбулентный параметр

			double fDS = Face_GetSquare(iFace) / (m_pFaceLP_[iFace] + m_pFaceLE_[iFace]);

			double fMf = m_pFaceRUN[iFace];

			double fM_max = 0.0;
			double fM_min = 0.0;
			if (fMf > 0.0) fM_max = fMf;
			else			 fM_min = fMf;

			REAL fKP = fDS * _max_(0.0, ((fML + fMT) + rConst.fCb2 * m_pRg[iP] * (fNT - pNTg[iP])) / rConst.fSigmaNT);
			REAL fKE = fDS * _max_(0.0, ((fML + fMT) + rConst.fCb2 * m_pRg[iE] * (fNT - pNTg[iE])) / rConst.fSigmaNT);

			if (rPatch.GetType() == PATCH_CYCLIC)
			{
				Solver_AddValToCyclic(iFace, fKP + fM_max, -fKP + fM_min);  	//*********************************************************
			}
			else {
				Solver_AddValToDU(iFace, fKP + fM_max, -fKP + fM_min); 		//*********************************************************
				Solver_AddValToDL(iFace, fKE - fM_min, -fKE - fM_max); 		//*********************************************************
			} // if
		}
		else {
			if (rPatch.GetType() == PATCH_CYCLIC) continue;		// пропустили CYCLIC
			if (rPatch.GetType() == PATCH_EXCHANGE) continue;		// пропустили фиктивные грани
			if (rPatch.GetType() == PATCH_WALL_SLIP) continue;	// конвективные и диффузионные потоки = 0
			if (rPatch.GetType() == PATCH_SYMMETRY) continue;	// конвективные и диффузионные потоки = 0

			if (rPatch.GetType() == PATCH_WALL_NO_SLIP) // стенка с прилипанием
			{
				if (m_eWallTreatment == LOW_Y_PLUS)
				{
					for (INDEX iFace = rPatch.GetFaceStart(); iFace < rPatch.GetFaceFinish(); ++iFace)
					{
						INDEX  iP = Face_GetCell_1(iFace);
						REAL   fDS = Face_GetSquare(iFace) / m_pFaceLP_[iFace];

						REAL fK = m_pML[iP] / rConst.fSigmaNT * fDS;

						Solver_AddValToDiag(iP, fK); 		//*********************************************************
					} // for( iFace )

				} // if

				// для стенки с прилипанием и ( m_eWallTreatment == HIGH_Y_PLUS ) в самый послелний момент

			}
			else {	//  if ( rPatch.GetType() != PATCH_WALL_NO_SLIP )

				PAR_LIST& rPar = rPatch.GetParam();		// параметры патча

				for (INDEX iFace = rPatch.GetFaceStart(); iFace < rPatch.GetFaceFinish(); ++iFace)
				{
					INDEX  iP = Face_GetCell_1(iFace);
					REAL   fDS = Face_GetSquare(iFace) / m_pFaceLP_[iFace];

					REAL fNTbound = rPatch.BC_PARAM(ID_NT, iFace);		// параметр на границе

					REAL fML = m_pML[iP];		// молекулярная вязкость
					REAL fMT = m_pRg[iP] * pNTg[iP];		// переносится ML + Ro * NTg а не ML + MT

					REAL fNT = pNTg[iP];		// турбулентный параметр

					REAL fMf = m_pFaceRUN[iFace];

					if (fMf > 0.0)
					{
						Solver_AddValToDiag(iP, fMf); 	//*********************************************************

					}
					else {
						REAL fQN = fMf * fNTbound;	// конвективный поток

						REAL fK = fDS * _max_(0.0, ((fML + fMT) + rConst.fCb2 * m_pRg[iP] * (fNTbound - fNT)) / rConst.fSigmaNT);

						Solver_AddValToDiag(iP, fK); 		//*********************************************************

						pBn[iP] += fK * fNTbound - fQN;
					} // if ( fMf > 0.0 )
				} // for( iFace )
			} // if
		}
	}

	// для HIGH_Y_PLUS и PATCH_WALL_NO_SLIP в самый послелний момент

	if (m_iSUB_SONIC) Solver_RelaxMatrix(0.9, m_pTP[0], m_pBnTP[0]); 	//*********************************************************
	else				Solver_RelaxMatrix(0.7, m_pTP[0], m_pBnTP[0]); 	//*********************************************************


	BadCellMatrix(ID_NT);

	m_pTP_iter[0] = Solver_Solve(m_pDTP[0], m_pBnTP[0]); 				//*********************************************************

	if (m_bSolverError) return;

	for (int iCell = 0; iCell < grid.cCount; iCell++)
	{
		double R = ro[iCell][0];		// плотность
		double ML = m_pML[iCell];		// молекулярная вязкость

		if (!_finite(m_pDTP[0][iCell]))
		{
			m_pDTP[0][iCell] = m_pDTP[0][iCell - 1];
		} // if

		double NTg = TP_relax * m_pDTP[0][iCell] + (1.0 - TP_relax) * pNT[iCell];

		double MLR = ML / R;

		double MT = fR * m_pFM[iCell] * fNTg;

		MT = _min_(MT, 1.0e+5 * ML);
		NTg = MT / (R * m_pFM[iCell] + 1e-20);
		NTg = _max_(NTg, 1.0e-12);

		m_pMT[iCell] = MT;
		m_pKT[iCell] = 0.0;

		NTg[iCell] = NTg;
	} // for( iCell )

	//Proc_Exchange(m_pMT);
	//Proc_Exchange(pNTg);
}
*/

void FEM_DG_IMPLICIT::calcMatrWithTau()
{
	for (int iCell = 0; iCell < grid.cCount; iCell++) {

		fillMtx(matrBig, 0.0, matrDim);
		

		for (int i = 0; i < BASE_FUNC_COUNT; i++) {
			for (int j = 0; j < BASE_FUNC_COUNT; j++) {
				matrSmall[i][j] = matrA[iCell][i][j] / cTau[iCell];
			}
		}

		for (int ii = 0; ii < FIELD_COUNT; ii++) {
			addSmallMatrToBigMatr(matrBig, matrSmall, ii, ii);
		}

		solverMtx->addMatrElement(iCell, iCell, matrBig);

	}
}

void FEM_DG_IMPLICIT::calcIntegral()
{
	double fRO, fRU, fRV, fRW, fRE;
	double FR, u, v, w, FE;
	double **mx = allocMtx5();
	double **my = allocMtx5();
	double **mz = allocMtx5();

	for (int iCell = 0; iCell < grid.cCount; iCell++) {
		
		fillMtx(matrBig, 0.0, matrDim);
		
		for (int iGP = 0; iGP < GP_CELL_COUNT; iGP++) {
			Point& p = cellGP[iCell][iGP];
			double w = cellGW[iCell][iGP];
			getFields(fRO, fRU, fRV, fRW, fRE, iCell, p);
			Param par;
			consToPar(fRO, fRU, fRV, fRW, fRE, par);
			Material& mat = getMaterial(iCell);
			mat.URS(par, 0); // p=p(r,e)
			double H = par.E + par.p / par.r;
			calcA(mx, par.cz, getGAM(iCell), par.u, 1.0, par.v, 0.0, par.w, 0.0, H);
			calcA(my, par.cz, getGAM(iCell), par.u, 0.0, par.v, 1.0, par.w, 0.0, H);
			calcA(mz, par.cz, getGAM(iCell), par.u, 0.0, par.v, 0.0, par.w, 1.0, H);
			//calcAx_(my, par, getGAM(iCell));
			for (int i = 0; i < FIELD_COUNT; i++) {
				for (int j = 0; j < FIELD_COUNT; j++) {
					for (int ii = 0; ii < BASE_FUNC_COUNT; ii++) {
						for (int jj = 0; jj < BASE_FUNC_COUNT; jj++) {
							matrSmall[ii][jj] = -mx[i][j] * getDfDx(ii, iCell, p);
							matrSmall[ii][jj] -= my[i][j] * getDfDy(ii, iCell, p);
							matrSmall[ii][jj] *= getF(jj, iCell, p);
							matrSmall[ii][jj] *= w;
						}
					}
					addSmallMatrToBigMatr(matrBig, matrSmall, i, j);
				}
			}

		}
		
		multMtxToVal(matrBig, cellJ[iCell] * SIGMA, matrDim);
		
		solverMtx->addMatrElement(iCell, iCell, matrBig);
	}

	freeMtx5(mx);
	freeMtx5(my);
	freeMtx5(mz);
}

void FEM_DG_IMPLICIT::calcMatrFlux()
{
	double **eigenMtx5, **rEigenVector5, **lEigenVector5;
	double **Amtx5P, **Amtx5M, **mtx5;
	
	mtx5 = allocMtx5();
	eigenMtx5 = allocMtx5();
	rEigenVector5 = allocMtx5();
	lEigenVector5 = allocMtx5();
	Amtx5P = allocMtx5();
	Amtx5M = allocMtx5();


	double fRO1, fRU1, fRV1, fRW1, fRE1, fRO2, fRU2, fRV2, fRW2, fRE2;
	Param par1, par2;

	//int mSize = BASE_FUNC_COUNT * 4;

	for (int iFace = 0; iFace < grid.fCount; iFace++) {
		Face& face = grid.faces[iFace];
		Vector&	n = grid.faces[iFace].n;
		int c1 = face.c1;
		int c2 = face.c2;
		
		if (c2 >= 0) {
			
			fillMtx(matrBig, 0.0, matrDim); 
			fillMtx(matrBig2, 0.0, matrDim); 
			
			for (int iGP = 0; iGP < GP_FACE_COUNT; iGP++) {
				Point& gp = faceGP[iFace][iGP];
				double gw = faceGW[iFace][iGP];
				getFields(fRO1, fRU1, fRV1, fRW1, fRE1, c1, gp);
				consToPar(fRO1, fRU1, fRV1, fRW1, fRE1, par1);
				Material& mat1 = getMaterial(c1);
				mat1.URS(par1, 0); // p=p(r,e)
				
				getFields(fRO2, fRU2, fRV2, fRW2, fRE2, c2, gp); 
				consToPar(fRO2, fRU2, fRV2, fRW2, fRE2, par2);
				Material& mat2 = getMaterial(c2);
				mat2.URS(par2, 0); // p=p(r,e)

				Param pF;
				double GAM = getGAM(c1);
				calcFluxFields(pF, par1, par2, n, GAM);
				//calcRoeAverage(average, par1, par2, getGAM(c1), n);
				
				double q2 = 0.5 * pF.U2();
				double H = pF.E + q2 + pF.p / pF.r;
				eigenValues(eigenMtx5, pF.cz, pF.u, n.x, pF.v, n.y, pF.w, n.z);
				rightEigenVector(rEigenVector5, pF.cz, pF.u, n.x, pF.v, n.y, pF.w, n.z, H);
				leftEigenVector(lEigenVector5, pF.cz, getGAM(c1), pF.u, n.x, pF.v, n.y, pF.w, n.z);
				calcAP(Amtx5P, rEigenVector5, eigenMtx5, lEigenVector5);
				calcAM(Amtx5M, rEigenVector5, eigenMtx5, lEigenVector5);
				for (int i = 0; i < FIELD_COUNT; i++) {
					for (int j = 0; j < FIELD_COUNT; j++) {
						for (int ii = 0; ii < BASE_FUNC_COUNT; ii++) {
							for (int jj = 0; jj < BASE_FUNC_COUNT; jj++) {
								matrSmall[ii][jj]  = Amtx5P[i][j] * getF(jj, c1, gp) * getF(ii, c1, gp) * gw;
								matrSmall2[ii][jj] = Amtx5M[i][j] * getF(jj, c2, gp) * getF(ii, c1, gp) * gw;
							}
						}
						addSmallMatrToBigMatr(matrBig,  matrSmall, i, j);
						addSmallMatrToBigMatr(matrBig2, matrSmall2, i, j);
					}
				}

			}
			
			multMtxToVal(matrBig, faceJ[iFace] * SIGMA, matrDim);
			multMtxToVal(matrBig2, faceJ[iFace] * SIGMA, matrDim);

			solverMtx->addMatrElement(c1, c1, matrBig);
			solverMtx->addMatrElement(c1, c2, matrBig2);


			
			fillMtx(matrBig, 0.0, matrDim);
			fillMtx(matrBig2, 0.0, matrDim);

			for (int iGP = 0; iGP < GP_FACE_COUNT; iGP++) {
				Point& gp = faceGP[iFace][iGP];
				double gw = faceGW[iFace][iGP];
				getFields(fRO1, fRU1, fRV1, fRW1, fRE1, c1, gp);
				consToPar(fRO1, fRU1, fRV1, fRW1, fRE1, par1);
				Material& mat1 = getMaterial(c1);
				mat1.URS(par1, 0); // p=p(r,e)

				getFields(fRO2, fRU2, fRV2, fRW2, fRE2, c2, gp);
				consToPar(fRO2, fRU2, fRV2, fRW2, fRE2, par2);
				Material& mat2 = getMaterial(c2);
				mat2.URS(par2, 0); // p=p(r,e)

				Param pF;
				double GAM = getGAM(c1);
				calcFluxFields(pF, par1, par2, n, GAM);
				//calcRoeAverage(average, par1, par2, getGAM(c1), n);
				double q2 = 0.5 * pF.U2();
				double H = pF.E + q2 + pF.p / pF.r;

				eigenValues(eigenMtx5, pF.cz, pF.u, -n.x, pF.v, -n.y, pF.w, -n.z);
				rightEigenVector(rEigenVector5, pF.cz, pF.u, -n.x, pF.v, -n.y, pF.w, -n.z, H);
				leftEigenVector(lEigenVector5, pF.cz, getGAM(c1), pF.u, -n.x, pF.v, -n.y, pF.w, -n.z);
				calcAP(Amtx5P, rEigenVector5, eigenMtx5, lEigenVector5);
				calcAM(Amtx5M, rEigenVector5, eigenMtx5, lEigenVector5);
				for (int i = 0; i < FIELD_COUNT; i++) {
					for (int j = 0; j < FIELD_COUNT; j++) {
						for (int ii = 0; ii < BASE_FUNC_COUNT; ii++) {
							for (int jj = 0; jj < BASE_FUNC_COUNT; jj++) {
								matrSmall[ii][jj]  = Amtx5P[i][j] * getF(jj, c2, gp) * getF(ii, c2, gp) * gw;
								matrSmall2[ii][jj] = Amtx5M[i][j] * getF(jj, c1, gp) * getF(ii, c2, gp) * gw;
							}
						}
						addSmallMatrToBigMatr(matrBig, matrSmall, i, j);
						addSmallMatrToBigMatr(matrBig2, matrSmall2, i, j);
					}
				}

			}

			multMtxToVal(matrBig, faceJ[iFace] * SIGMA, matrDim);
			multMtxToVal(matrBig2, faceJ[iFace] * SIGMA, matrDim);

			solverMtx->addMatrElement(c2, c2, matrBig);
			solverMtx->addMatrElement(c2, c1, matrBig2);

		}
		else {

			fillMtx(matrBig, 0.0, matrDim);

			for (int iGP = 0; iGP < GP_FACE_COUNT; iGP++) {
				Point& gp = faceGP[iFace][iGP];
				double gw = faceGW[iFace][iGP];
				getFields(fRO1, fRU1, fRV1, fRW1, fRE1, c1, gp);
				consToPar(fRO1, fRU1, fRV1, fRW1, fRE1, par1);
				Material& mat1 = getMaterial(c1);
                mat1.URS(par1, 0); // p=p(r,e)
                mat1.URS(par1, 1); // T=T(e)

				boundaryCond(iFace, par1, par2);

				Param average;
				calcRoeAverage(average, par1, par2, getGAM(c1), n);

				double q2 = 0.5 * average.U2();
				double H = average.E + q2 + average.p / average.r;

				eigenValues(eigenMtx5, average.cz, average.u, n.x, average.v, n.y, average.w, n.z);
				rightEigenVector(rEigenVector5, average.cz, average.u, n.x, average.v, n.y, average.w, n.z, H);
				leftEigenVector(lEigenVector5, average.cz, getGAM(c1), average.u, n.x, average.v, n.y, average.w, n.z);
				calcAP(Amtx5P, rEigenVector5, eigenMtx5, lEigenVector5);
				for (int i = 0; i < FIELD_COUNT; i++) {
					for (int j = 0; j < FIELD_COUNT; j++) {
						for (int ii = 0; ii < BASE_FUNC_COUNT; ii++) {
							for (int jj = 0; jj < BASE_FUNC_COUNT; jj++) {
								matrSmall[ii][jj] = Amtx5P[i][j] * getF(jj, c1, gp) * getF(ii, c1, gp) * gw;
							}
						}
						addSmallMatrToBigMatr(matrBig, matrSmall, i, j);
					}
				}

			}
			multMtxToVal(matrBig, faceJ[iFace] * SIGMA, matrDim);

			solverMtx->addMatrElement(c1, c1, matrBig);

		}
	}
	freeMtx5(eigenMtx5);
	freeMtx5(rEigenVector5);
	freeMtx5(lEigenVector5);
	freeMtx5(Amtx5P);
	freeMtx5(Amtx5M);
	freeMtx5(mtx5);
}

void FEM_DG_IMPLICIT::calcRHS()
{
	/* volume integral */
	for (int iCell = 0; iCell < grid.cCount; iCell++) {
		memset(tmpArr, 0, sizeof(double)*matrDim);
		for (int iBF = 0; iBF < BASE_FUNC_COUNT; iBF++) {
			double sRO = 0.0;
			double sRU = 0.0;
			double sRV = 0.0;
			double sRW = 0.0;
			double sRE = 0.0;
			for (int iGP = 0; iGP < GP_CELL_COUNT; iGP++) {
				Point& gp = cellGP[iCell][iGP];
				double gw = cellGW[iCell][iGP];
				double fRO, fRU, fRV, fRW, fRE;
				getFields(fRO, fRU, fRV, fRW, fRE, iCell, gp);
				Param par;
				consToPar(fRO, fRU, fRV, fRW, fRE, par);
				Material& mat = getMaterial(iCell);
				mat.URS(par, 0); // p=p(r,e)

				double FR = par.r * par.u;
				double FU = FR * par.u + par.p;
				double FV = FR * par.v;
				double FW = FR * par.w;
				double FE = (fRE + par.p)*par.u;

				double GR = par.r*par.v;
				double GU = GR*par.u;
				double GV = GR*par.v + par.p;
				double GW = GR * par.w;
				double GE = (fRE + par.p)*par.v;

				double HR = par.r * par.w;
				double HU = HR * par.u;
				double HV = HR * par.v;
				double HW = HR * par.w + par.p;
				double HE = (fRE + par.p) * par.w;

				double dFdx = getDfDx(iBF, iCell, gp) * gw;
				double dFdy = getDfDy(iBF, iCell, gp) * gw;
				double dFdz = getDfDz(iBF, iCell, gp) * gw;

				sRO += (FR * dFdx + GR * dFdy + HR * dFdz);
				sRU += (FU * dFdx + GU * dFdy + HU * dFdz);
				sRV += (FV * dFdx + GV * dFdy + HV * dFdz);
				sRW += (FW * dFdx + GW * dFdy + HW * dFdz);
				sRE += (FE * dFdx + GE * dFdy + HE * dFdz);
			}
			sRO *= cellJ[iCell];
			sRU *= cellJ[iCell];
			sRV *= cellJ[iCell];
			sRW *= cellJ[iCell];
			sRE *= cellJ[iCell];
			
			int shift = 0;
			tmpArr[shift + iBF] = sRO; shift += BASE_FUNC_COUNT;
			tmpArr[shift + iBF] = sRU; shift += BASE_FUNC_COUNT;
			tmpArr[shift + iBF] = sRV; shift += BASE_FUNC_COUNT;
			tmpArr[shift + iBF] = sRW; shift += BASE_FUNC_COUNT;
			tmpArr[shift + iBF] = sRE; //shift += BASE_FUNC_COUNT;
		}

		solverMtx->addRightElement(iCell, tmpArr);
	}

	/* surf integral */

	memset(tmpCFL, 0, grid.cCount*sizeof(double));
	
	for (int iFace = 0; iFace < grid.fCount; iFace++) {
		memset(tmpArr1, 0, sizeof(double)*matrDim);
		memset(tmpArr2, 0, sizeof(double)*matrDim);

		int c1 = grid.faces[iFace].c1;
		int c2 = grid.faces[iFace].c2;
		if (c2 >= 0) {
			for (int iBF = 0; iBF < BASE_FUNC_COUNT; iBF++) {
				double sRO1 = 0.0;
				double sRU1 = 0.0;
				double sRV1 = 0.0;
				double sRW1 = 0.0;
				double sRE1 = 0.0;

				double sRO2 = 0.0;
				double sRU2 = 0.0;
				double sRV2 = 0.0;
				double sRW2 = 0.0;
				double sRE2 = 0.0;
				for (int iGP = 0; iGP < GP_FACE_COUNT; iGP++) {
					double fRO, fRU, fRV, fRW, fRE;
					double FR, FU, FV, FW, FE;

					Point& gp = faceGP[iFace][iGP];
					double gw = faceGW[iFace][iGP];

					getFields(fRO, fRU, fRV, fRW, fRE, c1, gp);
					Param par1;
					consToPar(fRO, fRU, fRV, fRW, fRE, par1);
					Material& mat1 = getMaterial(c1);
					mat1.URS(par1, 0); // p=p(r,e)

					getFields(fRO, fRU, fRV, fRW, fRE, c2, gp);
					Param par2;
					consToPar(fRO, fRU, fRV, fRW, fRE, par2);
					Material& mat2 = getMaterial(c2);
					mat2.URS(par2, 0); // p=p(r,e)

					calcFlux(FR, FU, FV, FW, FE, par1, par2, grid.faces[iFace].n, getGAM(c1));
					
					// вычисляем спектральный радиус для вычисления шага по времени
					if (STEADY) {
						double u1 = sqrt(par1.U2()) + par1.cz;
						double u2 = sqrt(par2.U2()) + par2.cz;
						double lambda = _max_(u1, u2);
						lambda *= gw;
						lambda *= faceJ[iFace];
						tmpCFL[c1] += lambda;
						tmpCFL[c2] += lambda;
					}

					double cGP1 = gw * getF(iBF, c1, gp);
					double cGP2 = gw * getF(iBF, c2, gp);

					sRO1 += FR * cGP1;
					sRU1 += FU * cGP1;
					sRV1 += FV * cGP1;
					sRW1 += FW * cGP1;
					sRE1 += FE * cGP1;

					sRO2 += FR * cGP2;
					sRU2 += FU * cGP2;
					sRV2 += FV * cGP2;
					sRW2 += FW * cGP2;
					sRE2 += FE * cGP2;
				}

				sRO1 *= faceJ[iFace];
				sRU1 *= faceJ[iFace];
				sRV1 *= faceJ[iFace];
				sRW1 *= faceJ[iFace];
				sRE1 *= faceJ[iFace];

				int shift = 0;
				tmpArr1[shift + iBF] = -sRO1; shift += BASE_FUNC_COUNT;
				tmpArr1[shift + iBF] = -sRU1; shift += BASE_FUNC_COUNT;
				tmpArr1[shift + iBF] = -sRV1; shift += BASE_FUNC_COUNT;
				tmpArr1[shift + iBF] = -sRW1; shift += BASE_FUNC_COUNT;
				tmpArr1[shift + iBF] = -sRE1; //shift += BASE_FUNC_COUNT;

				sRO2 *= faceJ[iFace];
				sRU2 *= faceJ[iFace];
				sRV2 *= faceJ[iFace];
				sRW2 *= faceJ[iFace];
				sRE2 *= faceJ[iFace];

				shift = 0;
				tmpArr2[shift + iBF] = sRO2; shift += BASE_FUNC_COUNT;
				tmpArr2[shift + iBF] = sRU2; shift += BASE_FUNC_COUNT;
				tmpArr2[shift + iBF] = sRV2; shift += BASE_FUNC_COUNT;
				tmpArr2[shift + iBF] = sRW2; shift += BASE_FUNC_COUNT;
				tmpArr2[shift + iBF] = sRE2; //shift += BASE_FUNC_COUNT;
			}
			
			solverMtx->addRightElement(c1, tmpArr1);
			solverMtx->addRightElement(c2, tmpArr2);
		}
		else {
			for (int iBF = 0; iBF < BASE_FUNC_COUNT; iBF++) {
				double sRO1 = 0.0;
				double sRU1 = 0.0;
				double sRV1 = 0.0;
				double sRW1 = 0.0;
				double sRE1 = 0.0;
				for (int iGP = 0; iGP < GP_FACE_COUNT; iGP++) {
					double fRO, fRU, fRV, fRW, fRE;
					double FR, FU, FV, FW, FE;

					Point& gp = faceGP[iFace][iGP];
					double gw = faceGW[iFace][iGP];

					getFields(fRO, fRU, fRV, fRW, fRE, c1, gp);
					Param par1;
					consToPar(fRO, fRU, fRV, fRW, fRE, par1);
					Material& mat = getMaterial(c1);
                    mat.URS(par1, 0); // p=p(r,e)
                    mat.URS(par1, 1);

					Param par2;
					boundaryCond(iFace, par1, par2);

					calcFlux(FR, FU, FV, FW, FE, par1, par2, grid.faces[iFace].n, getGAM(c1));

					// вычисляем спектральный радиус для вычисления шага по времени
					if (STEADY) {
						double u1 = sqrt(par1.U2()) + par1.cz;
						double u2 = sqrt(par2.U2()) + par2.cz;
						double lambda = _max_(u1, u2);
						lambda *= gw;
						lambda *= faceJ[iFace];
						tmpCFL[c1] += lambda;
					}

					double cGP1 = gw * getF(iBF, c1, gp);

					sRO1 += FR * cGP1;
					sRU1 += FU * cGP1;
					sRV1 += FV * cGP1;
					sRW1 += FW * cGP1;
					sRE1 += FE * cGP1;;

				}

				sRO1 *= faceJ[iFace];
				sRU1 *= faceJ[iFace];
				sRV1 *= faceJ[iFace];
				sRW1 *= faceJ[iFace];
				sRE1 *= faceJ[iFace];

				int shift = 0;
				tmpArr1[shift + iBF] = -sRO1; shift += BASE_FUNC_COUNT;
				tmpArr1[shift + iBF] = -sRU1; shift += BASE_FUNC_COUNT;
				tmpArr1[shift + iBF] = -sRV1; shift += BASE_FUNC_COUNT;
				tmpArr1[shift + iBF] = -sRW1; shift += BASE_FUNC_COUNT;
				tmpArr1[shift + iBF] = -sRE1; //shift += BASE_FUNC_COUNT;

			}

			solverMtx->addRightElement(c1, tmpArr1);
		}
	}
}


void FEM_DG_IMPLICIT::run()
{
	int solverErr = 0;
	double t = 0.0;
	int step = 0;
	long totalCalcTime = 0;
	while (t < TMAX && step < STEP_MAX) {
		long timeStart, timeEnd;
		timeStart = clock();

		if (!solverErr) step++;
		
		if (STEADY) {
			calcTimeStep();
		}
		else {
			t += TAU;
		}

		long tmStart = clock();
		solverErr = 0;
		solverMtx->zero();

		/* Заполняем правую часть */
		calcRHS();
        if (VISCOSITY) {
            /* Заполняем правую часть от диффузионных членов уравнений */
            calcViscousRHS();
        }

		/* Вычисляем шаги по времени в ячейках по насчитанным ранее значениям спектра */
		if (STEADY) {
			double minTau = DBL_MAX;
			for (int iCell = 0; iCell < grid.cCount; iCell++)
			{
				cTau[iCell] = CFL*grid.cells[iCell].V / (tmpCFL[iCell] + 1.0e-100);
				if (cTau[iCell] < minTau) minTau = cTau[iCell];
			}
		}

		/* Заполняем элементы матрицы */
		calcMatrWithTau();		// вычисляем матрицы перед производной по времени
		calcIntegral();			// вычисляем интеграл от(dF / dU)*deltaU*dFi / dx
		calcMatrFlux();			// Вычисляем потоковые величины

		if (VISCOSITY) {
            /* Заполняем элементы матрицы от диффузионных членов уравнений */
            calcMatrTensor();            //!< Вычисляем матрицы перед компонентами тензора вязких напряжений
            calcViscousIntegral();    //!< Вычисляем интеграл от (dH / dU)*dFi / dx и (dG / dU)*dFi / dx
            calcMatrViscousFlux();    //!< Вычисляем потоковые величины от диффузионных членов
        }

		/* Решаем СЛАУ */
		int maxIter = SOLVER_ITER;
		double eps = SOLVER_EPS;

		solverErr = solverMtx->solve(eps, maxIter);

		if (solverErr == MatrixSolver::RESULT_OK) {

			//if (SMOOTHING) smoothingDelta(solverMtx->x);

			for (int cellIndex = 0, ind = 0; cellIndex < grid.cCount; cellIndex++)
			{
				Cell &cell = grid.cells[cellIndex];

				//if (cellIsLim(cellIndex))	continue;
				for (int iFld = 0; iFld < fieldCount; iFld++) {
					for (int iF = 0; iF < BASE_FUNC_COUNT; iF++) {
						fields[iFld][cellIndex][iF] += solverMtx->x[ind++];
					}
				}
			}
			
			
			if (limiter != NULL) {
				limiter->run();
			}

			for (int cellIndex = 0, ind = 0; cellIndex < grid.cCount; cellIndex++)
			{
				Cell &cell = grid.cells[cellIndex];

				Param par;
				convertConsToPar(cellIndex, par);
				if (par.r > limitRmax)			{ par.r = limitRmax;	setCellFlagLim(cellIndex); }
				if (par.r < limitRmin)			{ par.r = limitRmin;	setCellFlagLim(cellIndex); }
				if (fabs(par.u) > limitUmax)		{ par.u = limitUmax*par.u / fabs(par.u);	setCellFlagLim(cellIndex); }
				if (fabs(par.v) > limitUmax)		{ par.v = limitUmax*par.v / fabs(par.v);	setCellFlagLim(cellIndex); }
				if (fabs(par.w) > limitUmax)		{ par.w = limitUmax * par.w / fabs(par.w);	setCellFlagLim(cellIndex); }
				if (par.p > limitPmax)			{ par.p = limitPmax;	setCellFlagLim(cellIndex); }
				if (par.p < limitPmin)			{ par.p = limitPmin;	setCellFlagLim(cellIndex); }
				if (cellIsLim(cellIndex)) 		{ par.e = par.p / ((getGAM(cellIndex) - 1)*par.r); convertParToCons(cellIndex, par); }

				double fRO, fRU, fRV, fRW, fRE;
				for (int iGP = 0; iGP < GP_CELL_COUNT; iGP++) {

					fRO = getField(FIELD_RO, cellIndex, cellGP[cellIndex][iGP]);
					fRU = getField(FIELD_RU, cellIndex, cellGP[cellIndex][iGP]);
					fRV = getField(FIELD_RV, cellIndex, cellGP[cellIndex][iGP]);
					fRW = getField(FIELD_RW, cellIndex, cellGP[cellIndex][iGP]);
					fRE = getField(FIELD_RE, cellIndex, cellGP[cellIndex][iGP]);
					consToPar(fRO, fRU, fRV, fRW, fRE, par);
					if (par.r > limitRmax)			{ par.r = limitRmax;	setCellFlagLim(cellIndex); }
					if (par.r < limitRmin)			{ par.r = limitRmin;	setCellFlagLim(cellIndex); }
					if (fabs(par.u) > limitUmax)		{ par.u = limitUmax*par.u / fabs(par.u);	setCellFlagLim(cellIndex); }
					if (fabs(par.v) > limitUmax)		{ par.v = limitUmax*par.v / fabs(par.v);	setCellFlagLim(cellIndex); }
					if (fabs(par.w) > limitUmax)		{ par.w = limitUmax * par.w / fabs(par.w);	setCellFlagLim(cellIndex); }
					if (par.p > limitPmax)			{ par.p = limitPmax;	setCellFlagLim(cellIndex); }
					if (par.p < limitPmin)			{ par.p = limitPmin;	setCellFlagLim(cellIndex); }
					if (cellIsLim(cellIndex)) 		{ par.e = par.p / ((getGAM(cellIndex) - 1)*par.r); convertParToCons(cellIndex, par); }
				}

				for (int iFace = 0; iFace < cell.fCount; iFace++) {
					int edgInd = cell.facesInd[iFace];
					for (int iGP = 0; iGP < GP_FACE_COUNT; iGP++) {
						fRO = getField(FIELD_RO, cellIndex, faceGP[edgInd][iGP]);
						fRU = getField(FIELD_RU, cellIndex, faceGP[edgInd][iGP]);
						fRV = getField(FIELD_RV, cellIndex, faceGP[edgInd][iGP]);
						fRW = getField(FIELD_RW, cellIndex, faceGP[edgInd][iGP]);
						fRE = getField(FIELD_RE, cellIndex, faceGP[edgInd][iGP]);
						consToPar(fRO, fRU, fRV, fRW, fRE, par);
						if (par.r > limitRmax)			{ par.r = limitRmax;	setCellFlagLim(cellIndex); }
						if (par.r < limitRmin)			{ par.r = limitRmin;	setCellFlagLim(cellIndex); }
						if (fabs(par.u) > limitUmax)		{ par.u = limitUmax*par.u / fabs(par.u);	setCellFlagLim(cellIndex); }
						if (fabs(par.v) > limitUmax)		{ par.v = limitUmax*par.v / fabs(par.v);	setCellFlagLim(cellIndex); }
						if (fabs(par.w) > limitUmax)		{ par.w = limitUmax*par.w / fabs(par.v);	setCellFlagLim(cellIndex); }
						if (par.p > limitPmax)			{ par.p = limitPmax;	setCellFlagLim(cellIndex); }
						if (par.p < limitPmin)			{ par.p = limitPmin;	setCellFlagLim(cellIndex); }
						if (cellIsLim(cellIndex)) 		{ par.e = par.p / ((getGAM(cellIndex) - 1)*par.r); convertParToCons(cellIndex, par); }
					}
				}


			}
			remediateLimCells();

			int limCells = getLimitedCellsCount();
			if (STEADY && (limCells >= maxLimCells)) decCFL();

			timeEnd = clock();
			totalCalcTime += (timeEnd - timeStart);
			if (step % FILE_SAVE_STEP == 0)
			{
				save(step);
			}
			if (step % PRINT_STEP == 0)
			{
				calcLiftForce();
				if (!STEADY) {

					log("step: %6d  time step: %.16f\tmax iter: %5d\tlim: %4d\tlift force (Fx, Fy, Fz) = (%.16f, %.16f, %.16f)\ttime: %6d ms\ttotal calc time: %ld\n", step, t, maxIter, limCells, Fx, Fy, Fz, timeEnd - timeStart, totalCalcTime);
				}
				else {
					log("step: %6d  max iter: %5d\tlim: %4d\tlift force (Fx, Fy, Fz) = (%.16f, %.16f, %.16f)\ttime: %6d ms\ttotal calc time: %ld\n", step, maxIter, limCells, Fx, Fy, Fz, timeEnd - timeStart, totalCalcTime);
				}
			}

			if (STEADY && (step % stepCFL == 0)) incCFL();
		}
		else {
			if (solverErr & MatrixSolver::RESULT_ERR_CONVERG) {
				log("Solver error: residual condition not considered.\n");
			}
			if (solverErr & MatrixSolver::RESULT_ERR_MAX_ITER) {
				log("Solver error: max iterations done.\n");
			}
			if (STEADY) {
				decCFL();
			}
			else {
				solverErr = 0;
			}
		}
	}
}

void calcFluxHLLCx(double& RI, double& EI, double& PI, double& UI, double& VI, double& WI, Param pL, Param pR, double GAM)
{
	//int             i;
	//size_t          c_count = charm_get_comp_count(p4est);
	double AGAM = (GAM - 1.0);
	double sl, sr, p_star, s_star, p_pvrs, ql, qr, tmp;

	p_pvrs = 0.5 * (pL.p + pR.p) - 0.5 * (pR.u - pL.u) * 0.25 * (pL.r + pR.r) * (pL.cz + pR.cz);
	p_star = (p_pvrs > 0.) ? p_pvrs : 0.;

	ql = (p_star <= pL.p) ? 1 : sqrt(1. + (pL.GAM + 1.) * (p_star / pL.p - 1.) / (2. * pL.GAM));
	qr = (p_star <= pR.p) ? 1 : sqrt(1. + (pR.GAM + 1.) * (p_star / pR.p - 1.) / (2. * pR.GAM));

	sl = pL.u - pL.cz * ql;
	sr = pR.u + pR.cz * qr;

	if (sl > sr) {
		tmp = sl;
		sl = sr;
		sr = tmp;
	}

	s_star = pR.p - pL.p;
	s_star += pL.r * pL.u * (sl - pL.u);
	s_star -= pR.r * pR.u * (sr - pR.u);
	s_star /= (pL.r * (sl - pL.u) - pR.r * (sr - pR.u));

	if (s_star < sl) s_star = sl;
	if (s_star > sr) s_star = sr;


	if (!((sl <= s_star) && (s_star <= sr))) {
		log("HLLC: inequaluty SL <= S* <= SR is FALSE.\n");
		exit(1);
	}

	if (sl >= 0.) {
		RI = pL.r * pL.u;
		UI = pL.r * pL.u * pL.u + pL.p;
		VI = pL.r * pL.v * pL.u;
		WI = pL.r * pL.w * pL.u;
		EI = (pL.r * pL.E + pL.p) * pL.u;
		/*for (int i = 0; i < c_count; i++) {
			qc[i] = pL.r * pL.u * pL.c[i];
		}*/
	}
	else if (sr <= 0.) {
		RI = pR.r * pR.u;
		UI = pR.r * pR.u * pR.u + pR.p;
		VI = pR.r * pR.v * pR.u;
		WI = pR.r * pR.w * pR.u;
		EI = (pR.r * pR.E + pR.p) * pR.u;
		/*for (i = 0; i < c_count; i++) {
			qc[i] = pR.r * pR.u * pR.c[i];
		}*/
	}
	else {
		if (s_star >= 0) {
			RI = F_HLLC_V( /*  UK, FK, SK, SS, PK, RK, VK */
				pL.r,
				pL.r * pL.u,
				sl, s_star, pL.p, pL.r, pL.u
			);
			UI = F_HLLC_U( /*  UK, FK, SK, SS, PK, RK, VK */
				pL.r * pL.u,
				pL.r * pL.u * pL.u + pL.p,
				sl, s_star, pL.p, pL.r, pL.u
			);
			VI = F_HLLC_V( /*  UK, FK, SK, SS, PK, RK, VK */
				pL.r * pL.v,
				pL.r * pL.u * pL.v,
				sl, s_star, pL.p, pL.r, pL.u
			);
			WI = F_HLLC_V( /*  UK, FK, SK, SS, PK, RK, VK */
				pL.r * pL.w,
				pL.r * pL.u * pL.w,
				sl, s_star, pL.p, pL.r, pL.u
			);
			EI = F_HLLC_E( /*  UK, FK, SK, SS, PK, RK, VK */
				pL.r * pL.E,
				(pL.r * pL.E + pL.p) * pL.u,
				sl, s_star, pL.p, pL.r, pL.u
			);
			//for (i = 0; i < c_count; i++) {
			//	qc[i] = F_HLLC_V( /*  UK, FK, SK, SS, PK, RK, VK */
			//		pL.r * pL.c[i],
			//		pL.r * pL.c[i] * pL.u,
			//		sl, s_star, pL.p, pL.r, pL.u
			//	);
			//}
		}
		else {
			RI = F_HLLC_V( /*  UK, FK, SK, SS, PK, RK, VK */
				pR.r,
				pR.r * pR.u,
				sr, s_star, pR.p, pR.r, pR.u
			);
			UI = F_HLLC_U( /*  UK, FK, SK, SS, PK, RK, VK */
				pR.r * pR.u,
				pR.r * pR.u * pR.u + pR.p,
				sr, s_star, pR.p, pR.r, pR.u
			);
			VI = F_HLLC_V( /*  UK, FK, SK, SS, PK, RK, VK */
				pR.r * pR.v,
				pR.r * pR.u * pR.v,
				sr, s_star, pR.p, pR.r, pR.u
			);
			WI = F_HLLC_V( /*  UK, FK, SK, SS, PK, RK, VK */
				pR.r * pR.w,
				pR.r * pR.u * pR.w,
				sr, s_star, pR.p, pR.r, pR.u
			);
			EI = F_HLLC_E( /*  UK, FK, SK, SS, PK, RK, VK */
				pR.r * pR.E,
				(pR.r * pR.E + pR.p) * pR.u,
				sr, s_star, pR.p, pR.r, pR.u
			);
			//for (i = 0; i < c_count; i++) {
			//	qc[i] = F_HLLC_V( /*  UK, FK, SK, SS, PK, RK, VK */
			//		pR.r * pR.c[i],
			//		pR.r * pR.c[i] * pR.u,
			//		sr, s_star, pR.p, pR.r, pR.u
			//	);
			//}
		}
	}
	PI = AGAM * EI * RI;
}

void FEM_DG_IMPLICIT::calcFlux(double& fr, double& u, double& v, double& w, double& fe, Param pL, Param pR, Vector n, double GAM)
{
	if (FLUX == FLUX_GODUNOV) {	// GODUNOV FLUX
		double RI, EI, PI, UI, VI, WI, UN, UTy, UTz;

		Vector nty, ntz;

		double ri = sqrt(n.x * n.x + n.y * n.y);
		if (ri > DBL_EPSILON) {
			nty.x = -n.y / ri;
			nty.y = n.x / ri;
			nty.z = 0.;
		}
		else {
			ri = sqrt(n.y * n.y + n.z * n.z);
			if (ri > DBL_EPSILON) {
				nty.x = 0.;
				nty.y = -n.z / ri;
				nty.z = n.y / ri;
			}
			else {
				ri = sqrt(n.x * n.x + n.z * n.z);
				if (ri > DBL_EPSILON) {
					nty.x = -n.z / ri;
					nty.y = 0.;
					nty.z = n.x / ri;
				}
				else {
					log("Error in defining tangential vectors");
				}
			}
		}
		
		ntz = vector_prod(n, nty);

		double unl = pL.u * n.x + pL.v * n.y + pL.w * n.z;
		double unr = pR.u * n.x + pR.v * n.y + pR.w * n.z;
		double utly = pL.u * nty.x + pL.v * nty.y + pL.w * nty.z;
		double utlz = pL.u * ntz.x + pL.v * ntz.y + pL.w * ntz.z;
		double utry = pR.u * nty.x + pR.v * nty.y + pR.w * nty.z;
		double utrz = pR.u * ntz.x + pR.v * ntz.y + pR.w * ntz.z;

		rim_orig(RI, EI, PI, UN, UTy, UTz, pL.r, pL.p, unl, utly, utlz, pR.r, pR.p, unr, utry, utrz, GAM);

		UI = UN * n.x + UTy * nty.x + UTz * ntz.x;
		VI = UN * n.y + UTy * nty.y + UTz * ntz.y;
		WI = UN * n.z + UTy * nty.z + UTz * ntz.z;

		fr = RI * UN;
		u = fr * UI + PI * n.x;
		v = fr * VI + PI * n.y;;
		w = fr * WI + PI * n.z;
		fe = (RI * (EI + 0.5 * (UI * UI + VI * VI + WI * WI)) + PI) * UN;

		return;
	}
	if (FLUX == FLUX_LAX) {	// LAX-FRIEDRIX FLUX
		double unl = pL.u * n.x + pL.v * n.y + pL.w * n.z;
		double unr = pR.u * n.x + pR.v * n.y + pR.w * n.z;
		double rol, rul, rvl, rwl, rel, ror, rur, rvr, rwr, rer;
		double alpha = _max_(fabs(unl) + sqrt(GAM*pL.p / pL.r), fabs(unr) + sqrt(GAM*pR.p / pR.r));
		rol = pL.r;
		rul = pL.r * pL.u;
		rvl = pL.r * pL.v;
		rwl = pL.r * pL.w;
		rel = pL.p / (GAM - 1.0) + 0.5 * pL.r * (pL.u * pL.u + pL.v * pL.v + pL.w * pL.w);
		ror = pR.r;
		rur = pR.r * pR.u;
		rvr = pR.r * pR.v;
		rwr = pR.r * pR.w;
		rer = pR.p / (GAM - 1.0) + 0.5 * pR.r * (pR.u * pR.u + pR.v * pR.v + pR.w * pR.w);
		double frl = rol*unl;
		double frr = ror*unr;
		fr = 0.5*(frr + frl - alpha*(ror - rol));
		u = 0.5 * (frr * pR.u + frl * pL.u + (pR.p + pL.p) * n.x - alpha * (rur - rul));
		v = 0.5 * (frr * pR.v + frl * pL.v + (pR.p + pL.p) * n.y - alpha * (rvr - rvl));
		w = 0.5 * (frr * pR.w + frl * pL.w + (pR.p + pL.p) * n.z - alpha * (rwr - rwl));
		fe = 0.5 * ((rer + pR.p) * unr + (rel + pL.p) * unl - alpha * (rer - rel));
		
		return;
	}
	if (FLUX == FLUX_CD) {	// FLUX CD
		double unl = pL.u * n.x + pL.v * n.y + pL.w * n.z;
		double unr = pR.u * n.x + pR.v * n.y + pR.w * n.z;
		double rol, rul, rvl, rwl, rel, ror, rur, rvr, rwr, rer;
		rol = pL.r;
		rul = pL.r * pL.u;
		rvl = pL.r * pL.v;
		rwl = pL.r * pL.w;
		rel = pL.p / (GAM - 1.0) + 0.5 * pL.r * (pL.u * pL.u + pL.v * pL.v + pL.w * pL.w);
		ror = pR.r;
		rur = pR.r * pR.u;
		rvr = pR.r * pR.v;
		rwr = pR.r * pR.w;
		rer = pR.p / (GAM - 1.0) + 0.5 * pR.r * (pR.u * pR.u + pR.v * pR.v + pR.w * pR.w);
		double frl = rol * unl;
		double frr = ror * unr;
		fr = 0.5 * (frr + frl);
		u = 0.5 * (frr * pR.u + frl * pL.u + (pR.p + pL.p) * n.x);
		v = 0.5 * (frr * pR.v + frl * pL.v + (pR.p + pL.p) * n.y);
		w = 0.5 * (frr * pR.w + frl * pL.w + (pR.p + pL.p) * n.z);
		fe = 0.5 * ((rer + pR.p) * unr + (rel + pL.p) * unl);

		return;
	}
	if (FLUX == FLUX_KIR) { // FLUX KIR
		static Param pB;
		static VECTOR QFL(5), QFR(5);
		static VECTOR UL(5), UR(5);

		double RI, EI, PI, UI, VI, WI, UN, UTy, UTz;
		Vector nty, ntz;

		//double ri = sqrt(n.x * n.x + n.y * n.y);
		//if (ri > DBL_EPSILON) {
		//	nty.x = -n.y / ri;
		//	nty.y = n.x / ri;
		//	nty.z = 0.;
		//}
		//else {
		//	ri = sqrt(n.y * n.y + n.z * n.z);
		//	if (ri > DBL_EPSILON) {
		//		nty.x = 0.;
		//		nty.y = -n.z / ri;
		//		nty.z = n.y / ri;
		//	}
		//	else {
		//		ri = sqrt(n.x * n.x + n.z * n.z);
		//		if (ri > DBL_EPSILON) {
		//			nty.x = -n.z / ri;
		//			nty.y = 0.;
		//			nty.z = n.x / ri;
		//		}
		//		else {
		//			log("Error in defining tangential vectors");
		//		}
		//	}
		//}
		//ntz = vector_prod(n, nty);

		//double unl = pL.u * n.x + pL.v * n.y + pL.w * n.z;
		//double unr = pR.u * n.x + pR.v * n.y + pR.w * n.z;
		//double utly = pL.u * nty.x + pL.v * nty.y + pL.w * nty.z;
		//double utlz = pL.u * ntz.x + pL.v * ntz.y + pL.w * ntz.z;
		//double utry = pR.u * nty.x + pR.v * nty.y + pR.w * nty.z;
		//double utrz = pR.u * ntz.x + pR.v * ntz.y + pR.w * ntz.z;

		//pL.u = unl;
		//pL.v = utly;
		//pL.w = utlz;

		//pR.u = unr;
		//pR.v = utry;
		//pR.w = utrz;

		//n.x = 1.; n.y = 0.; n.z = 0.;

		pL.un = pL.u * n.x + pL.v * n.y + pL.w * n.z;
		pR.un = pR.u * n.x + pR.v * n.y + pR.w * n.z;

		pB = pL;

		{{ // Схема ROE
			pB.GAM = pL.GAM;

			double gamma = pL.GAM;

			double SL = sqrt(pL.r);
			double SR = sqrt(pR.r);
			double S_ = 1.0 / (SL + SR);

			pB.r = SL * SR;

			pB.u = (SL * pL.u + SR * pR.u) * S_;
			pB.v = (SL * pL.v + SR * pR.v) * S_;
			pB.w = (SL * pL.w + SR * pR.w) * S_;

			pB.e = (SL * pL.e + SR * pR.e) * S_;
			pB.T = (SL * pL.T + SR * pR.T) * S_;

			pL.H = pL.e + 0.5*pL.U2() + pL.p / pL.r;
			pR.H = pR.e + 0.5*pR.U2() + pR.p / pR.r;

			pB.H = (SL * pL.H + SR * pR.H) * S_;

			pB.p = (pB.H - 0.5*pB.U2()) * pB.r * (gamma - 1.0) / gamma;

			pB.cz = sqrt(gamma * pB.p / pB.r);

			pB.un = pB.u * n.x + pB.v * n.y + pB.w * n.z;
		}}

		double ur = scalar_prod(pL.omegaR, n);

		double** mtrxA = allocMtx5();
		double H = pB.H;
		calcAbsA(mtrxA, pB.cz, pB.GAM, pB.u, n.x, pB.v, n.y, pB.w, n.z, H);

		UL[0] = 1.0;
		UL[1] = pL.u;
		UL[2] = pL.v;
		UL[3] = pL.w;
		UL[4] = pL.e + 0.5 * pL.U2();
		UL *= pL.r;

		UR[0] = 1.0;
		UR[1] = pR.u;
		UR[2] = pR.v;
		UR[3] = pR.w;
		UR[4] = pR.e + 0.5 * pR.U2();
		UR *= pR.r;

		QFL = UL;
		QFL *= (pL.un - ur);
		QFL[1] += pL.p * n.x;
		QFL[2] += pL.p * n.y;
		QFL[3] += pL.p * n.z;
		QFL[4] += pL.p * pL.un;

		QFR = UR;
		QFR *= (pR.un - ur);
		QFR[1] += pR.p * n.x;
		QFR[2] += pR.p * n.y;
		QFR[3] += pR.p * n.z;
		QFR[4] += pR.p * pR.un;

		static VECTOR  QF(5);

		//double faL_Dis = 1.0;

		QF = UL;
		QF -= UR;
		/*QF = mA * QF;*/
		QF *= mtrxA;
		//QF *= faL_Dis;

		QF += QFR;
		QF += QFL;
		QF *= 0.5;

		/*UI = QF[1] * n.x + QF[2] * nty.x + QF[3] * ntz.x;
		VI = QF[1] * n.y + QF[2] * nty.y + QF[3] * ntz.y;
		WI = QF[1] * n.z + QF[2] * nty.z + QF[3] * ntz.z;*/

		/*RI = QF[0];
		UN = QF[1];
		UTy = QF[2];
		UTz = QF[3];
		EI = QF[4];
		PI = (1. - GAM) * EI * RI;

		UI = QF[1] * n.x + QF[2] * nty.x + QF[3] * ntz.x;
		VI = QF[1] * n.y + QF[2] * nty.y + QF[3] * ntz.y;
		WI = QF[1] * n.z + QF[2] * nty.z + QF[3] * ntz.z;

		fr = RI;
		u = UI;
		v = VI;
		w = WI;
		fe = EI;*/

		fr = QF[0];
		u = QF[1];
		v = QF[2];
		w = QF[3];
		fe = QF[4];

		freeMtx5(mtrxA);

		return;
	}
	if (FLUX == FLUX_AUSMP) { // FLUX AUSMP
		static Param pB;
		static VECTOR UL(5), UR(5);

		double Kp = 0.25;
		double Ku = 0.75;
		double Sigma = 0.4;

		pL.un = pL.u * n.x + pL.v * n.y + pL.w * n.z;
		pR.un = pR.u * n.x + pR.v * n.y + pR.w * n.z;

		double LUN = pL.un - scalar_prod(pL.omegaR, n);
		double RUN = pR.un - scalar_prod(pR.omegaR, n);

		pL.H = pL.e + 0.5*pL.U2() + pL.p / pL.r;
		pR.H = pR.e + 0.5*pR.U2() + pR.p / pR.r;

		double A_L = ::sqrt(pL.GAM * pL.p / pL.r);
		double A_R = ::sqrt(pL.GAM * pR.p / pR.r);

		double A_avg = 0.5 * (A_L + A_R);

		//	double fVmag_L = pL.Q();
		//	double fVmag_R = pR.Q();

		double mag_L = LUN;
		double mag_R = RUN;


		//	REAL fMp2 = 0.5 * ( _pow_2_( fVmag_L ) + _pow_2_( fVmag_R ) ) / _pow_2_( fA_avg );
		//	REAL fM2  = _min_( 1.0, _max_( fMp2, _pow_2_( fM_co ) ) );

		pL.ur = A_avg;
		double M2 = _min_(1.0, POW_2(pL.ur / A_avg));

		double M2_SQR = ::sqrt(M2);

		double AlphaMo = M2_SQR * (2.0 - M2_SQR);

		//fAlphaMo = 1.0;

		double M_L = LUN / A_L;
		double M_R = RUN / A_R;


		//	REAL fM_1L = 0.5 * ( fM_L + ::fabs( fM_L ) ); 
		//	REAL fM_1R = 0.5 * ( fM_R - ::fabs( fM_R ) ); 

		double M_1L = _max_(M_L, 0.0);
		double M_1R = _min_(M_R, 0.0);

		double M_2L = 0.25 * POW_2(M_L + 1.0);
		double M_2R = -0.25 * POW_2(M_R - 1.0);

		double alfa = 3.0 / 16.0 * (-4.0 + 5.0 * POW_2(AlphaMo));
		double beta_ = 1.0 / 8.0;


		double M_plus, M_minus, P_plus, P_minus;


		if (POW_2(M_L) < 1.0)
		{
			M_plus = M_2L * (1.0 - 16.0 * beta_ * M_2R);
			P_plus = M_2L * ((2.0 - M_L) - 16.0 * alfa * M_L * M_2R);
		}
		else {
			M_plus = M_1L; P_plus = M_1L / M_L;
		} // if

		if (POW_2(M_R) < 1.0)
		{
			M_minus = M_2R * (1.0 + 16.0 * beta_ * M_2L);
			P_minus = M_2R * ((-2.0 - M_R) + 16.0 * alfa * M_R * M_2L);
		}
		else {
			M_minus = M_1R; P_minus = M_1R / M_R;
		} // if


	//	if ( _pow_2_( fM_1L ) + _pow_2_( fM_1R ) < 1.0e-30 )
		if (M_1L == 0.0 && M_1R == 0.0)
		{
			M_plus = M_2L * (1.0 - 16.0 * beta_ * M_2R);
			M_minus = M_2R * (1.0 + 16.0 * beta_ * M_2L);
			P_plus = M_2L * ((2.0 - M_L) - 16.0 * alfa * M_L * M_2R);
			P_minus = M_2R * ((-2.0 - M_R) + 16.0 * alfa * M_R * M_2L);
		} // if

		//fP_plus = 1.0 - fP_minus;

		double Sum = P_plus + P_minus;

		P_plus = P_plus / Sum;
		P_minus = P_minus / Sum;

		double Ro_avg = 0.5 * (pL.r + pR.r);

		//double fM_k = fM_plus + fM_minus;
		double M_k = M_plus + M_minus - Kp / AlphaMo * _max_(1.0 - Sigma * M2, 0.0) * (pR.p - pL.p) / (Ro_avg * POW_2(A_avg));

		double P_k = P_plus * pL.p + P_minus * pR.p - Ku * P_plus * P_minus * (pL.r + pR.r) * AlphaMo * A_avg * (RUN - LUN);


		//  REAL fP_k = 0.5 * ( rP.P + rE.P);
		//fi
		UL[0] = 1.0;
		UL[1] = pL.u;
		UL[2] = pL.v;
		UL[3] = pL.w;
		UL[4] = pL.H;

		//vUP   *= rP.R;

		UR[0] = 1.0;
		UR[1] = pR.u;
		UR[2] = pR.v;
		UR[3] = pR.w;
		UR[4] = pR.H;

		//vUE   *= rE.R;

		pB.omegaR = pL.omegaR;
		double ur = scalar_prod(pL.omegaR, n);

		pB.un = M_k * A_avg;


		if (M_k >= 0.0) {

			//pB = pL;
			UL *= pL.r;

			//g_vQF = rB.UN * vUP;
			fr = pB.un * UL[0];
			u = pB.un * UL[1];
			v = pB.un * UL[2];
			w = pB.un * UL[3];
			fe = pB.un * UL[4];
			LUN = pL.r * pB.un;

		}
		else {

			//pB = pR;
			UR *= pR.r;

			fr = pB.un * UR[0];
			u = pB.un * UR[1];
			v = pB.un * UR[2];
			w = pB.un * UR[3];
			fe = pB.un * UR[4];
			RUN = pR.r * pB.un;
		} // if

		u += P_k * n.x;
		v += P_k * n.y;
		w += P_k * n.z;
		fe += P_k * scalar_prod(pL.omegaR, n);

		/*g_rB.UN = rB.UN;
		g_rB.P = fP_k;
		g_rB.CZ = fabs(rB.UN) + fA_avg;
		g_fP = fP_k;*/

		return;
	}
	if (FLUX == FLUX_AUSMPW) { // FLUX AUSMPW
		static Param pB;
		static VECTOR UL(5), UR(5);
		static VECTOR QF(5);

		double Kp = 0.25;
		double Ku = 0.75;
		double Sigma = 0.4;

		pL.un = pL.u * n.x + pL.v * n.y + pL.w * n.z;
		pR.un = pR.u * n.x + pR.v * n.y + pR.w * n.z;

		double LUN = pL.un - scalar_prod(pL.omegaR, n);
		double RUN = pR.un - scalar_prod(pR.omegaR, n);

		pL.H = pL.e + 0.5*pL.U2() + pL.p / pL.r;
		pR.H = pR.e + 0.5*pR.U2() + pR.p / pR.r;

		// double A_L = ::sqrt(pL.GAM * pL.p / pL.r);
		// double A_R = ::sqrt(pL.GAM * pR.p / pR.r);
		
		double A_L = pL.cz;
		double A_R = pR.cz;

		double A_avg = 0.5 * (A_L + A_R);

		//	double fVmag_L = pL.Q();
		//	double fVmag_R = pR.Q();

		double mag_L = LUN;
		double mag_R = RUN;


		//	REAL fMp2 = 0.5 * ( _pow_2_( fVmag_L ) + _pow_2_( fVmag_R ) ) / _pow_2_( fA_avg );
		//	REAL fM2  = _min_( 1.0, _max_( fMp2, _pow_2_( fM_co ) ) );

		pL.ur = A_avg;
		double M2 = _min_(1.0, POW_2(pL.ur / A_avg));

		double M2_SQR = ::sqrt(M2);

		double AlphaMo = M2_SQR * (2.0 - M2_SQR);

		//fAlphaMo = 1.0;

		double M_L = LUN / A_L;
		double M_R = RUN / A_R;


		double M_1L = 0.5 * ( M_L + ::fabs( M_L ) ); 
		double M_1R = 0.5 * ( M_R - ::fabs( M_R ) ); 

		// double M_1L = _max_(M_L, 0.0);
		// double M_1R = _min_(M_R, 0.0);

		double M_2L = 0.25 * POW_2(M_L + 1.0);
		double M_2R = -0.25 * POW_2(M_R - 1.0);

		double alfa = 3.0 / 16.0 * (-4.0 + 5.0 * POW_2(AlphaMo));
		double beta = 1.0 / 8.0;


		double M_plus_b0, M_minus_b0, M_plus_b18, M_minus_b18, P_plus, P_minus;


		if (POW_2(M_L) <= 1.0)
		{
			M_plus_b18 = M_2L + beta * POW_2( M_L * M_L - 1.0 ); 
			M_plus_b0 = M_2L ; 
			P_plus = M_2L * ( 2.0 - M_L ) + 3.0 / 16.0 * M_L * POW_2( M_L * M_L - 1.0 );
		}
		else {
			M_plus_b18 = M_1L; 
			//fP_plus = fM_1L / fM_L;
			P_plus = 0.5 * ( 1.0 + _sign_( M_L ) );
		} // if

		if (POW_2(M_R) <= 1.0)
		{
			M_minus_b18 = M_2R - beta * POW_2( M_R * M_R - 1.0 );
			M_minus_b0 = M_2R;

			P_minus = 0.25 * POW_2( M_R - 1.0 ) * ( 2.0 + M_R ) - 3.0 / 16.0 * M_R * POW_2( M_R * M_R - 1.0 ); 
		}
		else {
			M_minus_b18 = M_1R;
			//fP_minus = fM_1R / fM_R;
			P_minus = 0.5 * ( 1.0 - _sign_( M_R ) );
		} // if


		if ( P_plus == 0.0 && P_minus == 0.0 )
		{
			P_plus  = 0.5 * pL.p;
			P_minus = 0.5 * pL.e;
		} // if

		double Sum = P_plus + P_minus;

		P_plus = P_plus / Sum;
		P_minus = P_minus / Sum;

		//Нахожу касательные скорости
		Vector V_L( pL.u, pL.v, pL.w );
		Vector V_R( pR.u, pR.v, pR.w );

		Vector Vn_L = n; Vn_L *= pL.un; Vector Vt_L = V_L - Vn_L;
		Vector Vn_R = n; Vn_R *= pR.un; Vector Vt_R = V_R - Vn_R;


		double OMEGA = 1.0 - _pow_3_( _min_( pL.p / pR.p, pR.p / pL.p ) );

		double PL_xy_L = 0.0;
		double PL_xy_R = 0.0;
		
		double MinXY_L = _min_( pL.p / pR.p, pR.p / pL.p );
		double MinXY_R = _min_( pR.p / pL.p, pL.p / pR.p );

		if ( MinXY_L >= 0.75 || MinXY_L < 1.0 )  PL_xy_L = 4.0 * MinXY_L - 3.0;
		if ( MinXY_R >= 0.75 || MinXY_R < 1.0 )  PL_xy_R = 4.0 * MinXY_R - 3.0;

		// Судя по статье закомментированное слагаемое призвано вносить искусственную вязкость
		// по моим расчетам ее величина порядка 1-2% лобового сопротивления на задаче NACA0012
		// пока комментирую 
		double P_k = P_plus * pL.p + P_minus * pR.p;	// - fKu * fP_plus * fP_minus * ( rP.R + rE.R ) * fAlphaMo * fA_avg * ( frEUN - frPUN );

		double F_L = 0.0;
		double F_R = 0.0;

		if ( _pow_2_( M_L ) < 1.0 ) F_L = ( pL.p / P_k - 1.0 ) * PL_xy_L * M_plus_b0  * _min_( 1.0, ::pow( _mag_(Vn_L) / A_avg , 0.25 ) );
		if ( _pow_2_( M_R ) < 1.0 ) F_R = ( pR.p / P_k - 1.0 ) * PL_xy_R * M_minus_b0 * _min_( 1.0, ::pow( _mag_(Vn_R) / A_avg , 0.25 ) );
		

		double M_k = M_plus_b18 + M_minus_b18;


		//  REAL fP_k = 0.5 * ( rP.P + rE.P);
		//fi
		UL[0] = 1.0;
		UL[1] = pL.u;
		UL[2] = pL.v;
		UL[3] = pL.w;
		UL[4] = pL.H;
		UL   *= pL.r;

		//vUP   *= rP.R;

		UR[0] = 1.0;
		UR[1] = pR.u;
		UR[2] = pR.v;
		UR[3] = pR.w;
		UR[4] = pR.H;
		UR   *= pR.r;
		//vUE   *= rE.R;

		double M_plus, M_minus;
	
		pB.un = M_k * A_avg;
		
		double Ro_avg = 0.5 * ( pL.r + pR.r );

		double Stab = - Kp / AlphaMo * _max_ ( 1.0 - Sigma * M2, 0.0 ) * (pR.p - pL.p) / ( Ro_avg * _pow_2_( A_avg ) );

		if ( M_k >=0.0 )
		{
			M_plus  = M_plus_b18 + M_minus_b18 - ( M_minus_b18 * OMEGA * ( 1.0 + F_R ) ) + ( F_L * M_plus_b18 + F_R * M_minus_b18) + Stab;
			M_minus = M_minus_b18 * OMEGA * ( 1.0 + F_R ) + Stab;
		
			//g_fRUN = rP.R * rB.UN * fDS;
			//g_rB = rP;
		}
		else
		{
			M_plus  = M_plus_b18 * OMEGA * ( 1.0 + F_L ) + Stab;
			M_minus = M_plus_b18 + M_minus_b18 - ( M_plus_b18 * OMEGA * ( 1.0 + F_L ) ) + ( F_L * M_plus_b18 + F_R * M_minus_b18 ) + Stab;

			//g_fRUN = rE.R * rB.UN * fDS;
			//g_rB = rE;
		} // if

		UL *= M_plus * A_avg;
		UR *= M_minus * A_avg;
		
		QF = UL;
		QF += UR;

		QF[ 1 ] += P_k * n.x;
		QF[ 2 ] += P_k * n.y;
		QF[ 3 ] += P_k * n.z;

		QF[ 4 ] += P_k * scalar_prod( pL.omegaR, n );

		fr = QF[0];
		u = QF[1];
		v = QF[2];
		w = QF[3];
		fe = QF[4];

		/*g_rB.UN = rB.UN;
		g_rB.P = fP_k;
		g_rB.CZ = fabs(rB.UN) + fA_avg;
		g_fP = fP_k;*/

		return;
	}
	if (FLUX == FLUX_HLLC) { // FLUX HLLC
		//int i, j;
		//charm_real_t ri, ei, pi, uu[3], uv[3];
		//charm_real_t nt[3][3], vv[2][3], vn[2][3];
		//charm_real_t r_[2], u_[2], v_[2], w_[2], p_[2];
		//size_t          c_count = charm_get_comp_count(p4est);
		//charm_real_t _qu, _qv, _qw;
		double RI, EI, PI, UI, VI, WI, UN, UTy, UTz;
		Vector nty, ntz;

		double ri = sqrt(n.x * n.x + n.y * n.y);
		if (ri > DBL_EPSILON) {
			nty.x = -n.y / ri;
			nty.y = n.x / ri;
			nty.z = 0.;
		}
		else {
			ri = sqrt(n.y * n.y + n.z * n.z);
			if (ri > DBL_EPSILON) {
				nty.x = 0.;
				nty.y = -n.z / ri;
				nty.z = n.y / ri;
			}
			else {
				ri = sqrt(n.x * n.x + n.z * n.z);
				if (ri > DBL_EPSILON) {
					nty.x = -n.z / ri;
					nty.y = 0.;
					nty.z = n.x / ri;
				}
				else {
					log("Error in defining tangential vectors");
				}
			}
		}
		ntz = vector_prod(n, nty);

		double unl = pL.u * n.x + pL.v * n.y + pL.w * n.z;
		double unr = pR.u * n.x + pR.v * n.y + pR.w * n.z;
		double utly = pL.u * nty.x + pL.v * nty.y + pL.w * nty.z;
		double utlz = pL.u * ntz.x + pL.v * ntz.y + pL.w * ntz.z;
		double utry = pR.u * nty.x + pR.v * nty.y + pR.w * nty.z;
		double utrz = pR.u * ntz.x + pR.v * ntz.y + pR.w * ntz.z;

		pL.u = unl;
		pL.v = utly;
		pL.w = utlz;

		pR.u = unr;
		pR.v = utry;
		pR.w = utrz;

		calcFluxHLLCx(RI, EI, PI, UN, UTy, UTz, pL, pR, GAM);

		UI = UN * n.x + UTy * nty.x + UTz * ntz.x;
		VI = UN * n.y + UTy * nty.y + UTz * ntz.y;
		WI = UN * n.z + UTy * nty.z + UTz * ntz.z;

		/*fr = RI * UN;
		u = fr * UI + PI * n.x;
		v = fr * VI + PI * n.y;;
		w = fr * WI + PI * n.z;
		fe = (RI * (EI + 0.5 * (UI * UI + VI * VI + WI * WI)) + PI) * UN;*/
		fr = RI;
		u = UI;
		v = VI;
		w = WI;
		fe = EI;

		
		return;
	}
}

void FEM_DG_IMPLICIT::calcFluxFields(Param& pF, Param pL, Param pR, Vector n, double GAM)
{
	if (FLUX == FLUX_GODUNOV) {	// GODUNOV FLUX
		double RI, EI, PI, UI, VI, WI;

		rim_orig(RI, EI, PI, UI, VI, WI, pL.r, pL.p, pL.u, pL.v, pL.w, pR.r, pR.p, pR.u, pR.v, pR.w, GAM);

		pF.r = RI;
		pF.u = UI;
		pF.v = VI;
		pF.w = WI;
		pF.e = EI;
		pF.p = PI;
		pF.cz = sqrt(GAM * pF.p / pF.r);
		pF.E = pF.r * (pF.e + 0.5 * (pF.u * pF.u + pF.v * pF.v + pF.w * pF.w));

		return;
	}
	else {
		calcRoeAverage(pF, pL, pR, GAM, n);
	}
	return;
}

void FEM_DG_IMPLICIT::boundaryCond(int iFace, Param& pL, Param& pR)
{
    Face &face = grid.faces[iFace];
    int c1 = face.c1;
    Material& m = getMaterial(c1);
    if (face.bnd) {
        face.bnd->run(iFace, pL, pR);
        m.URS(pR, 2);
        m.URS(pR, 1);
        pR.E = pR.e + 0.5*(pR.U2());
        return;
    }
    else {
        char msg[128];
        sprintf(msg, "Not defined boundary condition for face %d\n", iFace);
        throw Exception(msg, Exception::TYPE_BOUND_UNKNOWN);
    }
}

void FEM_DG_IMPLICIT::addSmallMatrToBigMatr(double **mB, double **mS, int i, int j)
{

	int ii = i * BASE_FUNC_COUNT;
	int jj = j * BASE_FUNC_COUNT;
	for (int i1 = 0; i1 < BASE_FUNC_COUNT; i1++) {
		for (int j1 = 0; j1 < BASE_FUNC_COUNT; j1++) {
			mB[ii + i1][jj + j1] += mS[i1][j1];
		}
	}
}

void FEM_DG_IMPLICIT::defineOrientationOnFaces()
{
	for (int iFace = 0; iFace < grid.fCount; iFace++) {
		beta += grid.faces[iFace].n;
		beta /= grid.fCount;
	}
	beta.norm();
	for (int iFace = 0; iFace < grid.fCount; iFace++) {
		Face& face = grid.faces[iFace];
		Vector n = face.n;

		if (scalar_prod(n, beta) > 0) {
			face.L = face.c1;
			face.R = face.c2;
		}
		else {
			face.L = face.c2;
			face.R = face.c1;
		}
	}
}

void FEM_DG_IMPLICIT::done()
{
	// TODO
}

void FEM_DG_IMPLICIT::calcLiftForce()
{
	//const double width = 1.0; // предполагаемая ширина профиля по z.
	Param		 par;
	Fx = Fy = Fz = 0.0;
	for (int iFace = 0; iFace < grid.fCount; ++iFace)
	{
		if (!grid.faces[iFace].bnd) continue;
		if (grid.faces[iFace].bnd->faceType == CFDBoundary::TYPE_ID_WALL_NO_SLIP || grid.faces[iFace].bnd->faceType == CFDBoundary::TYPE_ID_WALL_SLIP)
		{

			int			cellIndex = grid.faces[iFace].c1;
			double		nx = grid.faces[iFace].n.x;
			double		ny = grid.faces[iFace].n.y;
			double		nz = grid.faces[iFace].n.z;
			if (cellIndex < 0) {
				cellIndex = grid.faces[iFace].c2;
				nx *= -1.0;
				ny *= -1.0;
				nz *= -1.0;
			}
			convertConsToPar(cellIndex, par);
			Fx += (par.p * P_ - PRef) * grid.faces[iFace].S * nx;
			Fy += (par.p * P_ - PRef) * grid.faces[iFace].S * ny;
			Fz += (par.p * P_ - PRef) * grid.faces[iFace].S * nz;

		}
	}
}

double **FEM_DG_IMPLICIT::allocMtx11() {
    double		**tempMtx11 = new double*[11];
    for (int i = 0; i < 11; ++i) tempMtx11[i] = new double[11];
    return tempMtx11;
}

void FEM_DG_IMPLICIT::freeMtx11(double **mtx11) {
    for (int i = 0; i < 11; ++i)
        delete[] mtx11[i];
    delete[] mtx11;
}

void FEM_DG_IMPLICIT::multMtx11(double **dst11, double **srcA11, double **srcB11) {
    double sum;
    for (int i = 0; i < 11; ++i)
    {
        for (int j = 0; j < 11; ++j)
        {
            sum = 0;
            for (int k = 0; k < 11; ++k)
                sum += srcA11[i][k] * srcB11[k][j];
            dst11[i][j] = sum;
        }
    }
}

void FEM_DG_IMPLICIT::clearMtx11(double **mtx11) {
    for (int i = 0; i < 11; ++i)
        for (int j = 0; j < 11; ++j)
            mtx11[i][j] = 0;
}

void FEM_DG_IMPLICIT::getTensorComponents(double& fTAU_XX, double& fTAU_YY, double& fTAU_ZZ, double& fTAU_XY, double& fTAU_XZ, double& fTAU_YZ, int iCell, Point p) {
    getTensorComponents(fTAU_XX, fTAU_YY, fTAU_ZZ, fTAU_XY, fTAU_XZ, fTAU_YZ, iCell, p.x, p.y, p.z);
}

void
FEM_DG_IMPLICIT::getTensorComponents(double& fTAU_XX, double& fTAU_YY, double& fTAU_ZZ, double& fTAU_XY, double& fTAU_XZ, double& fTAU_YZ, int iCell, double x, double y, double z) {
    fTAU_XX = getField(FIELD_TAU_XX, iCell, x, y, z);
    fTAU_YY = getField(FIELD_TAU_YY, iCell, x, y, z);
    fTAU_ZZ = getField(FIELD_TAU_ZZ, iCell, x, y, z);
	fTAU_XY = getField(FIELD_TAU_XY, iCell, x, y, z);
	fTAU_XZ = getField(FIELD_TAU_XZ, iCell, x, y, z);
	fTAU_YZ = getField(FIELD_TAU_YZ, iCell, x, y, z);
}

void FEM_DG_IMPLICIT::calcMatrTensor() {
    for (int iCell = 0; iCell < grid.cCount; iCell++) {

        fillMtx(matrBig, 0.0, MATR_DIM);


        for (int i = 0; i < BASE_FUNC_COUNT; i++) {
            for (int j = 0; j < BASE_FUNC_COUNT; j++) {
                matrSmall[i][j] = matrA[iCell][i][j];
            }
        }

        for (int ii = FIELD_COUNT; ii < FIELD_COUNT_EXT; ii++) {
            addSmallMatrToBigMatr(matrBig, matrSmall, ii, ii);
        }

        solverMtx->addMatrElement(iCell, iCell, matrBig);

    }
}

void FEM_DG_IMPLICIT::calcViscousIntegral() {
    double fRO, fRU, fRV, fRW, fRE, fTAU_XX, fTAU_YY, fTAU_ZZ, fTAU_XY, fTAU_XZ, fTAU_YZ;
    double FR, u, v, w, FE;
    double **mx11 = allocMtx11();
    double **my11 = allocMtx11();
	double** mz11 = allocMtx11();

    for (int iCell = 0; iCell < grid.cCount; iCell++) {

        fillMtx(matrBig, 0.0, MATR_DIM);

        for (int iGP = 0; iGP < GP_CELL_COUNT; iGP++) {
            Point& gp = cellGP[iCell][iGP];
            double gw = cellGW[iCell][iGP];
            getFields(fRO, fRU, fRV, fRW, fRE, iCell, gp);
            Param par;
            consToPar(fRO, fRU, fRV, fRW, fRE, par);
            getTensorComponents(fTAU_XX, fTAU_YY, fTAU_ZZ, fTAU_XY, fTAU_XZ, fTAU_YZ, iCell, gp);
            Material& mat = getMaterial(iCell);
            mat.URS(par, 0); // p=p(r,e)
            mat.getML(par);
            double H = par.E + par.p / par.r;

            calcJ(mx11, par.r, par.u, 1.0, par.v, 0.0, par.w, 0.0, par.ML, fTAU_XX, fTAU_YY, fTAU_ZZ, fTAU_XY, fTAU_XZ, fTAU_YZ);
            calcJ(my11, par.r, par.u, 0.0, par.v, 1.0, par.w, 0.0, par.ML, fTAU_XX, fTAU_YY, fTAU_ZZ, fTAU_XY, fTAU_XZ, fTAU_YZ);
			calcJ(mz11, par.r, par.u, 0.0, par.v, 0.0, par.w, 1.0, par.ML, fTAU_XX, fTAU_YY, fTAU_ZZ, fTAU_XY, fTAU_XZ, fTAU_YZ);
            for (int i = 0; i < FIELD_COUNT_EXT; i++) {
                for (int j = 0; j < FIELD_COUNT_EXT; j++) {
                    for (int ii = 0; ii < BASE_FUNC_COUNT; ii++) {
                        for (int jj = 0; jj < BASE_FUNC_COUNT; jj++) {
                            matrSmall[ii][jj] = mx11[i][j] * getDfDx(ii, iCell, gp);
                            matrSmall[ii][jj] += my11[i][j] * getDfDy(ii, iCell, gp);
							matrSmall[ii][jj] += mz11[i][j] * getDfDz(ii, iCell, gp);
                            matrSmall[ii][jj] *= getF(jj, iCell, gp); // basis unction
                            matrSmall[ii][jj] *= gw;
                        }
                    }
                    addSmallMatrToBigMatr(matrBig, matrSmall, i, j);
                }
            }

        }

        multMtxToVal(matrBig, cellJ[iCell] * SIGMA, MATR_DIM);

        solverMtx->addMatrElement(iCell, iCell, matrBig);
    }

    freeMtx11(mx11);
    freeMtx11(my11);
	freeMtx11(mz11);
}

void FEM_DG_IMPLICIT::calcMatrViscousFlux() {
    double  **Amtx11P, **Amtx11M;

    Amtx11P = allocMtx11();
    Amtx11M = allocMtx11();


	double fRO1, fRU1, fRV1, fRW1, fRE1, fRO2, fRU2, fRV2, fRW2, fRE2;
	double fTAU_XX1, fTAU_YY1, fTAU_ZZ1, fTAU_XY1, fTAU_XZ1, fTAU_YZ1, fTAU_XX2, fTAU_YY2, fTAU_ZZ2, fTAU_XY2, fTAU_XZ2, fTAU_YZ2;
    Param par1, par2;

    for (int iFace = 0; iFace < grid.fCount; iFace++) {
        Face& face = grid.faces[iFace];
        Vector&	n = grid.faces[iFace].n;
        int c1 = face.c1;
        int c2 = face.c2;

        if (c2 >= 0) {

            fillMtx(matrBig, 0.0, MATR_DIM);
            fillMtx(matrBig2, 0.0, MATR_DIM);

            for (int iGP = 0; iGP < GP_FACE_COUNT; iGP++) {
                Point& gp = faceGP[iFace][iGP];
                double gw = faceGW[iFace][iGP];
                getFields(fRO1, fRU1, fRV1, fRW1, fRE1, c1, gp);
                consToPar(fRO1, fRU1, fRV1, fRW1, fRE1, par1);
                Material& mat1 = getMaterial(c1);
                mat1.URS(par1, 0); // p=p(r,e)
                mat1.getML(par1);
                getTensorComponents(fTAU_XX1, fTAU_YY1, fTAU_ZZ1, fTAU_XY1, fTAU_XZ1, fTAU_YZ1, c1, gp);

                getFields(fRO2, fRU2, fRV2, fRW2, fRE2, c2, gp);
                consToPar(fRO2, fRU2, fRV2, fRW2, fRE2, par2);
                Material& mat2 = getMaterial(c2);
                mat2.URS(par2, 0); // p=p(r,e)
                mat2.getML(par2);
                getTensorComponents(fTAU_XX2, fTAU_YY2, fTAU_ZZ2, fTAU_XY2, fTAU_XZ2, fTAU_YZ2, c2, gp);

                calcJ(Amtx11P, par1.r, par1.u, n.x, par1.v, n.y, par1.w, n.z, par1.ML, fTAU_XX1, fTAU_YY1, fTAU_ZZ1, fTAU_XY1, fTAU_XZ1, fTAU_YZ1);
                calcJ(Amtx11M, par2.r, par2.u, n.x, par2.v, n.y, par2.w, n.z, par2.ML, fTAU_XX2, fTAU_YY2, fTAU_ZZ2, fTAU_XY2, fTAU_XZ2, fTAU_YZ2);
                for (int i = 0; i < FIELD_COUNT_EXT; i++) {
                    for (int j = 0; j < FIELD_COUNT_EXT; j++) {
                        for (int ii = 0; ii < BASE_FUNC_COUNT; ii++) {
                            for (int jj = 0; jj < BASE_FUNC_COUNT; jj++) {
								matrSmall[ii][jj] = 0.5 * Amtx11P[i][j] * getF(jj, c1, gp) * getF(ii, c1, gp) * gw;
								matrSmall2[ii][jj] = 0.5 * Amtx11M[i][j] * getF(jj, c2, gp) * getF(ii, c1, gp) * gw;
                            }
                        }
                        addSmallMatrToBigMatr(matrBig,  matrSmall, i, j);
                        addSmallMatrToBigMatr(matrBig2, matrSmall2, i, j);
                    }
                }

            }

            multMtxToVal(matrBig, faceJ[iFace] * SIGMA, MATR_DIM);
            multMtxToVal(matrBig2, faceJ[iFace] * SIGMA, MATR_DIM);

            solverMtx->addMatrElement(c1, c1, matrBig);
            solverMtx->addMatrElement(c1, c2, matrBig2);



            fillMtx(matrBig, 0.0, MATR_DIM);
            fillMtx(matrBig2, 0.0, MATR_DIM);

            for (int iGP = 0; iGP < GP_FACE_COUNT; iGP++) {
                Point& gp = faceGP[iFace][iGP];
                double gw = faceGW[iFace][iGP];
                getFields(fRO1, fRU1, fRV1, fRW1, fRE1, c1, gp);
                consToPar(fRO1, fRU1, fRV1, fRW1, fRE1, par1);
                Material& mat1 = getMaterial(c1);
                mat1.URS(par1, 0); // p=p(r,e)
                mat1.getML(par1);
                getTensorComponents(fTAU_XX1, fTAU_YY1, fTAU_ZZ1, fTAU_XY1, fTAU_XZ1, fTAU_YZ1, c1, gp);

                getFields(fRO2, fRU2, fRV2, fRW2, fRE2, c2, gp);
                consToPar(fRO2, fRU2, fRV2, fRW2, fRE2, par2);
                Material& mat2 = getMaterial(c2);
                mat2.URS(par2, 0); // p=p(r,e)
                mat2.getML(par2);
                getTensorComponents(fTAU_XX2, fTAU_YY2, fTAU_ZZ2, fTAU_XY2, fTAU_XZ2, fTAU_YZ2, c2, gp);

                calcJ(Amtx11P, par1.r, par1.u, -n.x, par1.v, -n.y, par1.w, -n.z, par1.ML, fTAU_XX1, fTAU_YY1, fTAU_ZZ1, fTAU_XY1, fTAU_XZ1, fTAU_YZ1);
                calcJ(Amtx11M, par2.r, par2.u, -n.x, par2.v, -n.y, par2.w, -n.z, par2.ML, fTAU_XX2, fTAU_YY2, fTAU_ZZ2, fTAU_XY2, fTAU_XZ2, fTAU_YZ2);
                for (int i = 0; i < FIELD_COUNT_EXT; i++) {
                    for (int j = 0; j < FIELD_COUNT_EXT; j++) {
                        for (int ii = 0; ii < BASE_FUNC_COUNT; ii++) {
                            for (int jj = 0; jj < BASE_FUNC_COUNT; jj++) {
								matrSmall[ii][jj] = 0.5 * Amtx11P[i][j] * getF(jj, c2, gp) * getF(ii, c2, gp) * gw;
								matrSmall2[ii][jj] = 0.5 * Amtx11M[i][j] * getF(jj, c1, gp) * getF(ii, c2, gp) * gw;
                            }
                        }
                        addSmallMatrToBigMatr(matrBig, matrSmall, i, j);
                        addSmallMatrToBigMatr(matrBig2, matrSmall2, i, j);
                    }
                }

            }

            multMtxToVal(matrBig, faceJ[iFace] * SIGMA, MATR_DIM);
            multMtxToVal(matrBig2, faceJ[iFace] * SIGMA, MATR_DIM);

            solverMtx->addMatrElement(c2, c2, matrBig);
            solverMtx->addMatrElement(c2, c1, matrBig2);

        }
    }

    freeMtx11(Amtx11P);
    freeMtx11(Amtx11M);
}

void FEM_DG_IMPLICIT::calcJ(double **dst11, double r, double u, double nx, double v, double ny, double w, double nz, double mu, double tau_xx, double tau_yy, double tau_zz, double tau_xy, double tau_xz, double tau_yz) {
    
	dst11[0][0] = 0.0;   
	dst11[0][1] = 0.0;   
	dst11[0][2] = 0.0;   
	dst11[0][3] = 0.0;   
	dst11[0][4] = 0.0;   
	dst11[0][5] = 0.0;   
	dst11[0][6] = 0.0;
	dst11[0][7] = 0.0;
	dst11[0][8] = 0.0;
	dst11[0][9] = 0.0;
	dst11[0][10] = 0.0;

	dst11[1][0] = 0.0;
	dst11[1][1] = 0.0;
	dst11[1][2] = 0.0;
	dst11[1][3] = 0.0;
	dst11[1][4] = 0.0;
	dst11[1][5] = nx;
	dst11[1][6] = 0.0;
	dst11[1][7] = 0.0;
	dst11[1][8] = ny;
	dst11[1][9] = nz;
	dst11[1][10] = 0.0;

	dst11[2][0] = 0.0;
	dst11[2][1] = 0.0;
	dst11[2][2] = 0.0;
	dst11[2][3] = 0.0;
	dst11[2][4] = 0.0;
	dst11[2][5] = 0.0;
	dst11[2][6] = ny;
	dst11[2][7] = 0.0;
	dst11[2][8] = nx;
	dst11[2][9] = 0.0;
	dst11[2][10] = nz;

	dst11[3][0] = 0.0;
	dst11[3][1] = 0.0;
	dst11[3][2] = 0.0;
	dst11[3][3] = 0.0;
	dst11[3][4] = 0.0;
	dst11[3][5] = 0.0;
	dst11[3][6] = 0.0;
	dst11[3][7] = nz;
	dst11[3][8] = 0.0;
	dst11[3][9] = nx;
	dst11[3][10] = ny;

	dst11[4][0] = -(nx * (u * tau_xx + v * tau_xy + w * tau_xz) + ny * (u * tau_xy + v * tau_yy + w * tau_yz) + nz * (u * tau_xz + v * tau_yz + w * tau_zz)) / r;
	dst11[4][1] = (tau_xx * nx + tau_xy * ny + tau_xz * nz) / r;
	dst11[4][2] = (tau_xy * nx + tau_yy * ny + tau_yz * nz) / r;
	dst11[4][3] = (tau_xz * nx + tau_yz * ny + tau_zz * nz) / r;;
	dst11[4][4] = 0.0;
	dst11[4][5] = u * nx;
	dst11[4][6] = v * ny;
	dst11[4][7] = w * nz;
	dst11[4][8] = v * nx + u * ny;
	dst11[4][9] = w * nx + u * nz;
	dst11[4][10] = w * ny + v * nz;

	dst11[5][0] = -mu * (4. * u * nx - 2. * v * ny - 2. * w * nz) / 3. / r;
	dst11[5][1] = 4. * mu * nx / 3. / r;
	dst11[5][2] = -2. * mu * ny / 3. / r;
	dst11[5][3] = -2. * mu * nz / 3. / r;
	dst11[5][4] = 0.0;
	dst11[5][5] = 0.0;
	dst11[5][6] = 0.0;
	dst11[5][7] = 0.0;
	dst11[5][8] = 0.0;
	dst11[5][9] = 0.0;
	dst11[5][10] = 0.0;

	dst11[6][0] = -mu * (-2. * u * nx + 4. * v * ny - 2. * w * nz) / 3. / r;
	dst11[6][1] = -2. * mu * nx / 3. / r;
	dst11[6][2] = 4. * mu * ny / 3. / r;
	dst11[6][3] = -2. * mu * nz / 3. / r;
	dst11[6][4] = 0.0;
	dst11[6][5] = 0.0;
	dst11[6][6] = 0.0;
	dst11[6][7] = 0.0;
	dst11[6][8] = 0.0;
	dst11[6][9] = 0.0;
	dst11[6][10] = 0.0;

	dst11[7][0] = -mu * (-2. * u * nx - 2. * v * ny + 4. * w * nz) / 3. / r;
	dst11[7][1] = -2. * mu * nx / 3. / r;
	dst11[7][2] = -2. * mu * ny / 3. / r;
	dst11[7][3] = 4. * mu * nz / 3. / r;
	dst11[7][4] = 0.0;
	dst11[7][5] = 0.0;
	dst11[7][6] = 0.0;
	dst11[7][7] = 0.0;
	dst11[7][8] = 0.0;
	dst11[7][9] = 0.0;
	dst11[7][10] = 0.0;

	dst11[8][0] = -mu * (v * nx + u * ny) / r;
	dst11[8][1] = mu * ny / r;
	dst11[8][2] = mu * nx / r;
	dst11[8][3] = 0.0;
	dst11[8][4] = 0.0;
	dst11[8][5] = 0.0;
	dst11[8][6] = 0.0;
	dst11[8][7] = 0.0;
	dst11[8][8] = 0.0;
	dst11[8][9] = 0.0;
	dst11[8][10] = 0.0;

	dst11[9][0] = -mu * (w * nx + u * nz) / r;
	dst11[9][1] = mu * nz / r;
	dst11[9][2] = 0.0;
	dst11[9][3] = mu * nx / r;
	dst11[9][4] = 0.0;
	dst11[9][5] = 0.0;
	dst11[9][6] = 0.0;
	dst11[9][7] = 0.0;
	dst11[9][8] = 0.0;
	dst11[9][9] = 0.0;
	dst11[9][10] = 0.0;

	dst11[10][0] = -mu * (w * ny + v * nz) / r;
	dst11[10][1] = 0.0;
	dst11[10][2] = mu * nz / r;
	dst11[10][3] = mu * ny / r;
	dst11[10][4] = 0.0;
	dst11[10][5] = 0.0;
	dst11[10][6] = 0.0;
	dst11[10][7] = 0.0;
	dst11[10][8] = 0.0;
	dst11[10][9] = 0.0;
	dst11[10][10] = 0.0;
}

void FEM_DG_IMPLICIT::calcViscousRHS() {

    /* volume integral */

    for (int iCell = 0; iCell < grid.cCount; iCell++) {
        memset(tmpArr, 0, sizeof(double)*MATR_DIM);
        for (int iBF = 0; iBF < BASE_FUNC_COUNT; iBF++) {
            double s1 = 0.0;
            double s2 = 0.0;
            double s3 = 0.0;
            double s4 = 0.0;
            double s5 = 0.0;
            double s6 = 0.0;
            double s7 = 0.0;
			double s8 = 0.0;
			double s9 = 0.0;
			double s10 = 0.0;
			double s11 = 0.0;
            for (int iGP = 0; iGP < GP_CELL_COUNT; iGP++) {
                Point& gp = cellGP[iCell][iGP];
                double gw = cellGW[iCell][iGP];
                double fRO, fRU, fRV, fRW, fRE, fTAU_XX, fTAU_YY, fTAU_ZZ, fTAU_XY, fTAU_XZ, fTAU_YZ;
                getFields(fRO, fRU, fRV, fRW, fRE, iCell, gp);
                Param par;
                consToPar(fRO, fRU, fRV, fRW, fRE, par);
                getTensorComponents(fTAU_XX, fTAU_YY, fTAU_ZZ, fTAU_XY, fTAU_XZ, fTAU_YZ, iCell, gp);
                Material& mat = getMaterial(iCell);
                mat.URS(par, 0); // p=p(r,e)
                mat.getML(par);

                double F1 = 0.;
                double F2 = fTAU_XX;
                double F3 = fTAU_XY;
				double F4 = fTAU_XZ;
				double F5 = par.u * fTAU_XX + par.v * fTAU_XY + par.w * fTAU_XZ;
                double F6 = par.ML*(4.*par.u / 3.);
                double F7 = par.ML*(-2.*par.u / 3.);
                double F8 = par.ML*(-2.*par.u / 3.);
				double F9 = par.ML * (par.v);
				double F10 = par.ML * (par.w);
				double F11 = 0.0;

                double G1 = 0.0;
                double G2 = fTAU_XY;
                double G3 = fTAU_YY;
				double G4 = fTAU_YZ;
				double G5 = par.u * fTAU_XY + par.v * fTAU_YY + par.w * fTAU_YZ;
                double G6 = par.ML*(-2.*par.v / 3.);
                double G7 = par.ML*(4.*par.v / 3.);
                double G8 = par.ML*(-2.*par.v / 3.);
				double G9 = par.ML * (par.u);
				double G10 = 0.0;
				double G11 = par.ML * (par.v);

				double H1 = 0.0;
				double H2 = fTAU_XZ;
				double H3 = fTAU_YZ;
				double H4 = fTAU_ZZ;
				double H5 = par.u * fTAU_XZ + par.v * fTAU_YZ + par.w * fTAU_ZZ;
				double H6 = par.ML * (-2. * par.w / 3.);
				double H7 = par.ML * (-2. * par.w / 3.);
				double H8 = par.ML * (4. * par.w / 3.);
				double H9 = 0.0;
				double H10 = par.ML * (par.u);
				double H11 = par.ML * (par.v);

                double dFdx = getDfDx(iBF, iCell, gp) * gw;
                double dFdy = getDfDy(iBF, iCell, gp) * gw;
				double dFdz = getDfDy(iBF, iCell, gp) * gw;

				s1 += (F1 * dFdx + G1 * dFdy + H1 * dFdz);
				s2 += (F2 * dFdx + G2 * dFdy + H2 * dFdz);
				s3 += (F3 * dFdx + G3 * dFdy + H3 * dFdz);
				s4 += (F4 * dFdx + G4 * dFdy + H4 * dFdz);
				s5 += (F5 * dFdx + G5 * dFdy + H5 * dFdz);
				s6 += (F6 * dFdx + G6 * dFdy + H6 * dFdz);
				s7 += (F7 * dFdx + G7 * dFdy + H7 * dFdz);
				s8 += (F8 * dFdx + G8 * dFdy + H8 * dFdz);
				s9 += (F9 * dFdx + G9 * dFdy + H9 * dFdz);
				s10 += (F10 * dFdx + G10 * dFdy + H10 * dFdz);
				s11 += (F11 * dFdx + G11 * dFdy + H11 * dFdz);
            }
            s1 *= cellJ[iCell];
            s2 *= cellJ[iCell];
            s3 *= cellJ[iCell];
            s4 *= cellJ[iCell];
            s5 *= cellJ[iCell];
            s6 *= cellJ[iCell];
            s7 *= cellJ[iCell];
			s8 *= cellJ[iCell];
			s9 *= cellJ[iCell];
			s10 *= cellJ[iCell];
			s11 *= cellJ[iCell];

            int shift = 0;
            tmpArr[shift + iBF] = -s1; shift += BASE_FUNC_COUNT;
            tmpArr[shift + iBF] = -s2; shift += BASE_FUNC_COUNT;
            tmpArr[shift + iBF] = -s3; shift += BASE_FUNC_COUNT;
            tmpArr[shift + iBF] = -s4; shift += BASE_FUNC_COUNT;
            tmpArr[shift + iBF] = -s5; shift += BASE_FUNC_COUNT;
            tmpArr[shift + iBF] = -s6; shift += BASE_FUNC_COUNT;
            tmpArr[shift + iBF] = -s7; shift += BASE_FUNC_COUNT;
			tmpArr[shift + iBF] = -s8; shift += BASE_FUNC_COUNT;
			tmpArr[shift + iBF] = -s9; shift += BASE_FUNC_COUNT;
			tmpArr[shift + iBF] = -s10; shift += BASE_FUNC_COUNT;
			tmpArr[shift + iBF] = -s11; //shift += BASE_FUNC_COUNT;
        }

        solverMtx->addRightElement(iCell, tmpArr);
    }
    // W^m contribution
    for (int iCell = 0; iCell < grid.cCount; iCell++) {
        memset(tmpArr, 0, sizeof(double)*MATR_DIM);
        for (int iBF = 0; iBF < BASE_FUNC_COUNT; iBF++) {
            double s6 = 0.0;
            double s7 = 0.0;
            double s8 = 0.0;
			double s9 = 0.0;
			double s10 = 0.0;
			double s11 = 0.0;
            for (int iGP = 0; iGP < GP_CELL_COUNT; iGP++) {
                Point& gp = cellGP[iCell][iGP];
                double gw = cellGW[iCell][iGP];
                double fTAU_XX, fTAU_YY, fTAU_ZZ, fTAU_XY, fTAU_XZ, fTAU_YZ;
                getTensorComponents(fTAU_XX, fTAU_YY, fTAU_ZZ, fTAU_XY, fTAU_XZ, fTAU_YZ, iCell, gp);

                double cGP = gw*getF(iBF, iCell, gp);

				s6 += fTAU_XX * cGP;
				s7 += fTAU_YY * cGP;
				s8 += fTAU_ZZ * cGP;
				s9 += fTAU_XY * cGP;
				s10 += fTAU_XZ * cGP;
				s11 += fTAU_YZ * cGP;
            }
            s6 *= cellJ[iCell];
            s7 *= cellJ[iCell];
            s8 *= cellJ[iCell];
			s9 *= cellJ[iCell];
			s10 *= cellJ[iCell];
			s11 *= cellJ[iCell];

            int shift = FIELD_COUNT*BASE_FUNC_COUNT;
            tmpArr[shift + iBF] = -s6; shift += BASE_FUNC_COUNT;
            tmpArr[shift + iBF] = -s7; shift += BASE_FUNC_COUNT;
            tmpArr[shift + iBF] = -s8; shift += BASE_FUNC_COUNT;
			tmpArr[shift + iBF] = -s9; shift += BASE_FUNC_COUNT;
			tmpArr[shift + iBF] = -s10; shift += BASE_FUNC_COUNT;
			tmpArr[shift + iBF] = -s11; //shift += BASE_FUNC_COUNT;
        }

        solverMtx->addRightElement(iCell, tmpArr);
    }

    /* surf integral */

    //memset(tmpCFL, 0, grid.cCount*sizeof(double));

    for (int iFace = 0; iFace < grid.fCount; iFace++) {
        memset(tmpArr1, 0, sizeof(double)*MATR_DIM);
        memset(tmpArr2, 0, sizeof(double)*MATR_DIM);

        int c1 = grid.faces[iFace].c1;
        int c2 = grid.faces[iFace].c2;
        Vector n = grid.faces[iFace].n;
        if (c2 >= 0) {
            for (int iBF = 0; iBF < BASE_FUNC_COUNT; iBF++) {
                //double s11 = 0.0;
                double s21 = 0.0;
                double s31 = 0.0;
                double s41 = 0.0;
                double s51 = 0.0;
                double s61 = 0.0;
                double s71 = 0.0;
				double s81 = 0.0;
				double s91 = 0.0;
				double s101 = 0.0;
				double s111 = 0.0;

                //double s12 = 0.0;
                double s22 = 0.0;
                double s32 = 0.0;
                double s42 = 0.0;
                double s52 = 0.0;
                double s62 = 0.0;
                double s72 = 0.0;
				double s82 = 0.0;
				double s92 = 0.0;
				double s102 = 0.0;
				double s112 = 0.0;

                for (int iGP = 0; iGP < GP_FACE_COUNT; iGP++) {
                    double fRO, fRU, fRV, fRW, fRE, fTAU_XX1, fTAU_YY1, fTAU_ZZ1, fTAU_XY1, fTAU_XZ1, fTAU_YZ1, fTAU_XX2, fTAU_YY2, fTAU_ZZ2, fTAU_XY2, fTAU_XZ2, fTAU_YZ2;
                    double FS2, FS3, FS4, FS5, FS6, FS7, FS8, FS9, FS10, FS11;

                    Point& gp = faceGP[iFace][iGP];
                    double gw = faceGW[iFace][iGP];

                    getFields(fRO, fRU, fRV, fRW, fRE, c1, gp);
                    Param par1;
                    consToPar(fRO, fRU, fRV, fRW, fRE, par1);
                    Material& mat1 = getMaterial(c1);
                    mat1.URS(par1, 0); // p=p(r,e)
                    mat1.getML(par1);
                    getTensorComponents(fTAU_XX1, fTAU_YY1, fTAU_ZZ1, fTAU_XY1, fTAU_XZ1, fTAU_YZ1, c1, gp);

                    getFields(fRO, fRU, fRV, fRW, fRE, c2, gp);
                    Param par2;
                    consToPar(fRO, fRU, fRV, fRW, fRE, par2);
                    Material& mat2 = getMaterial(c2);
                    mat2.URS(par2, 0); // p=p(r,e)
                    mat2.getML(par2);
                    getTensorComponents(fTAU_XX2, fTAU_YY2, fTAU_ZZ2, fTAU_XY2, fTAU_XZ2, fTAU_YZ2, c2, gp);

					double FS2L = fTAU_XX1 * n.x + fTAU_XY1 * n.y + fTAU_XZ1 * n.z;
					double FS3L = fTAU_XY1 * n.x + fTAU_YY1 * n.y + fTAU_YZ1 * n.z;
					double FS4L = fTAU_XZ1 * n.x + fTAU_YZ1 * n.y + fTAU_ZZ1 * n.z;
					double FS5L = (par1.u * fTAU_XX1 + par1.v * fTAU_XY1 + par1.w*fTAU_XZ1) * n.x + (par1.u * fTAU_XY1 + par1.v * fTAU_YY1 + par1.w * fTAU_YZ1) * n.y + (par1.u * fTAU_XZ1 + par1.v * fTAU_YZ1 + par1.w * fTAU_ZZ1) * n.z;
					double FS6L = par1.ML * (4. * par1.u * n.x - 2. * par1.v * n.y - 2. * par1.w * n.z) / 3.;
					double FS7L = par1.ML * (-2. * par1.u * n.x + 4. * par1.v * n.y - 2. * par1.w * n.z) / 3.;
					double FS8L = par1.ML * (-2. * par1.u * n.x - 2. * par1.v * n.y + 4. * par1.w * n.z) / 3.;
					double FS9L = par1.ML * (par1.v * n.x + par1.u * n.y);
					double FS10L = par1.ML * (par1.w * n.x + par1.u * n.z);
					double FS11L = par1.ML * (par1.w * n.y + par1.v * n.z);

					double FS2R = fTAU_XX2 * n.x + fTAU_XY2 * n.y + fTAU_XZ2 * n.z;
					double FS3R = fTAU_XY2 * n.x + fTAU_YY2 * n.y + fTAU_YZ2 * n.z;
					double FS4R = fTAU_XZ2 * n.x + fTAU_YZ2 * n.y + fTAU_ZZ2 * n.z;
					double FS5R = (par2.u * fTAU_XX2 + par2.v * fTAU_XY2 + par2.w * fTAU_XZ2) * n.x + (par2.u * fTAU_XY2 + par2.v * fTAU_YY2 + par2.w * fTAU_YZ2) * n.y + (par2.u * fTAU_XZ2 + par2.v * fTAU_YZ2 + par2.w * fTAU_ZZ2) * n.z;
					double FS6R = par2.ML * (4. * par2.u * n.x - 2. * par2.v * n.y - 2. * par2.w * n.z) / 3.;
					double FS7R = par2.ML * (-2. * par2.u * n.x + 4. * par2.v * n.y - 2. * par2.w * n.z) / 3.;
					double FS8R = par2.ML * (-2. * par2.u * n.x - 2. * par2.v * n.y + 4. * par2.w * n.z) / 3.;
					double FS9R = par2.ML * (par2.v * n.x + par2.u * n.y);
					double FS10R = par2.ML * (par2.w * n.x + par2.u * n.z);
					double FS11R = par2.ML * (par2.w * n.y + par2.v * n.z);

					FS2 = 0.5 * (FS2L + FS2R);
					FS3 = 0.5 * (FS3L + FS3R);
					FS4 = 0.5 * (FS4L + FS4R);
					FS5 = 0.5 * (FS5L + FS5R);
					FS6 = 0.5 * (FS6L + FS6R);
					FS7 = 0.5 * (FS7L + FS7R);
					FS8 = 0.5 * (FS8L + FS8R);
					FS9 = 0.5 * (FS9L + FS9R);
					FS10 = 0.5 * (FS10L + FS10R);
					FS11 = 0.5 * (FS11L + FS11R);

                    double cGP1 = gw * getF(iBF, c1, gp);
                    double cGP2 = gw * getF(iBF, c2, gp);

                    //s11 += 0.;
					s21 += FS2 * cGP1;
					s31 += FS3 * cGP1;
					s41 += FS4 * cGP1;
					s51 += FS5 * cGP1;
					s61 += FS6 * cGP1;
					s71 += FS7 * cGP1;
					s81 += FS8 * cGP1;
					s91 += FS9 * cGP1;
					s101 += FS10 * cGP1;
					s111 += FS11 * cGP1;

                    //s12 += 0.;
					s22 += FS2 * cGP2;
					s32 += FS3 * cGP2;
					s42 += FS4 * cGP2;
					s52 += FS5 * cGP2;
					s62 += FS6 * cGP2;
					s72 += FS7 * cGP2;
					s82 += FS8 * cGP2;
					s92 += FS9 * cGP2;
					s102 += FS10 * cGP2;
					s112 += FS11 * cGP2;
                }

                //s11 *= faceJ[iFace];
                s21 *= faceJ[iFace];
                s31 *= faceJ[iFace];
                s41 *= faceJ[iFace];
                s51 *= faceJ[iFace];
                s61 *= faceJ[iFace];
                s71 *= faceJ[iFace];
				s81 *= faceJ[iFace];
				s91 *= faceJ[iFace];
				s101 *= faceJ[iFace];
				s111 *= faceJ[iFace];

                int shift = BASE_FUNC_COUNT;
                tmpArr1[shift + iBF] = s21; shift += BASE_FUNC_COUNT;
                tmpArr1[shift + iBF] = s31; shift += BASE_FUNC_COUNT;
                tmpArr1[shift + iBF] = s41; shift += BASE_FUNC_COUNT;
                tmpArr1[shift + iBF] = s51; shift += BASE_FUNC_COUNT;
                tmpArr1[shift + iBF] = s61; shift += BASE_FUNC_COUNT;
                tmpArr1[shift + iBF] = s71; shift += BASE_FUNC_COUNT;
				tmpArr1[shift + iBF] = s81; shift += BASE_FUNC_COUNT;
				tmpArr1[shift + iBF] = s91; shift += BASE_FUNC_COUNT;
				tmpArr1[shift + iBF] = s101; shift += BASE_FUNC_COUNT;
				tmpArr1[shift + iBF] = s111; //shift += BASE_FUNC_COUNT;

                s22 *= faceJ[iFace];
                s32 *= faceJ[iFace];
                s42 *= faceJ[iFace];
                s52 *= faceJ[iFace];
                s62 *= faceJ[iFace];
                s72 *= faceJ[iFace];
				s82 *= faceJ[iFace];
				s92 *= faceJ[iFace];
				s102 *= faceJ[iFace];
				s112 *= faceJ[iFace];

                shift = BASE_FUNC_COUNT;
                tmpArr2[shift + iBF] = -s22; shift += BASE_FUNC_COUNT;
                tmpArr2[shift + iBF] = -s32; shift += BASE_FUNC_COUNT;
                tmpArr2[shift + iBF] = -s42; shift += BASE_FUNC_COUNT;
                tmpArr2[shift + iBF] = -s52; shift += BASE_FUNC_COUNT;
                tmpArr2[shift + iBF] = -s62; shift += BASE_FUNC_COUNT;
                tmpArr2[shift + iBF] = -s72; shift += BASE_FUNC_COUNT;
				tmpArr2[shift + iBF] = -s82; shift += BASE_FUNC_COUNT;
				tmpArr2[shift + iBF] = -s92; shift += BASE_FUNC_COUNT;
				tmpArr2[shift + iBF] = -s102; shift += BASE_FUNC_COUNT;
				tmpArr2[shift + iBF] = -s112; //shift += BASE_FUNC_COUNT;
            }

            solverMtx->addRightElement(c1, tmpArr1);
            solverMtx->addRightElement(c2, tmpArr2);
        }
        else {
			for (int iBF = 0; iBF < BASE_FUNC_COUNT; iBF++) {
				//double s11 = 0.0;
				double s21 = 0.0;
				double s31 = 0.0;
				double s41 = 0.0;
				double s51 = 0.0;
				double s61 = 0.0;
				double s71 = 0.0;
				double s81 = 0.0;
				double s91 = 0.0;
				double s101 = 0.0;
				double s111 = 0.0;

				for (int iGP = 0; iGP < GP_FACE_COUNT; iGP++) {
					double fRO, fRU, fRV, fRW, fRE, fTAU_XX, fTAU_YY, fTAU_ZZ, fTAU_XY, fTAU_XZ, fTAU_YZ;
					double FS2, FS3, FS4, FS5, FS6, FS7, FS8, FS9, FS10, FS11;

					Point& gp = faceGP[iFace][iGP];
					double gw = faceGW[iFace][iGP];

					getFields(fRO, fRU, fRV, fRW, fRE, c1, gp);
					Param par;
					consToPar(fRO, fRU, fRV, fRW, fRE, par);
					Material& mat1 = getMaterial(c1);
					mat1.URS(par, 0); // p=p(r,e)
					mat1.getML(par);
					getTensorComponents(fTAU_XX, fTAU_YY, fTAU_ZZ, fTAU_XY, fTAU_XZ, fTAU_YZ, c1, gp);


					FS2 = fTAU_XX * n.x + fTAU_XY * n.y + fTAU_XZ * n.z;
					FS3 = fTAU_XY * n.x + fTAU_YY * n.y + fTAU_YZ * n.z;
					FS4 = fTAU_XZ * n.x + fTAU_YZ * n.y + fTAU_ZZ * n.z;
					FS5 = (par.u * fTAU_XX + par.v * fTAU_XY + par.w * fTAU_XZ) * n.x + (par.u * fTAU_XY + par.v * fTAU_YY + par.w * fTAU_YZ) * n.y + (par.u * fTAU_XZ + par.v * fTAU_YZ + par.w * fTAU_ZZ) * n.z;
					FS6 = par.ML * (4. * par.u * n.x - 2. * par.v * n.y - 2. * par.w * n.z) / 3.;
					FS7 = par.ML * (-2. * par.u * n.x + 4. * par.v * n.y - 2. * par.w * n.z) / 3.;
					FS8 = par.ML * (-2. * par.u * n.x - 2. * par.v * n.y + 4. * par.w * n.z) / 3.;
					FS9 = par.ML * (par.v * n.x + par.u * n.y);
					FS10 = par.ML * (par.w * n.x + par.u * n.z);
					FS11 = par.ML * (par.w * n.y + par.v * n.z);

					double cGP1 = gw * getF(iBF, c1, gp);
					double cGP2 = gw * getF(iBF, c2, gp);

					//s11 += 0.;
					s21 += FS2 * cGP1;
					s31 += FS3 * cGP1;
					s41 += FS4 * cGP1;
					s51 += FS5 * cGP1;
					s61 += FS6 * cGP1;
					s71 += FS7 * cGP1;
					s81 += FS8 * cGP1;
					s91 += FS9 * cGP1;
					s101 += FS10 * cGP1;
					s111 += FS11 * cGP1;

				}

				//s11 *= faceJ[iFace];
				s21 *= faceJ[iFace];
				s31 *= faceJ[iFace];
				s41 *= faceJ[iFace];
				s51 *= faceJ[iFace];
				s61 *= faceJ[iFace];
				s71 *= faceJ[iFace];
				s81 *= faceJ[iFace];
				s91 *= faceJ[iFace];
				s101 *= faceJ[iFace];
				s111 *= faceJ[iFace];

				int shift = BASE_FUNC_COUNT;
				tmpArr1[shift + iBF] = s21; shift += BASE_FUNC_COUNT;
				tmpArr1[shift + iBF] = s31; shift += BASE_FUNC_COUNT;
				tmpArr1[shift + iBF] = s41; shift += BASE_FUNC_COUNT;
				tmpArr1[shift + iBF] = s51; shift += BASE_FUNC_COUNT;
				tmpArr1[shift + iBF] = s61; shift += BASE_FUNC_COUNT;
				tmpArr1[shift + iBF] = s71; shift += BASE_FUNC_COUNT;
				tmpArr1[shift + iBF] = s81; shift += BASE_FUNC_COUNT;
				tmpArr1[shift + iBF] = s91; shift += BASE_FUNC_COUNT;
				tmpArr1[shift + iBF] = s101; shift += BASE_FUNC_COUNT;
				tmpArr1[shift + iBF] = s111; //shift += BASE_FUNC_COUNT;

			}

			if (grid.faces[iFace].bnd->faceType == CFDBoundary::TYPE_ID_WALL_NO_SLIP) {
				//for (int iBF = 0; iBF < BASE_FUNC_COUNT; iBF++) {
					//double sRO1 = 0.0;
				double sRU1 = 0.0;
				double sRV1 = 0.0;
				double sRW1 = 0.0;
				double sRE1 = 0.0;

				double fRO, fRU, fRV, fRW, fRE;
				//double FR, u, v, FE;
				double fDS = grid.faces[iFace].S;
				Point  pC = grid.faces[iFace].c;

				getFields(fRO, fRU, fRV, fRW, fRE, c1, pC);
				Param par1;
				consToPar(fRO, fRU, fRV, fRW, fRE, par1);
				Material& mat = getMaterial(c1);
				mat.getML(par1);

				// вектор скорости в ячейке
				Vector vV(par1.u, par1.v, par1.w);
				//Определяем заданную скорость вращения на регионе
				Vector vOmega(0., 0., 0.);
				Vector vOmegaR = vOmega; //vector_prod( vOmega, pC ); TODO: в двумерном пространстве векторное произведение вырождается в скаляр
				Vector vVwall = vOmegaR;
				vV -= vVwall;

				// значение скорости по нормали
				double vn = scalar_prod(vV, n);

				// нормальная компонента скорости
				Vector vVn = n;  vVn *= vn;

				// тангенциальная компонента скорости
				Vector vVt = vV;  vVt -= vVn;

				double fMU = par1.ML;			// молекулярная вязкость
				double fKP = 0.0;			// коэффициент теплопроводности TODO: добавить в код

				Vector vRP = grid.faces[iFace].c;
				vRP -= grid.cells[c1].c;			// расстояние от центра ячейки до центра грани
				double fLP_ = scalar_prod(vRP, n);		// от точки P' до центра грани

				Vector vRP_ = n;  vRP_ *= -fLP_;  vRP_ += vRP;

				fLP_ = 1.0 / fLP_;

				Vector vTauN;

				vTauN.x = -fMU * vVt.x * fLP_;
				vTauN.y = -fMU * vVt.y * fLP_;
				vTauN.z = -fMU * vVt.z * fLP_;

				double qV = scalar_prod(vTauN, vVwall);			// работа вязких сил

				double qT = 0.0;		// тепловой поток на границе

				sRU1 = vTauN.x * fDS;
				sRV1 = vTauN.y * fDS;
				sRW1 = vTauN.z * fDS;
				sRE1 = (qV + qT) * fDS;

				int shift = BASE_FUNC_COUNT;
				tmpArr1[shift/* + iBF*/] += sRU1; shift += BASE_FUNC_COUNT;
				tmpArr1[shift/* + iBF*/] += sRV1; shift += BASE_FUNC_COUNT;
				tmpArr1[shift/* + iBF*/] += sRW1; shift += BASE_FUNC_COUNT;
				tmpArr1[shift/* + iBF*/] += sRE1; //shift += BASE_FUNC_COUNT;

			//}
			}

            solverMtx->addRightElement(c1, tmpArr1);
        }
    }
}

