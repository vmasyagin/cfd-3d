#include "fvm_tvd.h"
#include "tinyxml.h"
#include <string>
#include "global.h"
#include "MeshReader.h"

void FVM_TVD::init(char * xmlFileName)
{
	STEADY = true;
	
	
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
	node0->FirstChild("CFL")->ToElement()->Attribute("value", &CFL);
	node0->FirstChild("STEP_MAX")->ToElement()->Attribute("value", &STEP_MAX);
	node0->FirstChild("FILE_OUTPUT_STEP")->ToElement()->Attribute("value", &FILE_SAVE_STEP);
	node0->FirstChild("LOG_OUTPUT_STEP")->ToElement()->Attribute("value", &PRINT_STEP);

	if (steadyVal == 0) {
		STEADY = false;
	} else {
		STEADY = true;
	}

	// чтение параметров о ПРЕДЕЛЬНЫХ ЗНАЧЕНИЯХ
	node0 = task->FirstChild("limits");
	node0->FirstChild("ro")->ToElement()->Attribute("min", &limitRmin);
	node0->FirstChild("ro")->ToElement()->Attribute("max", &limitRmax);
	node0->FirstChild("p")->ToElement()->Attribute( "min", &limitPmin);
	node0->FirstChild("p")->ToElement()->Attribute( "max", &limitPmax);
	node0->FirstChild("u")->ToElement()->Attribute( "max", &limitUmax);

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

		boundaries.push_back(b);

		bNode = bNode->NextSibling("boundCond");
	}

	bCount = boundaries.size();

	/* Чтение данных сетки. */
	node0 = task->FirstChild("mesh");
	const char* fName = node0->FirstChild("name")->ToElement()->Attribute("value");
	const char* tName = node0->FirstChild("filesType")->ToElement()->Attribute("value");
	MeshReader* mr = MeshReader::create(MeshReader::getType((char*)tName), (char*)fName);
	mr->read(&grid);


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


	cTau = new double[grid.cCount];

	ro		= new double[grid.cCount];
	ru		= new double[grid.cCount];
	rv		= new double[grid.cCount];
	re		= new double[grid.cCount];

	ro_old	= new double[grid.cCount];
	ru_old	= new double[grid.cCount];
	rv_old	= new double[grid.cCount];
	re_old	= new double[grid.cCount];

	ro_int	= new double[grid.cCount];
	ru_int	= new double[grid.cCount];
	rv_int	= new double[grid.cCount];
	re_int	= new double[grid.cCount];

	gradR		= new Vector[grid.cCount];
	gradP		= new Vector[grid.cCount];
	gradU		= new Vector[grid.cCount];
	gradV		= new Vector[grid.cCount];

	for (int i = 0; i < grid.cCount; i++)
	{
		Cell & c = grid.cells[i];
		Region & reg = getRegion(c.typeName);
		convertParToCons(i, reg.par);
	}

	memcpy(ro_old, ro, grid.cCount*sizeof(double));
	memcpy(ru_old, ru, grid.cCount*sizeof(double));
	memcpy(rv_old, rv, grid.cCount*sizeof(double));
	memcpy(re_old, re, grid.cCount*sizeof(double));

	calcTimeStep();
	save(0);
}


void FVM_TVD::calcTimeStep()
{
	if (STEADY) {
		for (int iCell = 0; iCell < grid.cCount; iCell++)
		{
			Param p;
			convertConsToPar(iCell, p);
			cTau[iCell] = CFL * grid.cells[iCell].S / _max_(abs(p.u) + p.cz, abs(p.v) + p.cz);
		}
	}
	else {
		for (int iCell = 0; iCell < grid.cCount; iCell++)
		{
			Param p;
			convertConsToPar(iCell, p);
			if (TAU > CFL* grid.cells[iCell].S / _max_(abs(p.u) + p.cz, abs(p.v) + p.cz)) {
				TAU = CFL * grid.cells[iCell].S / _max_(abs(p.u) + p.cz, abs(p.v) + p.cz);
			}
		}
		for (int iCell = 0; iCell < grid.cCount; iCell++) {
			cTau[iCell] = TAU;
		}
		log("time step: %25.16E\n", TAU);
	}
}

void FVM_TVD::calcGrad() 
{
	int nc = grid.cCount;
	int ne = grid.eCount;
	

	memset(gradR, 0, nc*sizeof(Vector));
	memset(gradP, 0, nc*sizeof(Vector));
	memset(gradU, 0, nc*sizeof(Vector));
	memset(gradV, 0, nc*sizeof(Vector));
	
    for (int iEdge = 0; iEdge < ne; iEdge++)
	{
			
		int c1	= grid.edges[iEdge].c1;
		int c2	= grid.edges[iEdge].c2;
			
		Param pL, pR;
		convertConsToPar(c1, pL);
		if (c2 > -1) {
			convertConsToPar(c2, pR);
		} else {
			boundaryCond(iEdge, pL, pR);
		}

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
		gradR[iCell] /= si;
		gradP[iCell] /= si;
		gradU[iCell] /= si;
		gradV[iCell] /= si;
	}
}

void FVM_TVD::run() 
{
	int nc = grid.cCount;
	int ne = grid.eCount;

	double			t		= 0.0;
	unsigned int	step	= 0;
	while (t < TMAX && step < STEP_MAX) 
	{
		if (!STEADY) {
			t += TAU; 
		} else {
			calcTimeStep();
		}
		step++;
		memcpy(ro_old, ro, nc*sizeof(double));
		memcpy(ru_old, ru, nc*sizeof(double));
		memcpy(rv_old, rv, nc*sizeof(double));
		memcpy(re_old, re, nc*sizeof(double));

		// первый подшаг метода Р.-К.
		memset(ro_int, 0, nc*sizeof(double));
		memset(ru_int, 0, nc*sizeof(double));
		memset(rv_int, 0, nc*sizeof(double));
		memset(re_int, 0, nc*sizeof(double));
		calcGrad();
		for (int iEdge = 0; iEdge < ne; iEdge++)
		{
			double fr, fu, fv, fe;
			int c1	= grid.edges[iEdge].c1;
			int c2	= grid.edges[iEdge].c2;
			Vector n	= grid.edges[iEdge].n;
			double l	= grid.edges[iEdge].l*0.5;
			Param pL, pR;
			fr = 0.0;
			fu = 0.0;
			fv = 0.0;
			fe = 0.0;
			for (int iGP = 1; iGP < grid.edges[iEdge].cCount; iGP++)
			{
				double fr1, fu1, fv1, fe1;
				reconstruct(iEdge, pL, pR, grid.edges[iEdge].c[iGP]);
				double __GAM = 1.4; // TODO: сделать правильное вычисление показателя адиабаты
				calcFlux(fr1, fu1, fv1, fe1, pL, pR, n, __GAM);
				fr += fr1;
				fu += fu1;
				fv += fv1;
				fe += fe1;

			}
			ro_int[c1] -= fr*l;
			ru_int[c1] -= fu*l;
			rv_int[c1] -= fv*l;
			re_int[c1] -= fe*l;
			if (c2 > -1) 
			{
				ro_int[c2] += fr*l;
				ru_int[c2] += fu*l;
				rv_int[c2] += fv*l;
				re_int[c2] += fe*l;
			}

		}
		for (int iCell = 0; iCell < nc; iCell++)
		{
			if (cellIsLim(iCell)) continue;
			register double cfl = cTau[iCell]/grid.cells[iCell].S;
			ro[iCell] += cfl*ro_int[iCell];
			ru[iCell] += cfl*ru_int[iCell];
			rv[iCell] += cfl*rv_int[iCell];
			re[iCell] += cfl*re_int[iCell];
		}

		// второй подшаг метода Р.-К.
		memset(ro_int, 0, nc*sizeof(double));
		memset(ru_int, 0, nc*sizeof(double));
		memset(rv_int, 0, nc*sizeof(double));
		memset(re_int, 0, nc*sizeof(double));
		calcGrad();
		for (int iEdge = 0; iEdge < ne; iEdge++)
		{
			double fr, fu, fv, fe;
			int c1	= grid.edges[iEdge].c1;
			int c2	= grid.edges[iEdge].c2;
			Vector n	= grid.edges[iEdge].n;
			double l	= grid.edges[iEdge].l*0.5;
			Param pL, pR;
			fr = 0.0;
			fu = 0.0;
			fv = 0.0;
			fe = 0.0;
			for (int iGP = 1; iGP < grid.edges[iEdge].cCount; iGP++) 
            {
				double fr1, fu1, fv1, fe1;
				reconstruct(iEdge, pL, pR, grid.edges[iEdge].c[iGP]);
				double __GAM = 1.4; // TODO: сделать правильное вычисление показателя адиабаты
				calcFlux(fr1, fu1, fv1, fe1, pL, pR, n, __GAM);
				fr += fr1;
				fu += fu1;
				fv += fv1;
				fe += fe1;

			}
			ro_int[c1] -= fr*l;
			ru_int[c1] -= fu*l;
			rv_int[c1] -= fv*l;
			re_int[c1] -= fe*l;
			if (c2 > -1) 
			{
				ro_int[c2] += fr*l;
				ru_int[c2] += fu*l;
				rv_int[c2] += fv*l;
				re_int[c2] += fe*l;
			}

		}
		for (int iCell = 0; iCell < nc; iCell++)
		{
			if (cellIsLim(iCell)) continue;
			register double cfl = cTau[iCell]/grid.cells[iCell].S;
			ro[iCell] += cfl*ro_int[iCell];
			ru[iCell] += cfl*ru_int[iCell];
			rv[iCell] += cfl*rv_int[iCell];
			re[iCell] += cfl*re_int[iCell];
		}

		// полусумма: формула (4.10) из icase-1997-65.pdf
		for (int iCell = 0; iCell < nc; iCell++)
		{
			if (cellIsLim(iCell)) continue;

			ro[iCell] = 0.5*(ro_old[iCell]+ro[iCell]);
			ru[iCell] = 0.5*(ru_old[iCell]+ru[iCell]);
			rv[iCell] = 0.5*(rv_old[iCell]+rv[iCell]);
			re[iCell] = 0.5*(re_old[iCell]+re[iCell]);

			Param par;
			convertConsToPar(iCell, par);
			if (par.r < limitRmin)			{ par.r = limitRmin; setCellFlagLim(iCell); }
			if (par.r > limitRmax)			{ par.r = limitRmax; setCellFlagLim(iCell); }
			if (par.p < limitPmin)			{ par.p = limitPmin; setCellFlagLim(iCell); }
			if (par.p > limitPmax)			{ par.p = limitPmax; setCellFlagLim(iCell); }
			if (fabs(par.u) > limitUmax)	{ par.u = limitUmax; setCellFlagLim(iCell); }
			if (fabs(par.v) > limitUmax)	{ par.v = limitUmax; setCellFlagLim(iCell); }
		}

		remediateLimCells();


		if (step % FILE_SAVE_STEP == 0)
		{
			save(step);
		}
		if (step % PRINT_STEP == 0)
		{
			log("step: %d\t\ttime step: %.16f\n", step, t);
		}
	}

}

void FVM_TVD::remediateLimCells()
{
	for (int iCell = 0; iCell < grid.cCount; iCell++) 
	{
		if (cellIsLim(iCell)) 
		{
			// пересчитываем по соседям
			double sRO = 0.0;
			double sRU = 0.0;
			double sRV = 0.0;
			double sRE = 0.0;
			double S   = 0.0;
			for (int i = 0; i < grid.cells[iCell].eCount; i++)
			{
				int iEdge = grid.cells[iCell].edgesInd[i];
				int j = grid.edges[iEdge].c2;
				if (j >= 0) {
					double s = grid.cells[j].S;
					S += s;
					sRO += ro[j] * s;
					sRU += ru[j] * s;
					sRV += rv[j] * s;
					sRE += re[j] * s;
				} 
			}
			ro[iCell] = sRO/S;
			ru[iCell] = sRU/S;
			rv[iCell] = sRV/S;
			re[iCell] = sRE/S;

			// после 0x20 итераций пробуем вернуть ячейку в счет
			grid.cells[iCell].flag += 0x010000;
			if (grid.cells[iCell].flag & 0x200000) grid.cells[iCell].flag &= 0x001110;
		}
	}
}

void FVM_TVD::save(int step)
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
		fprintf(fp, "%25.16f ", (p.u*p.u+p.v*p.v)/p.cz);
		if (i+1 % 8 == 0 || i+1 == grid.cCount) fprintf(fp, "\n");
	}

	fprintf(fp, "VECTORS Velosity float\n");
	for (int i = 0; i < grid.cCount; i++)
	{
		Param p;
		convertConsToPar(i, p);
		fprintf(fp, "%25.16f %25.16f %25.16f ", p.u, p.v, 0.0);
		if (i+1 % 8 == 0 || i+1 == grid.cCount) fprintf(fp, "\n");
	}

	fprintf(fp, "SCALARS Total_pressure float 1\nLOOKUP_TABLE default\n");
	for (int i = 0; i < grid.cCount; i++)
	{
		Material &mat = getMaterial(i);
		double gam = mat.getGamma();
		double agam = gam - 1.0;
		Param p;
		convertConsToPar(i, p);
		double M2 = (p.u*p.u+p.v*p.v)/(gam*p.p/p.r);
		fprintf(fp, "%f ", p.p*::pow(1.0+0.5*M2*agam, gam/(gam-1.0)) );
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

void FVM_TVD::calcFlux(double& fr, double& fu, double& fv, double& fe, Param pL, Param pR, Vector n, double GAM)
{
	{	// GODUNOV FLUX
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
// 	{	// LAX-FRIEDRIX FLUX
// 		double unl = pL.u*n.x+pL.v*n.y;
// 		double unr = pR.u*n.x+pR.v*n.y;
// 		double rol, rul, rvl, rel,  ror, rur, rvr, rer;
// 		double alpha = _max_(fabs(unl)+sqrt(GAM*pL.p/pL.r), fabs(unr)+sqrt(GAM*pR.p/pR.r));
// 		rol = pL.r;
// 		rul = pL.r*pL.u;
// 		rvl = pL.r*pL.v;
// 		rel = pL.r*pL.E;
// 		ror = pR.r;
// 		rur = pR.r*pR.u;
// 		rvr = pR.r*pR.v;
// 		rer = pR.r*pR.E;
// 		double frl = rol*unl;
// 		double frr = ror*unr;
// 		fr = 0.5*(frr+frl								- alpha*(ror-rol));
// 		fu = 0.5*(frr*pR.u+frl*pL.u + (pR.p+pL.p)*n.x	- alpha*(rur-rul));
// 		fv = 0.5*(frr*pR.v+frl*pL.v + (pR.p+pL.p)*n.y	- alpha*(rvr-rvl));
// 		fe = 0.5*((rer+pR.p)*unr + (rel+pL.p)*unl		- alpha*(rer-rel));
// 	}
}


void FVM_TVD::reconstruct(int iEdge, Param& pL, Param& pR, Point p)
{
	if (grid.edges[iEdge].type == Edge::TYPE_INNER) 
	{
		int c1	= grid.edges[iEdge].c1;
		int c2	= grid.edges[iEdge].c2;
		convertConsToPar(c1, pL);
		convertConsToPar(c2, pR);
		//return;
		//Point PE = grid.edges[iEdge].c[0];
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
		int c1 = grid.edges[iEdge].c1;
		convertConsToPar(c1, pL);
		//return;
		Point &PE = p;
		Point P1 = grid.cells[c1].c;
		Vector DL1;
		DL1.x = PE.x - P1.x;
		DL1.y = PE.y - P1.y;
		pL.r += gradR[c1].x*DL1.x + gradR[c1].y*DL1.y;
		pL.p += gradP[c1].x*DL1.x + gradP[c1].y*DL1.y;
		pL.u += gradU[c1].x*DL1.x + gradU[c1].y*DL1.y;
		pL.v += gradV[c1].x*DL1.x + gradV[c1].y*DL1.y;

		boundaryCond(iEdge, pL, pR);
	}
}


void FVM_TVD::boundaryCond(int iEdge, Param& pL, Param& pR)
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
	

	//Edge &edge = grid.edges[iEdge];
	//int c1	= grid.edges[iEdge].c1;
	//Material& m = getMaterial(c1);
	//if (edge.bnd) {
	//	edge.bnd->run(iEdge, pL, pR);
	//	m.URS(pR, 2);
	//	m.URS(pR, 1);
	//	return;
	//}
	//else {
	//	char msg[128];
	//	sprintf(msg, "Not defined boundary condition for edge %d\n", iEdge);
	//	throw Exception(msg, Exception::TYPE_BOUND_UNKNOWN);
	//}

}


void FVM_TVD::done()
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

	delete[] gradR;
	delete[] gradP;
	delete[] gradU;
	delete[] gradV;
}




Region & FVM_TVD::getRegionByCellType(int type)
{
	for (int i = 0; i < regCount; i++)
	{
		if (regions[i].cellType == type) return regions[i];
	}
	log("ERROR: unknown cell type %d...\n", type);
	EXIT(1);
}


Region   &	FVM_TVD::getRegion	(int iCell)
{
	return getRegionByCellType( grid.cells[iCell].type );
}

Region & FVM_TVD::getRegionByName(char* name)
{
	for (int i = 0; i < regCount; i++)
	{
		if (strcmp(regions[i].name.c_str(), name) == 0) return regions[i];
	}
	log("ERROR: unknown cell name '%d'...\n", name);
	EXIT(1);
}

Region & FVM_TVD::getRegion(char * name)
{
	return getRegionByName(name);
}

Material &	FVM_TVD::getMaterial(int iCell)
{
	Region & reg = getRegion(grid.cells[iCell].typeName);
	return materials[reg.matId];
}


void FVM_TVD::convertParToCons(int iCell, Param & par)
{
	ro[iCell] = par.r;
	ru[iCell] = par.r*par.u;
	rv[iCell] = par.r*par.v;
	re[iCell] = par.r*(par.e+0.5*(par.u*par.u+par.v*par.v));
}

void FVM_TVD::convertConsToPar(int iCell, Param & par)
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



