#include "bnd_cond.h"
#include "grid.h"

const char* CFDBoundary::TYPE_INLET_SUB         = "BOUND_INLET_SUB";
const char* CFDBoundary::TYPE_INLET_SUPER       = "BOUND_INLET_SUPER";
const char* CFDBoundary::TYPE_OUTLET_SUB        = "BOUND_OUTLET_SUB";
const char* CFDBoundary::TYPE_OUTLET_SUPER      = "BOUND_INLET_SUPER";
const char* CFDBoundary::TYPE_PRESSURE          = "BOUND_PRESSURE";
const char* CFDBoundary::TYPE_WALL_SLIP         = "BOUND_WALL_SLIP";
const char* CFDBoundary::TYPE_WALL_NO_SLIP      = "BOUND_WALL_NO_SLIP";
const char* CFDBoundary::TYPE_SYMMETRY          = "BOUND_SYMMETRY";
const char* CFDBoundary::TYPE_MASSFLOW          = "BOUND_MASSFLOW";
const char* CFDBoundary::TYPE_FREE_STREAM       = "BOUND_FREE_STREAM";
const char* CFDBoundary::TYPE_STAGNATION        = "BOUND_STAGNATION";

const int CFDBoundary::TYPE_ID_INLET_SUB        = 0;
const int CFDBoundary::TYPE_ID_INLET_SUPER      = 1;
const int CFDBoundary::TYPE_ID_OUTLET_SUB       = 2;
const int CFDBoundary::TYPE_ID_OUTLET_SUPER     = 3;
const int CFDBoundary::TYPE_ID_PRESSURE         = 4;
const int CFDBoundary::TYPE_ID_WALL_SLIP        = 5;
const int CFDBoundary::TYPE_ID_WALL_NO_SLIP     = 6;
const int CFDBoundary::TYPE_ID_SYMMETRY         = 7;
const int CFDBoundary::TYPE_ID_MASSFLOW         = 8;
const int CFDBoundary::TYPE_ID_FREE_STREAM      = 9;
const int CFDBoundary::TYPE_ID_STAGNATION       = 10;

CFDBoundary* CFDBoundary::create(TiXmlNode* bNode, Grid * g)
{
	CFDBoundary * b = nullptr;
	TiXmlNode* node = nullptr;
	TiXmlNode* node1 = nullptr;

	const char * type = bNode->FirstChild("type")->ToElement()->GetText();

	if (strcmp(type, CFDBoundary::TYPE_INLET_SUB) == 0) {
		b = new CFDBndInletSub();
        b->faceType = CFDBoundary::TYPE_ID_INLET_SUB;
		b->parCount = 6;
		b->par = new double[6];
		node = bNode->FirstChild("parameters");

		node1 = node->FirstChild("RZ");
		if (!node1) throw Exception("Parameter 'RZ' isn't specified  for BOUND_INLET_SUB.", Exception::TYPE_BOUND_NOPAR);
		node1->ToElement()->Attribute("value", &b->par[0]);

		node1 = node->FirstChild("PZ");
		if (!node1) throw Exception("Parameter 'PZ' isn't specified  for BOUND_INLET_SUB.", Exception::TYPE_BOUND_NOPAR);
		node1->ToElement()->Attribute("value", &b->par[1]);

		node1 = node->FirstChild("TZ");
		if (!node1) throw Exception("Parameter 'TZ' isn't specified  for BOUND_INLET_SUB.", Exception::TYPE_BOUND_NOPAR);
		node1->ToElement()->Attribute("value", &b->par[2]);

		node1 = node->FirstChild("CosAx");
		if (!node1) throw Exception("Parameter 'CosAx' isn't specified  for BOUND_INLET_SUB.", Exception::TYPE_BOUND_NOPAR);
		node1->ToElement()->Attribute("value", &b->par[3]);

        node1 = node->FirstChild("CosAy");
        if (!node1) throw Exception("Parameter 'CosAy' isn't specified  for BOUND_INLET_SUB.", Exception::TYPE_BOUND_NOPAR);
        node1->ToElement()->Attribute("value", &b->par[4]);

        node1 = node->FirstChild("CosAz");
        if (!node1) throw Exception("Parameter 'CosAz' isn't specified  for BOUND_INLET_SUB.", Exception::TYPE_BOUND_NOPAR);
        node1->ToElement()->Attribute("value", &b->par[5]);
	}

    if (strcmp(type, CFDBoundary::TYPE_INLET_SUPER) == 0) {
        b = new CFDBndInletSuper();
        b->faceType = CFDBoundary::TYPE_ID_INLET_SUPER;
        b->parCount = 5;
        b->par = new double[5];

        node = bNode->FirstChild("parameters");

        node1 = node->FirstChild("U");
        if (!node1) throw Exception("Parameter 'U' isn't specified  for BOUND_INLET_SUPER.", Exception::TYPE_BOUND_NOPAR);
        node1->ToElement()->Attribute("value", &b->par[0]);

        node1 = node->FirstChild("V");
        if (!node1) throw Exception("Parameter 'V' isn't specified  for BOUND_INLET_SUPER.", Exception::TYPE_BOUND_NOPAR);
        node1->ToElement()->Attribute("value", &b->par[1]);

        node1 = node->FirstChild("W");
        if (!node1) throw Exception("Parameter 'W' isn't specified  for BOUND_INLET_SUPER.", Exception::TYPE_BOUND_NOPAR);
        node1->ToElement()->Attribute("value", &b->par[2]);

        node1 = node->FirstChild("P");
        if (!node1) throw Exception("Parameter 'P' isn't specified  for BOUND_INLET_SUPER.", Exception::TYPE_BOUND_NOPAR);
        node1->ToElement()->Attribute("value", &b->par[3]);

        node1 = node->FirstChild("T");
        if (!node1) throw Exception("Parameter 'T' isn't specified  for BOUND_INLET_SUPER.", Exception::TYPE_BOUND_NOPAR);
        node1->ToElement()->Attribute("value", &b->par[4]);
    }

	if (strcmp(type, CFDBoundary::TYPE_OUTLET_SUB) == 0) {
		b = new CFDBndOutletSub();
        b->faceType = CFDBoundary::TYPE_ID_OUTLET_SUB;
        b->parCount = 1;
        b->par = new double[1];
        node = bNode->FirstChild("parameters");

        node1 = node->FirstChild("P");
        if (!node1) throw Exception("Parameter 'P' isn't specified  for TYPE_OUTLET_SUB.", Exception::TYPE_BOUND_NOPAR);
        node1->ToElement()->Attribute("value", &b->par[0]);
	}

    if (strcmp(type, CFDBoundary::TYPE_OUTLET_SUPER) == 0) {
        b = new CFDBndOutletSuper();
        b->faceType = CFDBoundary::TYPE_ID_OUTLET_SUPER;
        b->parCount = 0;
        b->par = nullptr;
    }

    if (strcmp(type, CFDBoundary::TYPE_PRESSURE) == 0) {
        b = new CFDBndPressure();
        b->faceType = CFDBoundary::TYPE_ID_PRESSURE;
        b->parCount = 2;
        b->par = new double[2];

        node = bNode->FirstChild("parameters");
        
        node1 = node->FirstChild("Pout");
        if (!node1) throw Exception("Parameter 'Pout' isn't specified  for BOUND_PRESSURE.", Exception::TYPE_BOUND_NOPAR);
        node1->ToElement()->Attribute("value", &b->par[0]);

        node1 = node->FirstChild("Tout");
        if (!node1) throw Exception("Parameter 'Tout' isn't specified  for BOUND_PRESSURE.", Exception::TYPE_BOUND_NOPAR);
        node1->ToElement()->Attribute("value", &b->par[1]); 
    }

    if (strcmp(type, CFDBoundary::TYPE_WALL_SLIP) == 0) {
        b = new CFDBndWallSlip();
        b->faceType = CFDBoundary::TYPE_ID_WALL_SLIP;
        b->parCount = 0;
        b->par = nullptr;
    }

	if (strcmp(type, CFDBoundary::TYPE_WALL_NO_SLIP) == 0) {
		b = new CFDBndWallNoSlip();
        b->faceType = CFDBoundary::TYPE_ID_WALL_NO_SLIP;
        b->parCount = 1;
        b->par = new double[1];
        node = bNode->FirstChild("parameters");

        node1 = node->FirstChild("T");
        if (!node1) throw Exception("Parameter 'T' isn't specified  for TYPE_WALL_NO_SLIP.", Exception::TYPE_BOUND_NOPAR);
        node1->ToElement()->Attribute("value", &b->par[0]);
	}

    if (strcmp(type, CFDBoundary::TYPE_SYMMETRY) == 0) {
        b = new CFDBndSymmetry();
        b->faceType = CFDBoundary::TYPE_ID_SYMMETRY;
        b->parCount = 0;
        b->par = nullptr;
    }

    if (strcmp(type, CFDBoundary::TYPE_MASSFLOW) == 0) {
        b = new CFDBndMassFlow();
        b->faceType = CFDBoundary::TYPE_ID_MASSFLOW;
        b->parCount = 6;
        b->par = new double[6];
        node = bNode->FirstChild("parameters");

        // статическое давление (сверхзвук)
        node1 = node->FirstChild("P");
        if (!node1) throw Exception("Parameter 'P' isn't specified  for BOUND_MASSFLOW.", Exception::TYPE_BOUND_NOPAR);
        node1->ToElement()->Attribute("value", &b->par[0]);

        // расход
        node1 = node->FirstChild("MassFlow");
        if (!node1) throw Exception("Parameter 'MassFlow' isn't specified  for BOUND_MASSFLOW.", Exception::TYPE_BOUND_NOPAR);
        node1->ToElement()->Attribute("value", &b->par[1]);

        // температура
        node1 = node->FirstChild("T");
        if (!node1) throw Exception("Parameter 'T' isn't specified  for BOUND_MASSFLOW.", Exception::TYPE_BOUND_NOPAR);
        node1->ToElement()->Attribute("value", &b->par[2]);

        // косинусы угла входа потока
        node1 = node->FirstChild("CosAx");
        if (!node1) throw Exception("Parameter 'CosAx' isn't specified  for BOUND_MASSFLOW.", Exception::TYPE_BOUND_NOPAR);
        node1->ToElement()->Attribute("value", &b->par[3]);

        node1 = node->FirstChild("CosAy");
        if (!node1) throw Exception("Parameter 'CosAy' isn't specified  for BOUND_MASSFLOW.", Exception::TYPE_BOUND_NOPAR);
        node1->ToElement()->Attribute("value", &b->par[4]);

        node1 = node->FirstChild("CosAz");
        if (!node1) throw Exception("Parameter 'CosAz' isn't specified  for BOUND_MASSFLOW.", Exception::TYPE_BOUND_NOPAR);
        node1->ToElement()->Attribute("value", &b->par[5]);
    }

    if (strcmp(type, CFDBoundary::TYPE_FREE_STREAM) == 0) {
        b = new CFDBndFreeStream();
        b->faceType = CFDBoundary::TYPE_ID_FREE_STREAM;
        b->parCount = 6;
        b->par = new double[6];
        node = bNode->FirstChild("parameters");

        // давление в точке торможения
        node1 = node->FirstChild("Pz");
        if (!node1) throw Exception("Parameter 'Pz' isn't specified  for BOUND_FREE_STREAM.", Exception::TYPE_BOUND_NOPAR);
        node1->ToElement()->Attribute("value", &b->par[0]);

        // температура в точке торможения
        node1 = node->FirstChild("Tz");
        if (!node1) throw Exception("Parameter 'Tz' isn't specified  for BOUND_FREE_STREAM.", Exception::TYPE_BOUND_NOPAR);
        node1->ToElement()->Attribute("value", &b->par[1]);

        // число Маха
        node1 = node->FirstChild("M");
        if (!node1) throw Exception("Parameter 'M' isn't specified  for BOUND_FREE_STREAM.", Exception::TYPE_BOUND_NOPAR);
        node1->ToElement()->Attribute("value", &b->par[2]);

        // косинусы угла входа потока
        node1 = node->FirstChild("CosAx");
        if (!node1) throw Exception("Parameter 'CosAx' isn't specified  for BOUND_FREE_STREAM.", Exception::TYPE_BOUND_NOPAR);
        node1->ToElement()->Attribute("value", &b->par[3]);

        node1 = node->FirstChild("CosAy");
        if (!node1) throw Exception("Parameter 'CosAy' isn't specified  for BOUND_FREE_STREAM.", Exception::TYPE_BOUND_NOPAR);
        node1->ToElement()->Attribute("value", &b->par[4]);

        node1 = node->FirstChild("CosAz");
        if (!node1) throw Exception("Parameter 'CosAz' isn't specified  for BOUND_FREE_STREAM.", Exception::TYPE_BOUND_NOPAR);
        node1->ToElement()->Attribute("value", &b->par[5]);
    }

    if (strcmp(type, CFDBoundary::TYPE_STAGNATION) == 0) {
        b = new CFDBndStagnation();
        b->faceType = CFDBoundary::TYPE_ID_STAGNATION;
        b->parCount = 6;
        b->par = new double[6];
        node = bNode->FirstChild("parameters");

        // статическое давление (сверхзвук)
        node1 = node->FirstChild("P");
        if (!node1) throw Exception("Parameter 'P' isn't specified  for BOUND_STAGNATION.", Exception::TYPE_BOUND_NOPAR);
        node1->ToElement()->Attribute("value", &b->par[0]);

        // полное давление (дозвук)
        node1 = node->FirstChild("Pt");
        if (!node1) throw Exception("Parameter 'Pt' isn't specified  for BOUND_STAGNATION.", Exception::TYPE_BOUND_NOPAR);
        node1->ToElement()->Attribute("value", &b->par[1]);

        // полная температура
        node1 = node->FirstChild("Tt");
        if (!node1) throw Exception("Parameter 'Tt' isn't specified  for BOUND_STAGNATION.", Exception::TYPE_BOUND_NOPAR);
        node1->ToElement()->Attribute("value", &b->par[2]);

        // косинусы угла входа потока
        node1 = node->FirstChild("CosAx");
        if (!node1) throw Exception("Parameter 'CosAx' isn't specified  for BOUND_STAGNATION.", Exception::TYPE_BOUND_NOPAR);
        node1->ToElement()->Attribute("value", &b->par[3]);

        node1 = node->FirstChild("CosAy");
        if (!node1) throw Exception("Parameter 'CosAy' isn't specified  for BOUND_STAGNATION.", Exception::TYPE_BOUND_NOPAR);
        node1->ToElement()->Attribute("value", &b->par[4]);

        node1 = node->FirstChild("CosAz");
        if (!node1) throw Exception("Parameter 'CosAz' isn't specified  for BOUND_STAGNATION.", Exception::TYPE_BOUND_NOPAR);
        node1->ToElement()->Attribute("value", &b->par[5]);
    }

    if (!b) {
		throw Exception("Unknown boundary type '%s' specified.", Exception::TYPE_BOUND_UNKNOWN);
	}

	const char * name = bNode->FirstChild("name")->ToElement()->GetText();
	strcpy(b->name, name);
	b->g = g;
	return b;
}


void CFDBndInletSub::run(int iEdge, Param& pL, Param& pR)
{
    //double fRZ = rPar[0];
    //double fPZ = rPar[1];
    ////	double fTZ	= rPar[ 2 ];
    //double fCosAx = rPar[3];
    //double fCosAy = rPar[4];
    //double fCosAz = rPar[5];

    //double fG = rP.GAMA;
    //double fQ2 = rP.Q2();
    //double fC2 = fG * fPZ / fRZ - 0.5 * fQ2 * (fG - 1.0);
    //double fRZG = pow(fRZ, fG);

    //rF.R = ::pow(fC2 * fRZG / (fPZ * fG), 1.0 / (fG - 1.0));
    //rF.P = ::pow(rF.R, fG) * fPZ / fRZG;
    //rF.E = rF.P / (rF.R * (fG - 1.0));

    //double fQ = ::sqrt(fQ2);

    //rF.U = fQ * fCosAx;
    //rF.V = fQ * fCosAy;
    //rF.W = fQ * fCosAz;

    //double Pz = par[0];
    //double Tz = par[1];
    //double cosAx = par[2];
    //double cosAy = par[3];
    //double cosAz = par[4];

    //double Rz = Pz / (pR.GAM - 1.0) / pR.e; // вычислили R

    //double gamma = pR.GAM;
    //double fQ2 = rP.Q2();
    //double fC2 = fG * fPZ / fRZ - 0.5 * fQ2 * (fG - 1.0);
    //double fRZG = ::pow(fRZ, fG);

    //rF.R = ::pow(fC2 * fRZG / (fPZ * fG), 1.0 / (fG - 1.0));
    //rF.P = ::pow(rF.R, fG) * fPZ / fRZG;
    //rF.E = rF.P / (rF.R * (fG - 1.0));

    //double fQ = ::sqrt(fQ2);

    //rF.U = fQ * fCosAx;
    //rF.V = fQ * fCosAy;
    //rF.W = fQ * fCosAz;

    //for (INDEX ic = 0; ic < m_iCKCount; ++ic)
    //    rF.CK[ic] = rPar[ID_CK[ic]];

    //g_pMat->URS(rF, 3);	// вычислили T
    //g_pMat->ML(rF);		// вычислили ML

    //CalcTP(rF, rPar);	// для расчёта MT и TP[] необходимо задать в rF ( R, ML, U, V, W )
}

void CFDBndInletSuper::run(int iEdge, Param& pL, Param& pR)
{
    pR.u = par[0];
    pR.v = par[1];
    pR.w = par[2];
    pR.p = par[3];
    pR.T = par[4];

    double gamma = pL.GAM;
    pR.cv = pL.cv;
    pR.e = pR.cv * pR.T;

    pR.r = pR.p / ((gamma - 1.0) * pR.e);

    //g_pMat->URS(rF, 5);	// вычислили R E
    //g_pMat->ML(rF);		// вычислили ML

    //CalcTP(rF, rPar);	// для расчёта MT и TP[] необходимо задать в rF ( R, ML, U, V, W )
}

void CFDBndOutletSub::run(int iEdge, Param& pL, Param& pR)
{
    double Pout = par[0];
    //	REAL fTout	= rPar[ 1 ];

    pR = pL;	// скопировали GAMA, CKCount, TPCount, ...

    //double fG = ;// rP.GAMA;
    //double q = pL.magU();
    //double fR = fQ + 2.0 * rP.CZ / (fG - 1.0);
    //double fS = rP.P / ::pow(rP.R, fG);

    //rF.P = fPout;
    //rF.R = ::pow(rF.P / fS, 1.0 / fG);
    //rF.E = rF.P / (rF.R * (fG - 1.0));

    //REAL fBeta = (fR - 2.0 / (fG - 1.0) * sqrt(fG * rF.P / rF.R)) / fQ;

    //rF.U = fBeta * rP.U;
    //rF.V = fBeta * rP.V;
    //rF.W = fBeta * rP.W;
}

//done
void CFDBndOutletSuper::run(int iEdge, Param& pL, Param& pR)
{
    pR = pL;
}

//done
void CFDBndPressure::run(int iFace, Param& pL, Param& pR)
{
    double Pout = par[0];
    double Tout = par[1];
    pR = pL;
    Vector n = g->faces[iFace].n;
    double Un = pR.u * n.x + pR.v * n.y + pR.w * n.z;
    if (Un < 0.0) { // если поток втекает   
        pR.T = Tout;
        pR.p = Pout - 0.5 * pR.r * _sqr_(Un);
    }
    else {
        if (pR.U2() < _sqr_(pR.cz)) {
            pR.p = Pout;
        }
    }
}

//done
void CFDBndWallNoSlip::run(int iFace, Param& pL, Param& pR)
{
    pR = pL;

    Vector U(pL.u, pL.v, pL.w);	// вектор скорости
    double Un = scalar_prod(U, g->faces[iFace].n);
    Vector Vn = g->faces[iFace].n;  Vn *= Un;  U -= Vn;

    pR.u = U.x;
    pR.v = U.y;
    pR.w = U.z;

    pR.T = par[0]; // задана темепература на стенке

    pR.TP[0] = 0.0; // на стенке с прилипанием КТ = 0
}

//done
void CFDBndWallSlip::run(int iFace, Param& pL, Param& pR)
{
    pR = pL;

    Vector U(pL.u, pL.v, pL.w);	// вектор скорости
    double Un = scalar_prod(U, g->faces[iFace].n);
    Vector Vn = g->faces[iFace].n;  Vn *= Un;  U -= Vn;  U -= Vn;

    pR.u = U.x;
    pR.v = U.y;
    pR.w = U.z;
}

// done
void CFDBndSymmetry::run(int iFace, Param& pL, Param& pR)
{
    pR = pL;

    Vector U(pL.u, pL.v, pL.w);	// вектор скорости
    double Un = scalar_prod(U, g->faces[iFace].n);
    Vector Vn = g->faces[iFace].n;  Vn *= Un;  U -= Vn;  U -= Vn;

    pR.u = U.x;
    pR.v = U.y;
    pR.w = U.z;
}

void CFDBndMassFlow::run(int iFace, Param& pL, Param& pR)
{
    //REAL fPout = rPar[0];	// статическое давление (сверхзвук)
    //REAL fMassFlow = rPar[1];	// расход
    //REAL fT = rPar[2];	// температура
    //REAL fCosAx = rPar[3];
    //REAL fCosAy = rPar[4];	// косинусы угла входа потока
    //REAL fCosAz = rPar[5];

    //rF = rP;	// скопировали GAMA, CKCount, TPCount, ...

    //if (rF.Q2() > _pow_2_(rF.CZ))
    //{
    //    rF.P = fPout;
    //} // if


    //rF.E = fT * rF.CV;
    //rF.R = rF.P / (rF.E * (rF.GAMA - 1.0));

    //REAL fVmag = ::fabs(fMassFlow / (rF.R * (fCosAx * vN.x + fCosAy * vN.y + fCosAz * vN.z)));	// определили модуль скорости 


    //// или

    //// Направление задано
    //rF.U = fVmag * fCosAx;
    //rF.V = fVmag * fCosAy;  // компоненты скорости
    //rF.W = fVmag * fCosAz;
}

// done
void CFDBndFreeStream::run(int iFace, Param& pL, Param& pR)
{
    double p0 = par[0];	// давление в точке торможения
    double T0 = par[1];	// температура в точке торможения
    double M = par[2];	// число Маха
    double cosAx = par[3];
    double cosAy = par[4];	// косинусы угла входа потока
    double cosAz = par[5];

    pR = pL;
    double r0 = p0 / ((pL.GAM - 1.0) * pL.e);

    double gamma = pL.GAM;

    pR.GAM = gamma;

    double CZ0 = ::sqrt(gamma * p0 / r0);
    double q0 = M * CZ0;

    double U0 = q0 * cosAx;
    double V0 = q0 * cosAy;
    double W0 = q0 * cosAz;

    Vector n = g->faces[iFace].n;

    double UN0 = U0 * n.x + V0 * n.y + W0 * n.z;

    if (M < 1.0)
    {
        double UNF = pL.u * n.x + pL.v * n.y + pL.w * n.z;

        double I0 = UN0 - 2.0 * CZ0 / (gamma - 1.0);
        double IP = UNF + 2.0 * pL.cz / (gamma - 1.0);

        double UN = 0.5 * (IP + I0);
        pR.cz = 0.25 * (IP - I0) * (gamma - 1.0);

        double RG_P = ::pow(r0, gamma) / p0;
        pR.r = ::pow(pR.cz * pR.cz * RG_P / gamma, 1.0 / (gamma - 1.0));
        pR.p = ::pow(pR.r, gamma) / RG_P;

        pR.e = pR.p / (pR.r * (gamma - 1.0));

        static Vector UT;
        if (UN0 < 0.0) { // если поток втекает

            UT.x = U0 - UN0 * n.x;
            UT.y = V0 - UN0 * n.y;
            UT.z = W0 - UN0 * n.z;

            //for (INDEX ic = 0; ic < m_iCKCount; ++ic)  rF.CK[ic] = rPar[ID_CK[ic]];

        }
        else {			 // если поток вытекает

            UT.x = pL.u - UNF * n.x;
            UT.y = pL.v - UNF * n.y;
            UT.z = pL.w - UNF * n.z;

        } // if

        double  fUT = _mag_(UT);

        double q = ::sqrt(_sqr_(UN) + _sqr_(fUT));

        double K = q / q0;

        pR.u = U0 * K;
        pR.v = V0 * K;
        pR.w = W0 * K;

    }
    else {

        if (UN0 < 0.0) // сверхзвуковой вход
        {
            pR.u = U0;
            pR.v = V0;
            pR.w = W0;

            pR.p = p0;
            pR.r = r0;
            pR.e = pR.p / (pR.r * (gamma - 1.0));

            //			for ( INDEX ic=0; ic<m_iCKCount; ++ic )  rF.CK[ ic ] = rPar[ ID_CK[ ic ] ];

        } // if

    } // if
}

//done
void CFDBndStagnation::run(int iFace, Param& pL, Param& pR)
{
    double Pout = par[0];	// статическое давление (сверхзвук)
    double Pt = par[1];	// полное давление (дозвук)
    double Tt = par[2];	// полная температура
    double cosAx = par[3];
    double cosAy = par[4];	// косинусы угла входа потока
    double cosAz = par[5];

    pR = pL;	// скопировали GAMA, CKCount, TPCount, ...

    if (pR.U2() > _sqr_(pR.cz))
    {
        pR.p = Pout;
    } // if

    pR.T = Tt / ::pow(Pt / pR.p, (pR.CP - pR.cv) / pR.CP); // определили температуру
    pR.E = pR.T * pR.cv;
    pR.r = pR.p / (pR.E * (pR.GAM - 1.0));

    if (Tt < pR.T) pR.T = Tt;

    double Vmag = ::sqrt(2.0 * pR.CP * (Tt - pR.T));	// определили модуль скорости

    pR.u = Vmag * cosAx;
    pR.v = Vmag * cosAy;  // компоненты скорости
    pR.w = Vmag * cosAz;
}
