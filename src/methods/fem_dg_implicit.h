#pragma once
#include "fem_dg.h"
#include "MatrixSolver.h"
#include "Limiter.h"

class FEM_DG_IMPLICIT : public FEM_DG
{
	friend class Limiter;
public:
	void init(char * xmlFileName);
	void run();
	void done();

	virtual void getFields(double &fRO, double &fRU, double &fRV, double &fRW, double &fRE, int iCell, Point p);
	void getFieldsLinear(double& fRO, double& fRU, double& fRV, double& fRW, double& fRE, int iCell, double x, double y, double z);
	void getFieldsLinear(double& fRO, double& fRU, double& fRV, double& fRW, double& fRE, int iCell, Point p);
	virtual void getFields(double &fRO, double &fRU, double &fRV, double &fRW, double &fRE, int iCell, double x, double y, double z);

	void getTensorComponents(double& fTAU_XX, double& fTAU_YY, double& fTAU_ZZ, double& fTAU_XY, double& fTAU_XZ, double& fTAU_YZ, int iCell, Point p);
	void getTensorComponents(double& fTAU_XX, double& fTAU_YY, double& fTAU_ZZ, double& fTAU_XY, double& fTAU_XZ, double& fTAU_YZ, int iCell, double x, double y, double z);

	virtual double getFieldLinear(int fld, int iCell, Point p);
	virtual double getField(int fld, int iCell, Point p);
	virtual double getFieldLinear(int fldId, int iCell, double x, double y, double z);
	virtual double getField(int fld, int iCell, double x, double y, double z);

	virtual double getF(int id, int iCell, Point p);
	virtual double getF(int id, int iCell, double x, double y, double z);

private:
	void memAlloc();
	void memFree();

	void calcMassMatr(); //!< вычисляем матрицу масс
	void calcGaussPar(); //!< вычисляем узлы и коэффициенты квадратур

	void calcTimeStep();
	
	Region & getRegionByCellType(int type);
	Region& getRegionByName(char* name);

	Region   &	getRegion(int iCell);
	Region& getRegion(char* name);
	
	Material &	getMaterial(int iCell);

	inline void convertParToCons(int iCell, Param & par); //!< Преобразование примитивных переменных в консервативные

	inline void convertConsToPar(int iCell, Param & par); //!< Преобразование консервативных переменных в примитивные


	inline double getDfDx(int id, int iCell, Point p);
	inline double getDfDx(int id, int iCell, double x, double y, double z);

	inline double getDfDy(int id, int iCell, Point p);
	inline double getDfDy(int id, int iCell, double x, double y, double z);

	inline double getDfDz(int id, int iCell, Point p);
	inline double getDfDz(int id, int iCell, double x, double y, double z);

	void save(int step);

	inline void setCellFlagLim(int iCell){ grid.cells[iCell].flag |= CELL_FLAG_LIM; }
	inline bool cellIsLim(int iCell)		{ return (grid.cells[iCell].flag & CELL_FLAG_LIM) > 0; }
	int getLimitedCellsCount();
	void remediateLimCells();

	void incCFL();
	void decCFL();

	void calcIntegral();			//!< вычисляем интеграл от(dF / dU)*deltaU*dFi / dx
	void calcMatrWithTau();			//!< вычисляем матрицы перед производной по времени
	void calcMatrFlux();			//!< Вычисляем потоковые величины
	void calcRHS();					//!< Вычисляем столбец правых членов
	void calcMatrTensor();			//!< Вычисляем матрицы перед компонентами тензора вязких напряжений
	void calcViscousIntegral(); 	//!< Вычисляем интеграл от (dH / dU)*dFi / dx
	void calcMatrViscousFlux();	//!< Вычисляем потоковые величины от диффузионных членов
    void calcViscousRHS();					//!< Вычисляем столбец правых членов

	void calcLiftForce();

	double** allocMtx5();
	void freeMtx5(double **mtx5);
	void multMtx5(double **dst5, double **srcA5, double **srcB5);
	void clearMtx5(double **mtx5);
	double** allocMtx11();
	void freeMtx11(double **mtx11);
	void multMtx11(double **dst11, double **srcA11, double **srcB11);
	void clearMtx11(double **mtx11);
	void multMtxToVal(double **dst, double x, int N);
	void fillMtx(double** dst, double x, int N);

	void eigenValues(double** dst5, double c, double u, double nx, double v, double ny, double w, double nz);
	void eigenValuesAbs(double** dst5, double c, double u, double nx, double v, double ny, double w, double nz);
	void rightEigenVector(double **dst4, double c, double u, double nx, double v, double ny, double w, double nz, double H);
	void leftEigenVector(double **dst4, double c, double GAM, double u, double nx, double v, double ny, double w, double nz);
	void calcAP(double **dst5, double **rightEgnVecl5, double **egnVal5, double **leftEgnVecl5);
	void calcAM(double **dst5, double **rightEgnVecl5, double **egnVal5, double **leftEgnVecl5);
	void calcRoeAverage(Param& average, Param pL, Param pR, double GAM, Vector n);
	void _calcA(double **dst5, double **rightEgnVecl5, double **egnVal5, double **leftEgnVecl5);
	void calcA(double **dst5, double c, double GAM, double u, double nx, double v, double ny, double w, double nz, double H);
	void calcAbsA(double** dst5, double c, double GAM, double u, double nx, double v, double ny, double w, double nz, double H);
    void calcJ(double **dst11, double r, double u, double nx, double v, double ny, double w, double nz, double mu, double tau_xx, double tau_yy, double tau_zz, double tau_xy, double tau_xz, double tau_yz);
	void calcAx(double **dst5, double c, double GAM, double u, double v, double w, double H);
	void calcAy(double **dst5, double c, double GAM, double u, double v, double w, double H);
	void calcAz(double** dst5, double c, double GAM, double u, double v, double w, double H);
	void calcAx_(double **dst5, Param par, double GAM);
	void calcAy_(double **dst5, Param par, double GAM);
	void calcAz_(double** dst5, Param par, double GAM);
	void consToPar(double fRO, double fRU, double fRV, double fRW, double fRE, Param& par);
	void calcGeomPar();
	void calcGrad();
	void calcTurbulent();

	void calcFlux(double& fr, double& fu, double& fv, double& fw, double& fe, Param pL, Param pR, Vector n, double GAM);
	void calcFluxFields(Param& pF, Param pL, Param pR, Vector n, double GAM);
	void boundaryCond(int iFace, Param& pL, Param& pR);

	void addSmallMatrToBigMatr(double **mB, double **mS, int i, int j);
	void defineOrientationOnFaces();

private:
	double			TMAX;
	int				STEP_MAX;
	double			TAU;
	double			CFL;
	double			scaleCFL;
	double			maxCFL;
	int				stepCFL;
	double			maxLimCells;
	int				FILE_SAVE_STEP;
	int				PRINT_STEP;

	double			TAU_MIN;

	bool			STEADY;	// false - нестационарное течение, true - стационанрное течение.
	double			*cTau;  // локальный шаг по времени в ячейке.
	bool			SMOOTHING;
	double			SMOOTHING_PAR;
	int				FLUX;
	bool            VISCOSITY;
	bool            TURBULENCE; // false - ламинарное течение, true - турбулентное течение.
	int				TURBULENCE_MODEL;

	double			SIGMA;

	double			SOLVER_EPS		= 1.e-7;
	int				SOLVER_ITER		= 50;

	int				matCount;
	int				regCount;
	int				bCount;
	Material		*materials;
	Region			*regions;
	CFDBoundaries	boundaries;


	double			***fields;

	double			*tmpArr;
	double			*tmpArr1;
	double			*tmpArr2;
	int				*tmpArrInt;

	//! градиенты.
	Vector			*gradR;
	Vector			*gradP;
	Vector			*gradU;
	Vector			*gradV;
	Vector			*gradW;

	//! лимиты
	double			limitRmin;
	double			limitRmax;
	double			limitPmin;
	double			limitPmax;
	double			limitUmax;

	//! стандартные значения
	double			RRef;
	double			PRef;

	int				limCells;

	//! подъемная сила.
	double			Fx;
	double			Fy;
	double			Fz;

	Vector			beta; // фиксированный вектор

	//! узлы и коэффициенты квадратур
	/*Point			**cellGP;
	Point			**faceGP;
	double			**cellGW;
	double			*cellJ;
	double			**faceGW;
	double			*faceJ;*/

	double* faceLP_;	// массив растояний от точки P' до центра грани
	double* faceLE_;	// массив растояний от точки E' до центра грани
	Vector* faceRP;		// массив векторов  от цента ячейки 1(P) до центра грани
	Vector* faceRE;		// массив векторов  от цента ячейки 2(E) до центра грани

	// матрица масс
	double			***matrA;
	double			***matrInvA;

	MatrixSolver	*solverMtx;

	double			**matrBig;
	double			**matrSmall;
	double			**matrBig2;
	double			**matrSmall2;

	double			*tmpCFL;
	double			pCFL;

	Limiter			*limiter = NULL;

	// параметры обезразмеривания
	double			L_		= 1.0;
	double			U_		= 1.0;
	double			R_		= 1.0;
	double			P_		= 1.0;
	double			T_		= 1.0;
	double			E_		= 1.0;
	double			CV_ = 1.0;
	double			MU_ = 1.0;
	double			KP_ = 1.0;
	double			TIME_ = 1.0;
    double			TAU_ = 1.0;

	int				fieldCount;
	int				matrDim;

protected:
	//! параметры модели турбулентности SA
	double    SigmaNT = 0.6667;
	double    Cb1 = 0.1355;
	double    Cb2 = 0.622;
	double    Cv1 = 7.1;
	double    Cw2 = 0.3;
	double    Cw3 = 2.0;
	double    K = 0.41;
	double    Ct3 = 1.2;
	double    Ct4 = 0.5;
	double	  Cprod = 2.0;

	const static int MODEL_SA = 0;

	const static int FLUX_GODUNOV = 0;
	const static int FLUX_LAX = 1;
	const static int FLUX_CD = 2;
	const static int FLUX_HLLC = 3;
	const static int FLUX_AUSMP = 4;
	const static int FLUX_AUSMPW = 5;
	const static int FLUX_KIR = 6;

	//const static int GP_CELL_COUNT = 4;
	//const static int GP_FACE_COUNT = 3;

	const static int FIELD_COUNT = 5;
	const static int FIELD_COUNT_EXT = 11;
	const static int FIELD_RO = 0;
	const static int FIELD_RU = 1;
	const static int FIELD_RV = 2;
	const static int FIELD_RW = 3;
	const static int FIELD_RE = 4;
	const static int FIELD_TAU_XX = 5;
	const static int FIELD_TAU_YY = 6;
	const static int FIELD_TAU_ZZ = 7;
	const static int FIELD_TAU_XY = 8;
	const static int FIELD_TAU_XZ = 9;
	const static int FIELD_TAU_YZ = 10;

	const static int MATR_DIM = FIELD_COUNT_EXT * BASE_FUNC_COUNT;

	inline double getGAM(int iCell) { return getMaterial(iCell).getGamma(); } // TODO: сделать
};

