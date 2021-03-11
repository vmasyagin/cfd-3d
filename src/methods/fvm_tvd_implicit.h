#ifndef _FVM_TVD_IMPLICIT_H_
#define _FVM_TVD_IMPLICIT_H_

#include <method.h>
#include "MatrixSolver.h"
#include <algorithm>

class FVM_TVD_IMPLICIT: public Method
{
public:
	virtual void init(char * xmlFileName);
	virtual void run();
	virtual void done();
protected:
	Region & getRegionByCellType(int type);

	Region   &	getRegion	(int iCell);
	Material &	getMaterial	(int iCell);
	
	/**
	 *	ѕреобразование примитивных переменных в консервативные
	 */
	void convertParToCons(int iCell, Param & par);
	
	/**
	 *	ѕреобразование консервативных переменных в примитивные
	 */
	void convertConsToPar(int iCell, Param & par);
	
	/**
	 *	¬ычисление параметров с внешней стороны от границы €чейки согласно граничным услови€м
	 */
	void boundaryCond(int iFace, Param& pL, Param& pR);

	/**
	 *	¬ычисление шага по времени по значению CFL, если значение TAU из XML 
	 *	меньше вычисленного, то используетс€ значение, заданное в XML
	 */
	void calcTimeStep();

	/**
	 *	«апись значений газодинамических параметров в файл
	 */
	void save(int);

	/**
	 *	¬ычисление численного потока
	 */
	void calcFlux(double& fr, double& fu, double& fv, double& fe, Param pL, Param pR, Vector n, double GAM);

	void reconstruct(int iFace, Param& pL, Param& pR);
	void reconstruct(int iFace, Param& pL, Param& pR, Point p);

	void calcGrad();
private:
	double **allocMtx4();													//+.
	void freeMtx4(double **mtx4);											//+.
	void multMtx4(double **dst4, double **srcA4, double **srcB4);			//+.
	void clearMtx4(double **mtx4);											//+.
	void printMtx4(double **mtx4, char *msg = 0);							//+.
	
	// «аполн€ет в mtx4 собственные значени€ матрицы якоби функции потока.
	void eigenValues(double **dst4, double c, double u, double nx, double v, double ny);					//+.
	void rightEigenVector(double **dst4, double c, double u, double nx, double v, double ny, double H);		//+.
	void leftEigenVector(double **dst4, double c, double GAM, double u, double nx, double v, double ny);	//+.
	
	// «аполн€ет в mtx4 A+.
	void calcAP(double **dst4, double **rightEgnVecl4, double **egnVal4, double **leftEgnVecl4);			//+.
	// «аполн€ет в mtx4 A-.
	void calcAM(double **dst4, double **rightEgnVecl4, double **egnVal4, double **leftEgnVecl4);			//+.

	void calcRoeAverage(Param& average, Param pL, Param pR, double GAM, Vector n);	//+.
	
	void setCellFlagLim(int iCell)	{ grid.cells[iCell].flag |= CELL_FLAG_LIM; }
	bool cellIsLim(int iCell)		{ return (grid.cells[iCell].flag & CELL_FLAG_LIM) > 0; }
	int getLimitedCellsCount();
	void remediateLimCells();

	void smoothingDelta(double*);

	void incCFL();
	void decCFL();

	void calcLiftForce();
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
	double			*cTau;  // локальный шаг по времени в €чейке.
	bool			SMOOTHING;
	double			SMOOTHING_PAR;
	int				FLUX;


	int				matCount;
	int				regCount;
	int				bCount;
	Material	   *materials;
	Region		   *regions;
	Boundary	   *boundaries;

	//! консервативные переменные на текущем временном слое.
	double		   *ro;			 
	double		   *ru;			
	double		   *rv;			
	double		   *re;

	double			*tmpArr;
	int				*tmpArrInt;

	//! градиенты.
	Vector			*gradR;
	Vector			*gradP;
	Vector			*gradU;
	Vector			*gradV;	

	//! лимиты
	double			limitRmin;
	double			limitRmax;
	double			limitPmin;
	double			limitPmax;
	double			limitUmax;

	//! подъемна€ сила.
	double			Fx;
	double			Fy;
	
protected:
	const static int FLUX_GODUNOV		= 0;
	const static int FLUX_LAX			= 1;
};

#endif