#ifndef _BND_COND_H_
#define _BND_COND_H_

#include <cstdlib>
#include "global.h"
#include "tinyxml.h"

#include <vector>

class Face;
class Grid;
class Material;
struct CFDBoundary
{
	static CFDBoundary* create(TiXmlNode* bNode, Grid * g);
	virtual void run(int iFace, Param& pL, Param& pR) = 0;

	CFDBoundary() : par(NULL) {}
	~CFDBoundary() { if (par) delete[] par; par = NULL; }


	int			faceType;
	char		name[64];
	double *	par;
	int			parCount;
	Grid *		g;

    static const char* TYPE_INLET_SUB;
	static const char* TYPE_INLET_SUPER;
    static const char* TYPE_OUTLET_SUB;
	static const char* TYPE_OUTLET_SUPER;
    static const char* TYPE_PRESSURE;
    static const char* TYPE_WALL_SLIP;
    static const char* TYPE_WALL_NO_SLIP;
	static const char* TYPE_SYMMETRY;
	static const char* TYPE_MASSFLOW;
	static const char* TYPE_FREE_STREAM;
	static const char* TYPE_STAGNATION;

    static const int TYPE_ID_INLET_SUB;
	static const int TYPE_ID_INLET_SUPER;
    static const int TYPE_ID_OUTLET_SUB;
	static const int TYPE_ID_OUTLET_SUPER;
    static const int TYPE_ID_PRESSURE;
    static const int TYPE_ID_WALL_SLIP;
    static const int TYPE_ID_WALL_NO_SLIP;
	static const int TYPE_ID_SYMMETRY;
	static const int TYPE_ID_MASSFLOW;
	static const int TYPE_ID_FREE_STREAM;
	static const int TYPE_ID_STAGNATION;
};
typedef std::vector< CFDBoundary* > CFDBoundaries;

struct CFDBndInletSub : public CFDBoundary
{
	virtual void run(int iFace, Param& pL, Param& pR);
};

struct CFDBndInletSuper : public CFDBoundary
{
	virtual void run(int iFace, Param& pL, Param& pR);
};

struct CFDBndOutletSub : public CFDBoundary
{
	virtual void run(int iFace, Param& pL, Param& pR);
};

struct CFDBndOutletSuper : public CFDBoundary
{
	virtual void run(int iFace, Param& pL, Param& pR);
};

struct CFDBndPressure : public CFDBoundary
{
	virtual void run(int iFace, Param& pL, Param& pR);
};

struct CFDBndWallSlip : public CFDBoundary
{
	virtual void run(int iFace, Param& pL, Param& pR);
};

struct CFDBndWallNoSlip : public CFDBoundary
{
    virtual void run(int iFace, Param& pL, Param& pR);
};

struct CFDBndSymmetry : public CFDBoundary
{
	virtual void run(int iFace, Param& pL, Param& pR);
};

struct CFDBndMassFlow : public CFDBoundary
{
	virtual void run(int iFace, Param& pL, Param& pR);
};

struct CFDBndFreeStream : public CFDBoundary
{
	virtual void run(int iFace, Param& pL, Param& pR);
};

struct CFDBndStagnation : public CFDBoundary
{
	virtual void run(int iFace, Param& pL, Param& pR);
};

#endif
