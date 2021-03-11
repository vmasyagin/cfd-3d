#pragma once

#include "grid.h"

#define _DIST_2_ 

class MeshReader
{
public:
	static const int TYPE_WIAS_TET		= 1;
	static const int TYPE_SALOME		= 2;
	static const int TYPE_GMSH			= 3;

	virtual void read(Grid*) = 0;
	static MeshReader* create(int type, char* fileName);
	static int getType(char* name);
	//! расчёт квадрата расстояния между двумя точками
	inline double _calc_dist_2_(const Point& p1, const Point& p2)
	{
		return ((p2.x - p1.x)* (p2.x - p1.x) + (p2.y - p1.y)* (p2.y - p1.y) + (p2.z - p1.z)* (p2.z - p1.z));
	}
};

