#ifndef _GRID_H_
#define _GRID_H_

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <vector>
#include <map>
#include <set>
#include "global.h"
#include "bnd_cond.h"

const unsigned int CELL_FLAG_GOOD = 0x000000;
const unsigned int CELL_FLAG_BAD = 0x000001;
const unsigned int CELL_FLAG_LIM = 0x000002;



class Cell
{
public:
	Cell() : nodesInd(NULL), facesInd(NULL), nCount(0), fCount(0) {};
	~Cell();

	int				nCount;
	int				fCount;
	int*			nodesInd;
	int*			facesInd;
	int*			neigh;
	int				type;
	char			typeName[64];
	double			V;
	Point			c; // центр €чейки
	double			H; // характерный размер €чейки
	double			HX;
	double			HY;
	double			HZ;
	unsigned int	flag;

	friend class	Grid;
};

class Face
{
public:
	Face() : nodesInd(NULL), edgesInd(NULL), cCount(0), nCount(0) {};
	~Face();

	int*			nodesInd;      // узлы, образующие грань
	int*			edgesInd;
	int				nCount; // количество узлов грани
	int				c1;        // €чейка слева
	int				c2;        // €чейка справа, нормаль из с1 в с2
	int				L;         // €чейка слева (относительно вектора beta)
	int				R;         // €чейка справа (относительно вектора beta)
	Vector			n;         // нормаль к грани
	double			l;         // длина грани 
	double			cnl1;       // сумма длин проекций векторов из центра грани до центра смежной €чейки
	double			cnl2;       // сумма длин проекций векторов из центра грани до центра смежной €чейки
	int				cCount;    // количество точек на грани
	double			S;
	Point			c;         // центр грани
	int				type;      // тип грани (внутр., гранич.)
	char			typeName[64];
	CFDBoundary*	bnd = NULL;
	friend class	Grid;
public:
	static const int TYPE_INNER = 0;
	static const int TYPE_INLET = 1;
	static const int TYPE_OUTLET = 2;
	static const int TYPE_WALL = 3;
	static const int TYPE_DUAL = 4;

	static const int TYPE_NAMED = 0x100;
};

class Edge
{
public:
	Edge() : c(NULL), cCount(0) {};
	~Edge();

	int				n1;        // узел в начале
	int				n2;        // узел в конце
	int				c1;        // €чейка слева
	int				c2;        // €чейка справа, нормаль из с1 в с2
	Vector			n;         // нормаль к грани
	double			l;         // длина грани 
	double			cnl1;       // сумма длин проекций векторов из центра грани до центра смежной €чейки
	double			cnl2;       // сумма длин проекций векторов из центра грани до центра смежной €чейки
	int				cCount;    // количество точек на грани
	Point*			c;         // точки на грани
	int				type;      // тип грани (внутр., гранич.)
	char			typeName[64];
	CFDBoundary*	bnd = NULL;
	friend class	Grid;
public:
	static const int TYPE_INNER = 0;
	static const int TYPE_INLET = 1;
	static const int TYPE_OUTLET = 2;
	static const int TYPE_WALL = 3;
	static const int TYPE_DUAL = 4;

	static const int TYPE_NAMED = 0x100;
};

class Grid
{
public:
	Grid() : nodes(NULL), cells(NULL), faces(NULL), edges(NULL), dcells(NULL), dnodes(NULL), tetrahedrons(NULL), dfaces(NULL),
		nCount(0), cCount(0), fCount(0), eCount(0), dfCount(0) {};
	~Grid();

	void initFromFiles(char* fName);
	inline Point& getNode(int i) { return nodes[i]; };
	inline Cell& getCell(int i) { return cells[i]; };
	inline Face& getFace(int i) { return faces[i]; };
	//void replaceFaces(int if1, int if2);
	//virtual void reorderFaces();
	virtual void reorientNormals(int n, int p);
	int findEdge(int n1, int n2);
	int findFace(int n1, int n2, int n3);

	Point* nodes;
	Cell* cells;
	Face* faces;
	Edge* edges;
	Cell* dcells;
	Point* dnodes;
	Cell* tetrahedrons;
	Cell* pyramids;
	Face* dfaces;
	int    nCount;	//количество nodes.
	int    cCount;  //количество €чеек(cells).
	int    fCount;  //количество граней(faces).
	int    eCount;  //количество ребер(edges).
	int    tCount;	//количество tetrahedrons.
	int    pCount;	//количество pyramids.
	int    dnCount;	//количество dnodes.
	int    dfCount;  //количество граней двойственной сетки (dfaces).
	int    nCountEx;	//количество nodes с фиктивными €чйками.
	int    cCountEx;  //количество €чеек(cells) с фиктивными €чйками.
	int    fCountEx;  //количество граней(faces) с фиктивными €чйками.
	int    eCountEx;  //количество ребер(edges) с фиктивными €чйками.
	double* con11;

	std::vector<int>					recvCount;
	std::vector<int>					recvShift;
	std::vector< std::vector<int> >		sendInd;
	std::map<int, std::set<int> >		node_cells;
	std::map<int, std::set<int> >		dnode_tetrahedrons;
	std::map<int, std::set<int> >		dnode_pyramids;
	std::vector<int>					dface_cell;
	std::vector<int>					tetrahedron_cell;
	std::vector<int>					pyramid_cell;
	int* pyramid_to_cell;
	int* pyramid_to_node;

	void readMeshFiles();
	void saveMeshInfo();
};

inline double _max_(double a, double b)
{
	if (a > b) { return a; }
	else { return b; }
}

inline double _min_(double a, double b)
{
	if (a < b) { return a; }
	else { return b; }
}

inline double _max_(double a, double b, double c)
{
	return _max_(a, _max_(b, c));
}

inline double _min_(double a, double b, double c)
{
	return _min_(a, _min_(b, c));
}

inline double _sqr_(double a)
{
	return a * a;
}

#endif