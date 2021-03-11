#include "grid.h"
#include <map>
#include <set>
#include <algorithm>

typedef std::vector<int> idx_t;
typedef std::vector<idx_t> idx2_t;
typedef std::tuple<int, int, int> face_t;
typedef std::map<face_t, int> face_idx_t;
typedef std::pair<int, int> edge_t;
typedef std::map<edge_t, int> edge_idx_t;

edge_idx_t nodeToEdge;
face_idx_t nodeToFace;

int findEdgeByNodes(int n1, int n2)
{
	if (n1 > n2) {
		int tmp = n1;
		n1 = n2;
		n2 = tmp;
	}
	return nodeToEdge[edge_t(n1, n2)];
}

int findFaceByNodes(int n1, int n2, int n3)
{
	int _n1 = _min_(n1, n2, n3);
	int _n3 = _max_(n1, n2, n3);
	int _n2 = (n1 + n2 + n3) - _n1 - _n3;
	return nodeToFace[face_t(_n1, _n2, _n3)];
}

Cell::~Cell() 
{
	delete[] nodesInd;
	delete[] facesInd;
	delete[] neigh;
}

Face::~Face()
{
	delete[] nodesInd;
	delete[] edgesInd;
	delete[] bnd;
}

Edge::~Edge() 
{
	delete[] c;
	delete[] bnd;
}

Grid::~Grid() 
{
	delete[] nodes;
	delete[] cells;
	delete[] faces;
	delete[] edges;
	delete[] pyramids;
	delete[] tetrahedrons;
	delete[] dcells;
	delete[] dfaces;
	delete[] dnodes;
	delete[] con11;
}

int Grid::findEdge(int n1, int n2)
{
	for (int iEdge = 0; iEdge < eCount; iEdge++)
	{
		if ((edges[iEdge].n1 == n1 && edges[iEdge].n2 == n2) || (edges[iEdge].n1 == n2 && edges[iEdge].n2 == n1))
		{
			return iEdge;
		}
	}
	return -1;
}

int Grid::findFace(int n1, int n2, int n3)
{
	std::set<int> searchNodes{ n1, n2, n3 };
	for (int iFace = 0; iFace < fCountEx; iFace++)
	{
		std::set<int> faceNodes{ faces[iFace].nodesInd[0], faces[iFace].nodesInd[1], faces[iFace].nodesInd[2] };
		if (searchNodes == faceNodes)
		{
			return iFace;
			break;
		}
	}
	return -1;
}

//Grid::initFromFiles: загрузка сетки из файла fName.
/*
void Grid::initFromFiles(char* fName) 
{
	char str[50];
	FILE *fp;
	int tmp; 

	// читаем данные об ”«Ћј’
	sprintf(str, "%s.node", fName);
	fp = fopen(str, "r");
	if (!fp) {
		log("Can not open file '%s'\n", str);
		EXIT(1);
	}
	fscanf(fp, "%d %d %d %d", &nCount, &tmp, &tmp, &tmp);
	nodes = new Point[nCount];
	for (int i = 0; i < nCount; i++) 
	{
		fscanf(fp, "%d %lf %lf %d", &tmp, &(nodes[i].x), &(nodes[i].y), &tmp);
	}
	fclose(fp);

	// читаем данные о я„≈… ј’
	sprintf(str, "%s.ele", fName);
	fp = fopen(str, "r");
	if (!fp) {
		log("Can not open file '%s'\n", str);
		EXIT(1);
	}
	fscanf(fp, "%d %d %d", &cCount, &tmp, &tmp);
	cells = new Cell[cCount];
	for (int i = 0; i < cCount; i++) 
	{
		cells[i].nCount = 3;
		cells[i].nodesInd = new int[cells[i].nCount];
		fscanf(fp, "%d %d %d %d %d", &tmp, &(cells[i].nodesInd[0]), &(cells[i].nodesInd[1]), &(cells[i].nodesInd[2]), &(cells[i].type));
		cells[i].nodesInd[0]--;
		cells[i].nodesInd[1]--;
		cells[i].nodesInd[2]--;
		cells[i].c.x = (nodes[cells[i].nodesInd[0]].x+nodes[cells[i].nodesInd[1]].x+nodes[cells[i].nodesInd[2]].x)/3.0;
		cells[i].c.y = (nodes[cells[i].nodesInd[0]].y+nodes[cells[i].nodesInd[1]].y+nodes[cells[i].nodesInd[2]].y)/3.0;
		cells[i].HX = _max_( fabs(nodes[cells[i].nodesInd[0]].x-nodes[cells[i].nodesInd[1]].x), 
			                 fabs(nodes[cells[i].nodesInd[1]].x-nodes[cells[i].nodesInd[2]].x),
							 fabs(nodes[cells[i].nodesInd[0]].x-nodes[cells[i].nodesInd[2]].x) );
		cells[i].HY = _max_( fabs(nodes[cells[i].nodesInd[0]].y-nodes[cells[i].nodesInd[1]].y), 
			                 fabs(nodes[cells[i].nodesInd[1]].y-nodes[cells[i].nodesInd[2]].y),
							 fabs(nodes[cells[i].nodesInd[0]].y-nodes[cells[i].nodesInd[2]].y) );
		cells[i].eCount = 3;
		cells[i].edgesInd = new int[cells[i].eCount];
	}
	fclose(fp);

	// формируем данные о –≈Ѕ–ј’
	sprintf(str, "%s.neigh", fName);
	fp = fopen(str, "r");
	if (!fp) {
		log("Can not open file '%s'\n", str);
		EXIT(1);
	}
	fscanf(fp, "%d %d", &tmp, &tmp);
	int** neigh;
	neigh = new int*[cCount]; 
	for (int i = 0; i < cCount; i++) 
	{
		neigh[i] = new int[3];
		fscanf(fp, "%d %d %d %d", &tmp, &(neigh[i][0]), &(neigh[i][1]), &(neigh[i][2]));
		neigh[i][0]--;
		neigh[i][1]--;
		neigh[i][2]--;
		cells[i].neigh[0] = neigh[i][0];
		cells[i].neigh[1] = neigh[i][1];
		cells[i].neigh[2] = neigh[i][2];
	}
	fclose(fp);
	eCount = 0;
	for (int i = 0; i < cCount; i++) 
	{
		for (int j = 0; j < 3; j++) 
		{
			int p = neigh[i][j];
			if (p > -1) 
			{
				for (int k = 0; k < 3; k++) 
				{ // убираем у соседа номер этой €чейки, чтобы грань не повтор€лась
					if (neigh[p][k] == i) neigh[p][k] = -1;
				}
				eCount++;
			}
			if (p == -2) eCount++;
		}
	}
	edges = new Edge[eCount];

	int iEdge = 0;
	int * cfi = new int[cCount];
	for (int i = 0; i < cCount; i++) 
	{
		cfi[i] = 0;
	}
	// ::memset(cfi, 0, cCount*sizeof(int));
	for (int i = 0; i < cCount; i++) 
	{
		for (int j = 0; j < 3; j++) 
		{
			int p = neigh[i][j];
			if (p != -1) 
			{
				edges[iEdge].n1     = cells[i].nodesInd[(j+1)%3];
				edges[iEdge].n2     = cells[i].nodesInd[(j+2)%3];
				
				edges[iEdge].cCount = 3;
				edges[iEdge].c      = new Point[edges[iEdge].cCount];
				double _sqrt3 = 1.0/sqrt(3.0);
				// центр ребра
				edges[iEdge].c[0].x = (nodes[edges[iEdge].n1].x+nodes[edges[iEdge].n2].x)/2.0;
				edges[iEdge].c[0].y = (nodes[edges[iEdge].n1].y+nodes[edges[iEdge].n2].y)/2.0;
				// перва§ точка vаусса
				edges[iEdge].c[1].x = (nodes[edges[iEdge].n1].x+nodes[edges[iEdge].n2].x)/2.0-_sqrt3*(nodes[edges[iEdge].n2].x-nodes[edges[iEdge].n1].x)/2.0;
				edges[iEdge].c[1].y = (nodes[edges[iEdge].n1].y+nodes[edges[iEdge].n2].y)/2.0-_sqrt3*(nodes[edges[iEdge].n2].y-nodes[edges[iEdge].n1].y)/2.0;
				// втора§ точка vаусса
				edges[iEdge].c[2].x = (nodes[edges[iEdge].n1].x+nodes[edges[iEdge].n2].x)/2.0+_sqrt3*(nodes[edges[iEdge].n2].x-nodes[edges[iEdge].n1].x)/2.0;
				edges[iEdge].c[2].y = (nodes[edges[iEdge].n1].y+nodes[edges[iEdge].n2].y)/2.0+_sqrt3*(nodes[edges[iEdge].n2].y-nodes[edges[iEdge].n1].y)/2.0;

				edges[iEdge].n.x    = nodes[edges[iEdge].n2].y-nodes[edges[iEdge].n1].y;
				edges[iEdge].n.y    = nodes[edges[iEdge].n1].x-nodes[edges[iEdge].n2].x;
				edges[iEdge].l      = sqrt(edges[iEdge].n.x*edges[iEdge].n.x+edges[iEdge].n.y*edges[iEdge].n.y);
				edges[iEdge].n.x    /= edges[iEdge].l;
				edges[iEdge].n.y    /= edges[iEdge].l;
				edges[iEdge].c1     = i;
				cells[i].edgesInd[cfi[i]] = iEdge;
				cfi[i]++;
				edges[iEdge].cnl1 = fabs(edges[iEdge].n.x*(edges[iEdge].c[0].x-cells[edges[iEdge].c1].c.x)+edges[iEdge].n.y*(edges[iEdge].c[0].y-cells[edges[iEdge].c1].c.y) );

				if (p > -1) 
				{

					edges[iEdge].c2 = p;
					cells[p].edgesInd[cfi[p]] = iEdge;
					cfi[p]++;
					edges[iEdge].cnl2 = fabs(edges[iEdge].n.x*(cells[edges[iEdge].c2].c.x-edges[iEdge].c[0].x)+edges[iEdge].n.y*(cells[edges[iEdge].c2].c.y-edges[iEdge].c[0].y) );
					edges[iEdge].type = Edge::TYPE_INNER;
				}
				if (p == -2) 
				{
					edges[iEdge].c2 = -1;
					edges[iEdge].cnl2 = 0;
					edges[iEdge].type = -1;
				}
				iEdge++;
			}
		}
		
	}
	
	// чтение данных о граничных гран€х
	sprintf(str, "%s.poly", fName);
	fp = fopen(str, "r");
	if (!fp) {
		log("Can not open file '%s'\n", str);
		EXIT(1);
	}
	int bndCount;
	fscanf(fp, "%d %d %d %d", &tmp, &tmp, &tmp, &tmp);
	fscanf(fp, "%d %d", &bndCount, &tmp);
	for (int i = 0; i < bndCount; i++) 
	{
		int n, n1, n2, type;
		fscanf(fp, "%d %d %d %d", &n, &n1, &n2, &type);
		n1--;
		n2--;
		int iEdge = findEdge(n1, n2);
		if (iEdge >= 0) edges[iEdge].type = type;
	}
	fclose(fp);

	for (int i = 0; i < cCount; i++) 
	{
		double a = edges[cells[i].edgesInd[0]].l;
		double b = edges[cells[i].edgesInd[1]].l;
		double c = edges[cells[i].edgesInd[2]].l;
		double p = (a+b+c)/2.0;
		cells[i].S = sqrt(p*(p-a)*(p-b)*(p-c));
	}
	for (int i = 0; i < cCount; i++)		
	{		
		cells[i].flag = CELL_FLAG_GOOD;		
	}

	for (int i = 0; i < cCount; i++) 
	{
		cells[i].flag = 0;
	}
	
	for (int i = 0; i < cCount; i++) 
	{
		delete[] neigh[i];
	}
	delete[] neigh;
	delete[] cfi;


	for (int i = 0; i < eCount; i++) {

	}


}
*/
//typedef std::vector<int> idx_t;
//typedef std::vector<idx_t> idx2_t;
//typedef std::pair<int, int> edge_t;
//typedef std::map<edge_t, int> edge_idx_t;
//
//edge_idx_t nodeToEdge;
//int findEdgeByNodes(int n1, int n2)
//{
//	if (n1 > n2) {
//		int tmp = n1;
//		n1 = n2;
//		n2 = tmp;
//	}
//	return nodeToEdge[edge_t(n1, n2)];
//}


void Grid::readMeshFiles()
{
	int p = Parallel::procId;
	char fName[64];
	int tmp;
	sprintf(fName, "mesh/mesh.%04d.proc", p);
	FILE* fp = fopen(fName, "r");
	if (!fp) {
		log("Can not open file '%s'\n", fName);
		EXIT(1);
	}

	log("Reading mesh part structure:\n");

	// Nodes
	log("\t- nodes;\n");
	fscanf(fp, "%d %d", &nCount, &nCountEx);
	nodes = new Point[nCountEx];
	for (int i = 0; i < nCountEx; i++)
	{
		fscanf(fp, "%d %lf %lf %lf", &tmp, &(nodes[i].x), &(nodes[i].y), &(nodes[i].z));
	}

	// Cells
	log("\t- cells;\n");
	fscanf(fp, "%d %d", &cCount, &cCountEx);
	cells = new Cell[cCountEx];
	for (int i = 0; i < cCountEx; i++)
	{
		cells[i].nCount = 3;
		cells[i].nodesInd = new int[cells[i].nCount];
		fscanf(fp, "%d %d %d %d %d %s", &tmp, &(cells[i].nodesInd[0]), &(cells[i].nodesInd[1]), &(cells[i].nodesInd[2]), &(cells[i].nodesInd[3]), &(cells[i].typeName));
		cells[i].c.x = (nodes[cells[i].nodesInd[0]].x + nodes[cells[i].nodesInd[1]].x + nodes[cells[i].nodesInd[2]].x + nodes[cells[i].nodesInd[3]].x) / 4.0;
		cells[i].c.y = (nodes[cells[i].nodesInd[0]].y + nodes[cells[i].nodesInd[1]].y + nodes[cells[i].nodesInd[2]].y + nodes[cells[i].nodesInd[3]].y) / 4.0;
		cells[i].c.z = (nodes[cells[i].nodesInd[0]].z + nodes[cells[i].nodesInd[1]].z + nodes[cells[i].nodesInd[2]].z + nodes[cells[i].nodesInd[3]].z) / 4.0;
		cells[i].HX = _max_(_max_(nodes[cells[i].nodesInd[0]].x, nodes[cells[i].nodesInd[1]].x), _max_(nodes[cells[i].nodesInd[2]].x, nodes[cells[i].nodesInd[3]].x)) -
			_min_(_min_(nodes[cells[i].nodesInd[0]].x, nodes[cells[i].nodesInd[1]].x), _min_(nodes[cells[i].nodesInd[2]].x, nodes[cells[i].nodesInd[3]].x));
		cells[i].HY = _max_(_max_(nodes[cells[i].nodesInd[0]].y, nodes[cells[i].nodesInd[1]].y), _max_(nodes[cells[i].nodesInd[2]].y, nodes[cells[i].nodesInd[3]].y)) -
			_min_(_min_(nodes[cells[i].nodesInd[0]].y, nodes[cells[i].nodesInd[1]].y), _min_(nodes[cells[i].nodesInd[2]].y, nodes[cells[i].nodesInd[3]].y));
		cells[i].HZ = _max_(_max_(nodes[cells[i].nodesInd[0]].z, nodes[cells[i].nodesInd[1]].z), _max_(nodes[cells[i].nodesInd[2]].z, nodes[cells[i].nodesInd[3]].z)) -
			_min_(_min_(nodes[cells[i].nodesInd[0]].z, nodes[cells[i].nodesInd[1]].z), _min_(nodes[cells[i].nodesInd[2]].z, nodes[cells[i].nodesInd[3]].z));
		cells[i].fCount = 4;
		cells[i].facesInd = new int[cells[i].fCount];
	}
	// €чейки, относ€щиес€ к данному узлу
	//std::map<int, std::set<int> > node_cells;
	for (int i = 0; i < cCountEx; i++) {
		node_cells[cells[i].nodesInd[0]].insert(i);
		node_cells[cells[i].nodesInd[1]].insert(i);
		node_cells[cells[i].nodesInd[2]].insert(i);
		node_cells[cells[i].nodesInd[3]].insert(i);
	}

	// Faces
	log("\t- faces;\n");
	fscanf(fp, "%d %d", &fCount, &fCountEx);
	faces = new Face[fCountEx];
	for (int i = 0; i < fCountEx; i++)
	{
		faces[i].nCount = 3;
		faces[i].nodesInd = new int[faces[i].nCount];
		faces[i].edgesInd = new int[faces[i].nCount];
		fscanf(fp, "%d %d %d %d %d", &tmp, &(faces[i].nodesInd[0]), &(faces[i].nodesInd[1]), &(faces[i].nodesInd[2]), &(faces[i].type));
		if (faces[i].type != 0) {
			fscanf(fp, "%s", &(faces[i].typeName));
		}
		else {
			strcpy(faces[i].typeName, "");
		}
		int n1, n2, n3;
		n1 = _min_(faces[i].nodesInd[0], faces[i].nodesInd[1], faces[i].nodesInd[2]);
		n3 = _max_(faces[i].nodesInd[0], faces[i].nodesInd[1], faces[i].nodesInd[2]);
		n2 = (faces[i].nodesInd[0] + faces[i].nodesInd[1] + faces[i].nodesInd[2]) - n1 - n3;
		nodeToFace[face_t(n1, n2, n3)] = i;
	}

	idx2_t face_cells;
	face_cells.resize(fCountEx);
	for (int i = 0; i < cCountEx; i++) {
		for (int j = 0; j < 4; j++) {
			int n1 = cells[i].nodesInd[(j + 1) % 4];
			int n2 = cells[i].nodesInd[(j + 2) % 4];
			int n3 = cells[i].nodesInd[(j + 3) % 4];
			int iFace = findFaceByNodes(n1, n2, n3);
			face_cells[iFace].push_back(i);
			cells[i].facesInd[j] = iFace;
		}
	}

	for (int iFace = 0; iFace < fCountEx; iFace++) {
		idx_t& fc = face_cells[iFace];
		Face& f = faces[iFace];
		if (fc.size() == 2) {
			f.c1 = fc[0];
			f.c2 = fc[1];
		}
		else if (fc.size() == 1) {
			f.c1 = fc[0];
			f.c2 = -1;
		}
		else {
			log("ERROR: face #%d is not assigned to any cell. Sourse: %s; line: %d\n", iFace, __FILE__, __LINE__);
			EXIT(1);
		}
		f.cCount = 1;
		//f.c = new Point[f.cCount];
		double x1 = nodes[f.nodesInd[0]].x; double y1 = nodes[f.nodesInd[0]].y; double z1 = nodes[f.nodesInd[0]].z;
		double x2 = nodes[f.nodesInd[1]].x; double y2 = nodes[f.nodesInd[1]].y; double z2 = nodes[f.nodesInd[1]].z;
		double x3 = nodes[f.nodesInd[2]].x; double y3 = nodes[f.nodesInd[2]].y; double z3 = nodes[f.nodesInd[2]].z;
		// центр грани
		f.c.x = (x1 + x2 + x3) / 3.0;
		f.c.y = (y1 + y2 + y3) / 3.0;
		f.c.z = (z1 + z2 + z3) / 3.0;
		// нормаль и площадь грани
		double a1 = x1 - x3; double a2 = x2 - x3; double a3 = x3;
		double b1 = y1 - y3; double b2 = y2 - y3; double b3 = y3;
		double c1 = z1 - z3; double c2 = z2 - z3; double c3 = z3;

		f.n.x = b1 * c2 - c1 * b2;
		f.n.y = -(a1 * c2 - c1 * a2);
		f.n.z = a1 * b2 - b1 * a2;
		f.n.x /= sqrt((b1 * c2 - c1 * b2) * (b1 * c2 - c1 * b2) + (a1 * c2 - c1 * a2) * (a1 * c2 - c1 * a2) + (a1 * b2 - b1 * a2) * (a1 * b2 - b1 * a2));
		f.n.y /= sqrt((b1 * c2 - c1 * b2) * (b1 * c2 - c1 * b2) + (a1 * c2 - c1 * a2) * (a1 * c2 - c1 * a2) + (a1 * b2 - b1 * a2) * (a1 * b2 - b1 * a2));
		f.n.z /= sqrt((b1 * c2 - c1 * b2) * (b1 * c2 - c1 * b2) + (a1 * c2 - c1 * a2) * (a1 * c2 - c1 * a2) + (a1 * b2 - b1 * a2) * (a1 * b2 - b1 * a2));

		// коррекци€ направлений нормалей
		Vector vc;

		vc.x = cells[f.c1].c.x - f.c.x;
		vc.y = cells[f.c1].c.y - f.c.y;
		vc.z = cells[f.c1].c.z - f.c.z;
		if (scalar_prod(vc, f.n) > 0) {
			f.n.x *= -1.0;
			f.n.y *= -1.0;
			f.n.z *= -1.0;
		}
	}

	// формируем данные о –≈Ѕ–ј’
	std::map<std::pair<int, int>, std::set<int>> pair_faces;
	for (int i = 0; i < cCountEx; i++)
	{
		Cell& cell = cells[i];
		for (int j = 0; j < 4; j++)
		{
			Face& face = faces[cell.facesInd[j]];
			for (int k = 0; k < 3; k++) {
				int n1 = _min_(face.nodesInd[k], face.nodesInd[(k + 1) % 3]);
				int n2 = _max_(face.nodesInd[k], face.nodesInd[(k + 1) % 3]);
				pair_faces[std::pair<int, int>(n1, n2)].insert(cell.facesInd[j]);
			}
		}
	}
	eCountEx = pair_faces.size();
	edges = new Edge[eCountEx];

	int iEdge = 0;
	for (std::map<std::pair<int, int>, std::set<int>>::iterator it = pair_faces.begin(); it != pair_faces.end(); ++it) {
		edges[iEdge].n1 = it->first.first;
		edges[iEdge].n2 = it->first.second;
		for (std::set<int>::iterator it_set = it->second.begin(); it_set != it->second.end(); ++it_set) {
			for (int j = 0; j < 3; j++) {
				int n1 = _min_(faces[*it_set].nodesInd[j], faces[*it_set].nodesInd[(j + 1) % 3]);
				int n2 = _max_(faces[*it_set].nodesInd[j], faces[*it_set].nodesInd[(j + 1) % 3]);
				if (n1 == edges[iEdge].n1 && n2 == edges[iEdge].n2) {
					faces[*it_set].edgesInd[(j + 2) % 3] = iEdge;
				}
			}
		}
		iEdge++;
	}

	// dedges
	log("\t- dfaces;\n");

	dfCount = 0;
	dnCount = 0;

	for (int iCell = 0; iCell < cCountEx; iCell++) {
		//dnCount++; // центр тетраэдра
		Cell& cell = cells[iCell];
		// цикл по гран€м
		for (int i = 0; i < 4; i++) {
			Face& face = faces[cell.facesInd[i]];
			if (face.nodesInd[0] < nCount || face.nodesInd[1] < nCount || face.nodesInd[2] < nCount) {
				//dnCount++; // центр грани
				// цикл по ребрам
				for (int j = 0; j < 3; j++) {
					Edge& edge = edges[face.edgesInd[j]];
					if (edge.n1 < nCount || edge.n2 < nCount) {
						//dnCount++; // середина ребра
						if (face.type != Face::TYPE_INNER) {
							if (edge.n1 < nCount && edge.n2 < nCount)
								dfCount += 3;
							else
								dfCount += 2;
						}
						else if (face.type == Face::TYPE_INNER) {
							dfCount++;
						}
					}
				}
			}
		}
	}
	for (int i = 0; i < fCountEx; i++) {
		Face& face = faces[i];
		if (face.nodesInd[0] < nCount || face.nodesInd[1] < nCount || face.nodesInd[2] < nCount) {
			dnCount++; // центр грани
		}
	}
	for (int i = 0; i < eCountEx; i++) {
		Edge& edge = edges[i];
		if (edge.n1 < nCount || edge.n2 < nCount) dnCount++; // центр ребра
	}
	dnCount += nCountEx + cCountEx;
	dfaces = new Face[dfCount];
	dnodes = new Point[dnCount];

	pCount = 0;
	for (int i = 0; i < nCount; i++) {
		for (std::set<int>::iterator it = node_cells[i].begin(); it != node_cells[i].end(); ++it) {
			pCount += 3;
		}
	}

	memcpy(dnodes, nodes, nCountEx * sizeof(Point));
	std::map<int, int> face_dnode;
	int iNode = nCountEx;
	// точки двойственной сетки на гран€х
	for (int i = 0; i < fCountEx; i++) {
		Face& face = faces[i];
		if (face.nodesInd[0] < nCount || face.nodesInd[1] < nCount || face.nodesInd[2] < nCount) {
			dnodes[iNode] = face.c;
			face_dnode.insert(std::pair<int, int>(i, iNode));
			iNode++;
		}
	}
	// точки двойственной сетки на ребрах
	std::map<int, int> edge_dnode;
	for (int i = 0; i < eCountEx; i++) {
		Edge& edge = edges[i];
		if (edge.n1 < nCount || edge.n2 < nCount) {
			dnodes[iNode] = nodes[edge.n1] + nodes[edge.n2];
			dnodes[iNode] *= 0.5;
			edge_dnode.insert(std::pair<int, int>(i, iNode));
			iNode++;
		}
	}
	con11 = new double[fCountEx];
	for (int i = 0; i < fCountEx; i++) con11[i] = 1.0;

	int iPyramid = 0;
	int iFace = 0;

	//tetrahedrons = new Cell[tCount];
	pyramids = new Cell[pCount];
	pyramid_to_node = new int[pCount];
	pyramid_to_cell = new int[pCount];
	for (int iCell = 0; iCell < cCountEx; iCell++) {
		Cell& cell = cells[iCell];
		dnodes[iNode] = cell.c;
		// цикл по гран€м
		for (int i = 0; i < 4; i++) {
			Face& face = faces[cell.facesInd[i]];
			// цикл по вершинам грани
			for (int j = 0; j < 3; j++) {
				if (face.nodesInd[j] < nCount) {
					/*tetrahedrons[iCell].nCount = 4;
					tetrahedrons[iCell].nodesInd = new int[tetrahedrons[i].nCount];
					tetrahedrons[iCell].nodesInd[0] = cell.nodesInd[j];
					tetrahedrons[iCell].nodesInd[1] = edge_dnode[findEdge(cell.nodesInd[j], cell.nodesInd[(j + 1) % 3])];
					tetrahedrons[iCell].nodesInd[2] = iNode;
					tetrahedrons[iCell].nodesInd[3] = edge_dnode[findEdge(cell.nodesInd[j], cell.nodesInd[(j + 2) % 3])];
					dnode_tetrahedrons[cell.nodesInd[j]].insert(iCell);
					tetrahedron_cell.push_back(i);*/
					pyramids[iPyramid].nCount = 5;
					pyramids[iPyramid].nodesInd = new int[pyramids[iPyramid].nCount];
					pyramids[iPyramid].nodesInd[0] = face.nodesInd[j];
					pyramids[iPyramid].nodesInd[1] = edge_dnode[face.edgesInd[(j + 1) % 3]];
					pyramids[iPyramid].nodesInd[2] = face_dnode[cell.facesInd[i]];
					pyramids[iPyramid].nodesInd[3] = edge_dnode[face.edgesInd[(j + 2) % 3]];
					pyramids[iPyramid].nodesInd[4] = iNode;
					dnode_pyramids[face.nodesInd[j]].insert(iPyramid);
					//pyramid_cell.push_back(iCell);
					pyramid_to_cell[iPyramid] = iCell;
					pyramid_to_node[iPyramid] = face.nodesInd[j];
					iPyramid++;
				}
			}
			// dfaces
			// цикл по ребрам грани
			for (int j = 0; j < 3; j++) {
				Edge& edge = edges[face.edgesInd[j]];
				if (edge.n1 < nCount || edge.n2 < nCount) {
					if (face.type == Face::TYPE_INNER) {
						dfaces[iFace].nCount = 3;
						dfaces[iFace].nodesInd = new int[dfaces[iFace].nCount];
						dfaces[iFace].nodesInd[0] = edge_dnode[face.edgesInd[j]];
						dfaces[iFace].nodesInd[1] = iNode;
						dfaces[iFace].nodesInd[2] = face_dnode[cell.facesInd[i]];

						dfaces[iFace].cCount = 1;
						//dfaces[iFace].c = new Point[dfaces[iFace].cCount];
						double x1 = dnodes[dfaces[iFace].nodesInd[0]].x; double y1 = dnodes[dfaces[iFace].nodesInd[0]].y; double z1 = dnodes[dfaces[iFace].nodesInd[0]].z;
						double x2 = dnodes[dfaces[iFace].nodesInd[1]].x; double y2 = dnodes[dfaces[iFace].nodesInd[1]].y; double z2 = dnodes[dfaces[iFace].nodesInd[1]].z;
						double x3 = dnodes[dfaces[iFace].nodesInd[2]].x; double y3 = dnodes[dfaces[iFace].nodesInd[2]].y; double z3 = dnodes[dfaces[iFace].nodesInd[2]].z;
						// центр грани
						dfaces[iFace].c.x = (x1 + x2 + x3) / 3.0;
						dfaces[iFace].c.y = (y1 + y2 + y3) / 3.0;
						dfaces[iFace].c.z = (z1 + z2 + z3) / 3.0;
						// нормаль и площадь грани
						double a1 = x1 - x3; double a2 = x2 - x3; double a3 = x3;
						double b1 = y1 - y3; double b2 = y2 - y3; double b3 = y3;
						double c1 = z1 - z3; double c2 = z2 - z3; double c3 = z3;

						dfaces[iFace].n.x = b1 * c2 - c1 * b2;
						dfaces[iFace].n.y = -(a1 * c2 - c1 * a2);
						dfaces[iFace].n.z = a1 * b2 - b1 * a2;
						dfaces[iFace].n.x /= sqrt((b1 * c2 - c1 * b2) * (b1 * c2 - c1 * b2) + (a1 * c2 - c1 * a2) * (a1 * c2 - c1 * a2) + (a1 * b2 - b1 * a2) * (a1 * b2 - b1 * a2));
						dfaces[iFace].n.y /= sqrt((b1 * c2 - c1 * b2) * (b1 * c2 - c1 * b2) + (a1 * c2 - c1 * a2) * (a1 * c2 - c1 * a2) + (a1 * b2 - b1 * a2) * (a1 * b2 - b1 * a2));
						dfaces[iFace].n.z /= sqrt((b1 * c2 - c1 * b2) * (b1 * c2 - c1 * b2) + (a1 * c2 - c1 * a2) * (a1 * c2 - c1 * a2) + (a1 * b2 - b1 * a2) * (a1 * b2 - b1 * a2));

						dfaces[iFace].c1 = edge.n1;
						dfaces[iFace].c2 = edge.n2;
						// коррекци€ направлений нормалей
						Vector vc;

						vc.x = nodes[dfaces[iFace].c1].x - dfaces[iFace].c.x;
						vc.y = nodes[dfaces[iFace].c1].y - dfaces[iFace].c.y;
						vc.z = nodes[dfaces[iFace].c1].z - dfaces[iFace].c.z;
						if (scalar_prod(vc, dfaces[iFace].n) > 0) {
							dfaces[iFace].n *= -1.0;
						}

						if (dfaces[iFace].c1 < nCount && dfaces[iFace].c2 < nCount) {
							dfaces[iFace].type = Face::TYPE_INNER;
						}
						else {
							dfaces[iFace].type = Face::TYPE_DUAL;
						}
						strcpy(dfaces[iFace].typeName, "");
						dface_cell.push_back(iCell);
						iFace++;
					}
					else {
						dfaces[iFace].nCount = 3;
						dfaces[iFace].nodesInd = new int[dfaces[iFace].nCount];
						dfaces[iFace].nodesInd[0] = edge_dnode[face.edgesInd[j]];
						dfaces[iFace].nodesInd[1] = iNode;
						dfaces[iFace].nodesInd[2] = face_dnode[cell.facesInd[i]];

						dfaces[iFace].cCount = 1;
						//dfaces[iFace].c = new Point[dfaces[iFace].cCount];
						double x1 = dnodes[dfaces[iFace].nodesInd[0]].x; double y1 = dnodes[dfaces[iFace].nodesInd[0]].y; double z1 = dnodes[dfaces[iFace].nodesInd[0]].z;
						double x2 = dnodes[dfaces[iFace].nodesInd[1]].x; double y2 = dnodes[dfaces[iFace].nodesInd[1]].y; double z2 = dnodes[dfaces[iFace].nodesInd[1]].z;
						double x3 = dnodes[dfaces[iFace].nodesInd[2]].x; double y3 = dnodes[dfaces[iFace].nodesInd[2]].y; double z3 = dnodes[dfaces[iFace].nodesInd[2]].z;
						// центр грани
						dfaces[iFace].c.x = (x1 + x2 + x3) / 3.0;
						dfaces[iFace].c.y = (y1 + y2 + y3) / 3.0;
						dfaces[iFace].c.z = (z1 + z2 + z3) / 3.0;
						// нормаль и площадь грани
						double a1 = x1 - x3; double a2 = x2 - x3; double a3 = x3;
						double b1 = y1 - y3; double b2 = y2 - y3; double b3 = y3;
						double c1 = z1 - z3; double c2 = z2 - z3; double c3 = z3;

						dfaces[iFace].n.x = b1 * c2 - c1 * b2;
						dfaces[iFace].n.y = -(a1 * c2 - c1 * a2);
						dfaces[iFace].n.z = a1 * b2 - b1 * a2;
						dfaces[iFace].n.x /= sqrt((b1 * c2 - c1 * b2) * (b1 * c2 - c1 * b2) + (a1 * c2 - c1 * a2) * (a1 * c2 - c1 * a2) + (a1 * b2 - b1 * a2) * (a1 * b2 - b1 * a2));
						dfaces[iFace].n.y /= sqrt((b1 * c2 - c1 * b2) * (b1 * c2 - c1 * b2) + (a1 * c2 - c1 * a2) * (a1 * c2 - c1 * a2) + (a1 * b2 - b1 * a2) * (a1 * b2 - b1 * a2));
						dfaces[iFace].n.z /= sqrt((b1 * c2 - c1 * b2) * (b1 * c2 - c1 * b2) + (a1 * c2 - c1 * a2) * (a1 * c2 - c1 * a2) + (a1 * b2 - b1 * a2) * (a1 * b2 - b1 * a2));

						dfaces[iFace].c1 = edge.n1;
						dfaces[iFace].c2 = edge.n2;
						// коррекци€ направлений нормалей
						Vector vc;

						vc.x = nodes[dfaces[iFace].c1].x - dfaces[iFace].c.x;
						vc.y = nodes[dfaces[iFace].c1].y - dfaces[iFace].c.y;
						vc.z = nodes[dfaces[iFace].c1].z - dfaces[iFace].c.z;
						if (scalar_prod(vc, dfaces[iFace].n) > 0) {
							dfaces[iFace].n *= -1.0;
						}

						if (dfaces[iFace].c1 < nCount && dfaces[iFace].c2 < nCount) {
							dfaces[iFace].type = Face::TYPE_INNER;
						}
						else {
							dfaces[iFace].type = Face::TYPE_DUAL;
						}
						strcpy(dfaces[iFace].typeName, "");
						dface_cell.push_back(iCell);
						iFace++;

						if (edge.n1 < nCount) {
							// гранична€ грань #1
							dfaces[iFace].nCount = 3;
							dfaces[iFace].nodesInd = new int[dfaces[iFace].nCount];
							dfaces[iFace].nodesInd[0] = edge.n1;
							dfaces[iFace].nodesInd[1] = edge_dnode[face.edgesInd[j]];
							dfaces[iFace].nodesInd[2] = face_dnode[cell.facesInd[i]];

							dfaces[iFace].cCount = 1;
							//dfaces[iFace].c = new Point[dfaces[iFace].cCount];
							double x1 = dnodes[dfaces[iFace].nodesInd[0]].x; double y1 = dnodes[dfaces[iFace].nodesInd[0]].y; double z1 = dnodes[dfaces[iFace].nodesInd[0]].z;
							double x2 = dnodes[dfaces[iFace].nodesInd[1]].x; double y2 = dnodes[dfaces[iFace].nodesInd[1]].y; double z2 = dnodes[dfaces[iFace].nodesInd[1]].z;
							double x3 = dnodes[dfaces[iFace].nodesInd[2]].x; double y3 = dnodes[dfaces[iFace].nodesInd[2]].y; double z3 = dnodes[dfaces[iFace].nodesInd[2]].z;
							// центр грани
							dfaces[iFace].c.x = (x1 + x2 + x3) / 3.0;
							dfaces[iFace].c.y = (y1 + y2 + y3) / 3.0;
							dfaces[iFace].c.z = (z1 + z2 + z3) / 3.0;
							// нормаль и площадь грани
							double a1 = x1 - x3; double a2 = x2 - x3; double a3 = x3;
							double b1 = y1 - y3; double b2 = y2 - y3; double b3 = y3;
							double c1 = z1 - z3; double c2 = z2 - z3; double c3 = z3;

							dfaces[iFace].n = face.n;
							/*dfaces[iFace].n.x = b1*c2 - c1*b2;
							dfaces[iFace].n.y = -(a1*c2 - c1*a2);
							dfaces[iFace].n.z = a1*b2 - b1*a2;
							dfaces[iFace].n.x /= sqrt((b1*c2 - c1*b2)*(b1*c2 - c1*b2) + (a1*c2 - c1*a2)*(a1*c2 - c1*a2) + (a1*b2 - b1*a2)*(a1*b2 - b1*a2));
							dfaces[iFace].n.y /= sqrt((b1*c2 - c1*b2)*(b1*c2 - c1*b2) + (a1*c2 - c1*a2)*(a1*c2 - c1*a2) + (a1*b2 - b1*a2)*(a1*b2 - b1*a2));
							dfaces[iFace].n.z /= sqrt((b1*c2 - c1*b2)*(b1*c2 - c1*b2) + (a1*c2 - c1*a2)*(a1*c2 - c1*a2) + (a1*b2 - b1*a2)*(a1*b2 - b1*a2));*/

							dfaces[iFace].c1 = edge.n1;
							dfaces[iFace].c2 = -1;
							//// коррекци€ направлений нормалей
							//Vector vc;

							//vc.x = nodes[dfaces[iFace].c1].x - dfaces[iFace].c[0].x;
							//vc.y = nodes[dfaces[iFace].c1].y - dfaces[iFace].c[0].y;
							//vc.z = nodes[dfaces[iFace].c1].z - dfaces[iFace].c[0].z;
							//if (scalar_prod(vc, dfaces[iFace].n) > 0) {
							//	dfaces[iFace].n.x *= -1.0;
							//	dfaces[iFace].n.y *= -1.0;
							//	dfaces[iFace].n.z *= -1.0;
							//}
							dfaces[iFace].type = face.type;
							strcpy(dfaces[iFace].typeName, face.typeName);
							dface_cell.push_back(iCell);
							iFace++;
						}

						if (edge.n2 < nCount) {
							// граничное ребро #2
							dfaces[iFace].nCount = 3;
							dfaces[iFace].nodesInd = new int[dfaces[iFace].nCount];
							dfaces[iFace].nodesInd[0] = edge.n2;
							dfaces[iFace].nodesInd[1] = edge_dnode[face.edgesInd[j]];
							dfaces[iFace].nodesInd[2] = face_dnode[cell.facesInd[i]];

							dfaces[iFace].cCount = 1;
							//dfaces[iFace].c = new Point[dfaces[iFace].cCount];
							double x1 = dnodes[dfaces[iFace].nodesInd[0]].x; double y1 = dnodes[dfaces[iFace].nodesInd[0]].y; double z1 = dnodes[dfaces[iFace].nodesInd[0]].z;
							double x2 = dnodes[dfaces[iFace].nodesInd[1]].x; double y2 = dnodes[dfaces[iFace].nodesInd[1]].y; double z2 = dnodes[dfaces[iFace].nodesInd[1]].z;
							double x3 = dnodes[dfaces[iFace].nodesInd[2]].x; double y3 = dnodes[dfaces[iFace].nodesInd[2]].y; double z3 = dnodes[dfaces[iFace].nodesInd[2]].z;
							// центр грани
							dfaces[iFace].c.x = (x1 + x2 + x3) / 3.0;
							dfaces[iFace].c.y = (y1 + y2 + y3) / 3.0;
							dfaces[iFace].c.z = (z1 + z2 + z3) / 3.0;
							// нормаль и площадь грани
							double a1 = x1 - x3; double a2 = x2 - x3; double a3 = x3;
							double b1 = y1 - y3; double b2 = y2 - y3; double b3 = y3;
							double c1 = z1 - z3; double c2 = z2 - z3; double c3 = z3;

							dfaces[iFace].n = face.n;
							/*dfaces[iFace].n.x = b1*c2 - c1*b2;
							dfaces[iFace].n.y = -(a1*c2 - c1*a2);
							dfaces[iFace].n.z = a1*b2 - b1*a2;
							dfaces[iFace].n.x /= sqrt((b1*c2 - c1*b2)*(b1*c2 - c1*b2) + (a1*c2 - c1*a2)*(a1*c2 - c1*a2) + (a1*b2 - b1*a2)*(a1*b2 - b1*a2));
							dfaces[iFace].n.y /= sqrt((b1*c2 - c1*b2)*(b1*c2 - c1*b2) + (a1*c2 - c1*a2)*(a1*c2 - c1*a2) + (a1*b2 - b1*a2)*(a1*b2 - b1*a2));
							dfaces[iFace].n.z /= sqrt((b1*c2 - c1*b2)*(b1*c2 - c1*b2) + (a1*c2 - c1*a2)*(a1*c2 - c1*a2) + (a1*b2 - b1*a2)*(a1*b2 - b1*a2));*/

							dfaces[iFace].c1 = edge.n2;
							dfaces[iFace].c2 = -1;
							// коррекци€ направлений нормалей
							/*Vector vc;

							vc.x = nodes[dfaces[iFace].c1].x - dfaces[iFace].c[0].x;
							vc.y = nodes[dfaces[iFace].c1].y - dfaces[iFace].c[0].y;
							vc.z = nodes[dfaces[iFace].c1].z - dfaces[iFace].c[0].z;
							if (scalar_prod(vc, dfaces[iFace].n) > 0) {
								dfaces[iFace].n.x *= -1.0;
								dfaces[iFace].n.y *= -1.0;
								dfaces[iFace].n.z *= -1.0;
							}*/
							dfaces[iFace].type = face.type;
							strcpy(dfaces[iFace].typeName, face.typeName);
							dface_cell.push_back(iCell);
							iFace++;
						}
					}
				}
			}
		}
		iNode++;
	}
	/*pyramid_to_cell = new double[pCount];

	for (int i = 0; i < pyramid_cell.size(); i++) {
		pyramid_to_cell[i] = pyramid_cell[i];
	}*/
	// мен€ем ориентацию граничных граней двойственной сетки
	for (int i = 0; i < dfCount; i++) {
		if (dfaces[i].type == Face::TYPE_DUAL) {
			if (dfaces[i].c1 > nCount - 1) {
				int tmp = dfaces[i].c1;
				dfaces[i].c1 = dfaces[i].c2;
				dfaces[i].c2 = tmp;
				dfaces[i].n *= -1.0;
			}
		}
	}

	// dcells
	log("\t- dcells;\n");
	dcells = new Cell[nCount];
	for (int i = 0; i < nCount; i++) {
		dcells[i].c = nodes[i];
		dcells[i].V = 0.0;
		double x_max = 0.0, y_max = 0.0, z_max = 0.0, x_min = DBL_MAX, y_min = DBL_MAX, z_min = DBL_MAX;
		for (std::set<int>::iterator it = dnode_pyramids[i].begin(); it != dnode_pyramids[i].end(); ++it) {
			for (int j = 0; j < 5; j++) {
				if (dnodes[pyramids[*it].nodesInd[j]].x > x_max) x_max = dnodes[pyramids[*it].nodesInd[j]].x;
				if (dnodes[pyramids[*it].nodesInd[j]].x < x_min) x_min = dnodes[pyramids[*it].nodesInd[j]].x;
				if (dnodes[pyramids[*it].nodesInd[j]].y > y_max) y_max = dnodes[pyramids[*it].nodesInd[j]].y;
				if (dnodes[pyramids[*it].nodesInd[j]].y < y_min) y_min = dnodes[pyramids[*it].nodesInd[j]].y;
				if (dnodes[pyramids[*it].nodesInd[j]].z > z_max) z_max = dnodes[pyramids[*it].nodesInd[j]].z;
				if (dnodes[pyramids[*it].nodesInd[j]].z < z_min) z_min = dnodes[pyramids[*it].nodesInd[j]].z;
			}
			double x1 = dnodes[pyramids[*it].nodesInd[0]].x; double y1 = dnodes[pyramids[*it].nodesInd[0]].y; double z1 = dnodes[pyramids[*it].nodesInd[0]].z;
			double x2 = dnodes[pyramids[*it].nodesInd[1]].x; double y2 = dnodes[pyramids[*it].nodesInd[1]].y; double z2 = dnodes[pyramids[*it].nodesInd[1]].z;
			double x3 = dnodes[pyramids[*it].nodesInd[2]].x; double y3 = dnodes[pyramids[*it].nodesInd[2]].y; double z3 = dnodes[pyramids[*it].nodesInd[2]].z;
			double x4 = dnodes[pyramids[*it].nodesInd[3]].x; double y4 = dnodes[pyramids[*it].nodesInd[3]].y; double z4 = dnodes[pyramids[*it].nodesInd[3]].z;
			double x5 = dnodes[pyramids[*it].nodesInd[4]].x; double y5 = dnodes[pyramids[*it].nodesInd[4]].y; double z5 = dnodes[pyramids[*it].nodesInd[4]].z;
			dcells[i].V += fabs((x2 - x1) * ((y3 - y1) * (z5 - z1) - (z3 - z1) * (y5 - y1)) - (y2 - y1) * ((x3 - x1) * (z5 - z1) - (z3 - z1) * (x5 - x1)) + (z2 - z1) * ((x3 - x1) * (y5 - y1) - (y3 - y1) * (x5 - x1))) / 6.0;
			dcells[i].V += fabs((x3 - x1) * ((y4 - y1) * (z5 - z1) - (z4 - z1) * (y5 - y1)) - (y3 - y1) * ((x4 - x1) * (z5 - z1) - (z4 - z1) * (x5 - x1)) + (z3 - z1) * ((x4 - x1) * (y5 - y1) - (y4 - y1) * (x5 - x1))) / 6.0;
		}
		dcells[i].HX = x_max - x_min;
		dcells[i].HY = y_max - y_min;
		dcells[i].HZ = z_max - z_min;
	}

	for (int i = 0; i < cCountEx; i++) {
		double x1 = nodes[cells[i].nodesInd[0]].x; double y1 = nodes[cells[i].nodesInd[0]].y; double z1 = nodes[cells[i].nodesInd[0]].z;
		double x2 = nodes[cells[i].nodesInd[1]].x; double y2 = nodes[cells[i].nodesInd[1]].y; double z2 = nodes[cells[i].nodesInd[1]].z;
		double x3 = nodes[cells[i].nodesInd[2]].x; double y3 = nodes[cells[i].nodesInd[2]].y; double z3 = nodes[cells[i].nodesInd[2]].z;
		double x4 = nodes[cells[i].nodesInd[3]].x; double y4 = nodes[cells[i].nodesInd[3]].y; double z4 = nodes[cells[i].nodesInd[3]].z;
		cells[i].V = fabs((x2 - x1) * ((y3 - y1) * (z4 - z1) - (z3 - z1) * (y4 - y1)) - (y2 - y1) * ((x3 - x1) * (z4 - z1) - (z3 - z1) * (x4 - x1)) + (z2 - z1) * ((x3 - x1) * (y4 - y1) - (y3 - y1) * (x4 - x1))) / 6.0;
	}

	recvCount.clear();
	int sh = cCount;
	for (int i = 0; i < Parallel::procCount; i++) {
		fscanf(fp, "%d", &tmp);
		recvCount.push_back(tmp);
		recvShift.push_back(sh);
		sh += tmp;
	}

	sendInd.clear();
	for (int i = 0; i < Parallel::procCount; i++) {
		//int np;
		int nFirstLevel;
		int nSecondLevel;
		fscanf(fp, "%d %d %d", &tmp, &nFirstLevel, &nSecondLevel);
		std::vector<int> ind;
		ind.clear();
		for (int j = 0; j < nFirstLevel + nSecondLevel; j++) {
			fscanf(fp, "%d", &tmp);
			ind.push_back(tmp);
		}
		sendInd.push_back(ind);
		if (i > p && nFirstLevel > 0) reorientNormals(nFirstLevel, i);
	}
	fclose(fp);

	log("  complete...\n");


}

void Grid::reorientNormals(int n, int p) {
	int* adjFaces = new int[n];
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < 4; j++) {
			if (faces[cells[sendInd[p][i]].facesInd[j]].c2 >= cCount)
				adjFaces[i] = cells[sendInd[p][i]].facesInd[j];
		}
	}

	for (int i = 0; i < n; i++) {
		faces[adjFaces[i]].n *= -1.0;
		int tmp = faces[adjFaces[i]].c1;
		faces[adjFaces[i]].c1 = faces[adjFaces[i]].c2;
		faces[adjFaces[i]].c2 = tmp;
	}
	delete[]adjFaces;
}

void Grid::saveMeshInfo() {
	char name[64];
	FILE * fp;

	sprintf(name, "mesh.%04d.nodes.info", Parallel::procId);
	fp = fopen(name, "w");
	fprintf(fp, "COUNT: %6d    COUNT_EX: %6d\n", nCount, nCountEx);
	for (int i = 0; i < nCount; i++) {
		Point & p = nodes[i];
		fprintf(fp, "============================================================\n");
		fprintf(fp, "#%06d\n", i);
		fprintf(fp, "POINT: %25.15e %25.15e\n", p.x, p.y);
		fprintf(fp, "============================================================\n\n\n");
	}
	fprintf(fp, "\n********************* Extended *********************\n", cCount, cCountEx);
	for (int i = nCount; i < nCountEx; i++) {
		Point & p = nodes[i];
		fprintf(fp, "============================================================\n");
		fprintf(fp, "#%06d\n", i);
		fprintf(fp, "POINT: %25.15e %25.15e\n", p.x, p.y);
		fprintf(fp, "============================================================\n\n\n");
	}
	fclose(fp);

	sprintf(name, "mesh.%04d.cells.info", Parallel::procId);
	fp = fopen(name, "w");
	fprintf(fp, "COUNT: %6d    COUNT_EX: %6d\n", cCount, cCountEx);
	for (int i = 0; i < cCount; i++) {
		Cell & c = cells[i];
		fprintf(fp, "============================================================\n");
		fprintf(fp, "#%06d\n", i);
		fprintf(fp, "NODES: %6d %6d %6d %6d\n", c.nodesInd[0], c.nodesInd[1], c.nodesInd[2], c.nodesInd[3]);
		fprintf(fp, "EDGES: %6d %6d %6d %6d\n", c.facesInd[0], c.facesInd[1], c.facesInd[2], c.facesInd[3]);
		fprintf(fp, "CENTER: %25.15e %25.15e %25.15e\n", c.c.x, c.c.y, c.c.z);
		fprintf(fp, "HX: %25.15e   HY: %25.15e   HZ: %25.15e\n", c.HX, c.HY, c.HZ);
		fprintf(fp, "VOLUME: %25.15e\n", c.V);
		fprintf(fp, "TYPE: %s\n", c.typeName);
		fprintf(fp, "============================================================\n\n\n");
	}
	fprintf(fp, "\n********************* Extended *********************\n", cCount, cCountEx);
	for (int i = cCount; i < cCountEx; i++) {
		Cell & c = cells[i];
		fprintf(fp, "============================================================\n");
		fprintf(fp, "#%06d\n", i);
		fprintf(fp, "NODES: %6d %6d %6d %6d\n", c.nodesInd[0], c.nodesInd[1], c.nodesInd[2], c.nodesInd[3]);
		fprintf(fp, "EDGES: %6d %6d %6d %6d\n", c.facesInd[0], c.facesInd[1], c.facesInd[2], c.facesInd[3]);
		fprintf(fp, "CENTER: %25.15e %25.15e %25.15e\n", c.c.x, c.c.y, c.c.z);
		fprintf(fp, "HX: %25.15e   HY: %25.15e   HZ: %25.15e\n", c.HX, c.HY, c.HZ);
		fprintf(fp, "VOLUME: %25.15e\n", c.V);
		fprintf(fp, "TYPE: %s\n", c.typeName);
		fprintf(fp, "============================================================\n\n\n");
	}
	fclose(fp);

	sprintf(name, "mesh.%04d.faces.info", Parallel::procId);
	fp = fopen(name, "w");
	fprintf(fp, "COUNT: %6d    COUNT_EX: %6d\n", fCount, fCountEx);
	for (int i = 0; i < fCount; i++) {
		Face & f = faces[i];
		fprintf(fp, "============================================================\n");
		fprintf(fp, "#%06d\n", i);
		fprintf(fp, "NODES: %6d %6d %6d\n", f.nodesInd[0], f.nodesInd[1], f.nodesInd[2]);
		fprintf(fp, "CELLS: %6d %6d\n", f.c1, f.c2);
		fprintf(fp, "CENTER: %25.15e %25.15e %25.15e\n", f.c.x, f.c.y, f.c.z);
		fprintf(fp, "SQUARE: %25.15e\n", f.S);
		fprintf(fp, "NORMAL: %25.15e %25.15e %25.15e\n", f.n.x, f.n.y, f.n.z);
		fprintf(fp, "TYPE: %s\n", f.typeName);
		fprintf(fp, "============================================================\n\n\n");
	}
	fprintf(fp, "\n********************* Extended *********************\n", cCount, cCountEx);
	for (int i = fCount; i < fCountEx; i++) {
		Face& f = faces[i];
		fprintf(fp, "============================================================\n");
		fprintf(fp, "#%06d\n", i);
		fprintf(fp, "NODES: %6d %6d %6d\n", f.nodesInd[0], f.nodesInd[1], f.nodesInd[2]);
		fprintf(fp, "CELLS: %6d %6d\n", f.c1, f.c2);
		fprintf(fp, "CENTER: %25.15e %25.15e %25.15e\n", f.c.x, f.c.y, f.c.z);
		fprintf(fp, "SQUARE: %25.15e\n", f.S);
		fprintf(fp, "NORMAL: %25.15e %25.15e %25.15e\n", f.n.x, f.n.y, f.n.z);
		fprintf(fp, "TYPE: %s\n", f.typeName);
		fprintf(fp, "============================================================\n\n\n");
	}
	fclose(fp);
}