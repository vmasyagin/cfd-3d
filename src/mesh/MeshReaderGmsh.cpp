#include "MeshReaderGmsh.h"
#include <iostream>
#include <fstream> // подключаем файлы
#include <algorithm>

std::vector<int> MeshReaderGmsh::getIntersection(index_list &sets)
{
	std::vector <int> result;  // To store the reaultant set 
	int smallSetInd = 0;  // Initialize index of smallest set 
	int minSize = sets[0].size(); // Initialize size of smallest set 

	// sort all the sets, and also find the smallest set 
	for (int i = 1; i < sets.size(); i++)
	{
		// sort this set 
		std::sort(sets[i].begin(), sets[i].end());

		// update minSize, if needed 
		if (minSize > sets[i].size())
		{
			minSize = sets[i].size();
			smallSetInd = i;
		}
	}

	std::map<int, int> elementsMap;

	// Add all the elements of smallest set to a map, if already present, 
	// update the frequency 
	for (int i = 0; i < sets[smallSetInd].size(); i++)
	{
		if (elementsMap.find(sets[smallSetInd][i]) == elementsMap.end())
			elementsMap[sets[smallSetInd][i]] = 1;
		else
			elementsMap[sets[smallSetInd][i]]++;
	}

	// iterate through the map elements to see if they are present in 
	// remaining sets 
	std::map<int, int>::iterator it;
	for (it = elementsMap.begin(); it != elementsMap.end(); ++it)
	{
		int elem = it->first;
		int freq = it->second;

		bool bFound = true;

		// Iterate through all sets 
		for (int j = 0; j < sets.size(); j++)
		{
			// If this set is not the smallest set, then do binary search in it 
			if (j != smallSetInd)
			{
				// If the element is found in this set, then find its frequency 
				if (binary_search(sets[j].begin(), sets[j].end(), elem))
				{
					int lInd = lower_bound(sets[j].begin(), sets[j].end(), elem)
						- sets[j].begin();
					int rInd = upper_bound(sets[j].begin(), sets[j].end(), elem)
						- sets[j].begin();

					// Update the minimum frequency, if needed 
					if ((rInd - lInd) < freq)
						freq = rInd - lInd;
				}
				// If the element is not present in any set, then no need  
				// to proceed for this element. 
				else
				{
					bFound = false;
					break;
				}
			}
		}

		// If element was found in all sets, then add it to result 'freq' times 
		if (bFound)
		{
			for (int k = 0; k < freq; k++)
				result.push_back(elem);
		}
	}
	return result;
}

int MeshReaderGmsh::findFace(int n1, int n2, int n3)
{
	std::set <int> s1 = { n1, n2, n3 };
	std::set <int> s2;
	for (int i = 0; i < bnd_faces.size(); i++) {
		s2 = { bnd_faces[i][1], bnd_faces[i][2], bnd_faces[i][3] };
		if (s1 == s2) {
			return bnd_faces[i][0];
		}
		s2.clear();
	}
	return -1;
}

void MeshReaderGmsh::read(Grid* g)
{
	char str[64];
	long long int tmp, elem, cell, face, type;
	double tmp1;
	int numPatches, numElements, retval;
	std::string line;

	// читаем сеточные данные
	std::sprintf(str, "%s.msh", fileName);
	std::ifstream file(str);
	string_list patches;
	while (getline(file, line)) {
		if (line[0] == '$') {
			if (line == "$PhysicalNames") {
				getline(file, line);
				sscanf(line.c_str(), "%d", &numPatches);
				patches.resize(numPatches);
				for (int i = 0; i < numPatches; i++) {
					getline(file, line);
					int dim, id;
					char name[128];
					sscanf(line.c_str(), "%d %d \"%[^\"]", &(dim), &(id), name);
					patches[--id] = name;
				}
			}
			else if (line == "$Nodes") {
				getline(file, line);
				sscanf(line.c_str(), "%d", &(g->nCount));
				g->nodes = new Point[g->nCount];
				for (int i = 0; i < g->nCount; i++) {
					getline(file, line);
					sscanf(line.c_str(), "%lld %lf %lf %lf", &tmp, &g->nodes[i].x, &g->nodes[i].y, &g->nodes[i].z);
				}
			}
			else if (line == "$Elements") {
				getline(file, line);
				sscanf(line.c_str(), "%d", &numElements);
				for (int i = 0; i < numElements; i++) {
					getline(file, line);
					retval = sscanf(line.c_str(), "%lld %lld", &elem, &type);
					if (retval != 2) {
						log("Premature end of file");
						exit(2);
					}

					if (type == 2) {
						indexes v;
						v.clear();
						int node1, node2, node3;
						int tagc, tag1, tag2;
						retval = sscanf(line.c_str(), "%lld %lld %d %d %d %d %d %d", &face, &type, &tagc, &tag1, &tag2,
							&node1, &node2, &node3);
						if (retval != 8) {
							log("Premature end of file");
							exit(2);
						}
						v.push_back(type);
						v.push_back(tagc);
						v.push_back(tag1);
						v.push_back(tag2);
						v.push_back(node1);
						v.push_back(node2);
						v.push_back(node3);
						faces.push_back(v);
					}

					if (type == 4) {
						indexes v;
						v.clear();
						int node1, node2, node3, node4;
						int tagc, tag1, tag2;
						retval = sscanf(line.c_str(), "%lld %lld %d %d %d %d %d %d %d", &cell, &type, &tagc, &tag1, &tag2,
							&node1, &node2, &node3, &node4);
						if (retval != 9) {
							log("Premature end of file");
							exit(2);
						}
						v.push_back(type);
						v.push_back(tagc);
						v.push_back(tag1);
						v.push_back(tag2);
						v.push_back(node1);
						v.push_back(node2);
						v.push_back(node3);
						v.push_back(node4);
						cells.push_back(v);
					}
				}
			}
		}
	}
	file.close();
	
	std::map<int, ind_set> node_cells;
	// читаем данные о ЯЧЕЙКАХ
	int iCell = 0;
	g->cCount = cells.size();
	g->cells = new Cell[g->cCount];
	for (index_list::iterator it = cells.begin(); it != cells.end(); it++)
	{
		g->cells[iCell].nCount = 4;
		g->cells[iCell].neigh = new int[4];
		g->cells[iCell].nodesInd = new int[g->cells[iCell].nCount];
		g->cells[iCell].nodesInd[0] = --(*it)[4];
		g->cells[iCell].nodesInd[1] = --(*it)[5];
		g->cells[iCell].nodesInd[2] = --(*it)[6];
		g->cells[iCell].nodesInd[3] = --(*it)[7];
		g->cells[iCell].type = (*it)[2];

		node_cells[g->cells[iCell].nodesInd[0]].insert(iCell);
		node_cells[g->cells[iCell].nodesInd[1]].insert(iCell);
		node_cells[g->cells[iCell].nodesInd[2]].insert(iCell);
		node_cells[g->cells[iCell].nodesInd[3]].insert(iCell);
	
		strcpy(g->cells[iCell].typeName, patches[g->cells[iCell].type - 1].c_str());
		g->cells[iCell].c.x = 0.25 * (g->nodes[g->cells[iCell].nodesInd[0]].x + g->nodes[g->cells[iCell].nodesInd[1]].x + g->nodes[g->cells[iCell].nodesInd[2]].x + g->nodes[g->cells[iCell].nodesInd[3]].x);
		g->cells[iCell].c.y = 0.25 * (g->nodes[g->cells[iCell].nodesInd[0]].y + g->nodes[g->cells[iCell].nodesInd[1]].y + g->nodes[g->cells[iCell].nodesInd[2]].y + g->nodes[g->cells[iCell].nodesInd[3]].y);
		g->cells[iCell].c.z = 0.25 * (g->nodes[g->cells[iCell].nodesInd[0]].z + g->nodes[g->cells[iCell].nodesInd[1]].z + g->nodes[g->cells[iCell].nodesInd[2]].z + g->nodes[g->cells[iCell].nodesInd[3]].z);

		g->cells[iCell].HX = _max_(_max_(g->nodes[g->cells[iCell].nodesInd[0]].x, g->nodes[g->cells[iCell].nodesInd[1]].x), _max_(g->nodes[g->cells[iCell].nodesInd[2]].x, g->nodes[g->cells[iCell].nodesInd[3]].x)) -
			_min_(_min_(g->nodes[g->cells[iCell].nodesInd[0]].x, g->nodes[g->cells[iCell].nodesInd[1]].x), _min_(g->nodes[g->cells[iCell].nodesInd[2]].x, g->nodes[g->cells[iCell].nodesInd[3]].x));
		g->cells[iCell].HY = _max_(_max_(g->nodes[g->cells[iCell].nodesInd[0]].y, g->nodes[g->cells[iCell].nodesInd[1]].y), _max_(g->nodes[g->cells[iCell].nodesInd[2]].y, g->nodes[g->cells[iCell].nodesInd[3]].y)) -
			_min_(_min_(g->nodes[g->cells[iCell].nodesInd[0]].y, g->nodes[g->cells[iCell].nodesInd[1]].y), _min_(g->nodes[g->cells[iCell].nodesInd[2]].y, g->nodes[g->cells[iCell].nodesInd[3]].y));
		g->cells[iCell].HZ = _max_(_max_(g->nodes[g->cells[iCell].nodesInd[0]].z, g->nodes[g->cells[iCell].nodesInd[1]].z), _max_(g->nodes[g->cells[iCell].nodesInd[2]].z, g->nodes[g->cells[iCell].nodesInd[3]].z)) -
			_min_(_min_(g->nodes[g->cells[iCell].nodesInd[0]].z, g->nodes[g->cells[iCell].nodesInd[1]].z), _min_(g->nodes[g->cells[iCell].nodesInd[2]].z, g->nodes[g->cells[iCell].nodesInd[3]].z));
		g->cells[iCell].fCount = 4;
		g->cells[iCell].facesInd = new int[g->cells[iCell].fCount];

		iCell++;
	}
	
	// определяем соседние ячейки
	int** neigh;
	neigh = new int*[g->cCount];
	for (int i = 0; i < g->cCount; i++)
	{
		neigh[i] = new int[4];
		Cell & c = g->cells[i];
		for (int k = 0; k < 4; k++) {
			index_list sets;
			
			sets.push_back(indexes(node_cells[c.nodesInd[(k + 1) % 4]].begin(), node_cells[c.nodesInd[(k + 1) % 4]].end()));
			sets.push_back(indexes(node_cells[c.nodesInd[(k + 2) % 4]].begin(), node_cells[c.nodesInd[(k + 2) % 4]].end()));
			sets.push_back(indexes(node_cells[c.nodesInd[(k + 3) % 4]].begin(), node_cells[c.nodesInd[(k + 3) % 4]].end()));

			indexes res = getIntersection(sets);
			sets.clear();

			c.neigh[k] = -2;
			for (indexes::iterator rit = res.begin(); rit != res.end(); rit++) {
				if (*rit != i) {
					c.neigh[k] = *rit;
				}
			}
			neigh[i][k] = c.neigh[k];
		}
		g->cells[i].neigh[0] = neigh[i][0];
		g->cells[i].neigh[1] = neigh[i][1];
		g->cells[i].neigh[2] = neigh[i][2];
		g->cells[i].neigh[3] = neigh[i][3];
	}

	// формируем данные о ГРАНЯХ
	g->fCount = 0;
	for (int i = 0; i < g->cCount; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			int p = neigh[i][j];
			if (p > -1)
			{
				for (int k = 0; k < 4; k++)
				{ // убираем у соседа номер этой ячейки, чтобы грань не повторялась
					if (neigh[p][k] == i) neigh[p][k] = -1;
				}
				g->fCount++;
			}
			if (p == -2) g->fCount++;
		}
	}
	g->fCountEx = g->fCount;
	g->faces = new Face[g->fCount];

	int iFace = 0;
	int* cfi = new int[g->cCount];
	for (int i = 0; i < g->cCount; i++)
	{
		cfi[i] = 0;
	}
	// ::memset(cfi, 0, cCount*sizeof(int));
	std::vector<int> v;

	for (int i = 0; i < g->cCount; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			int p = neigh[i][j];
			if (p != -1)
			{
				g->faces[iFace].nCount = 3;
				g->faces[iFace].nodesInd = new int[g->faces[iFace].nCount];
				g->faces[iFace].edgesInd = new int[g->faces[iFace].nCount];
				g->faces[iFace].nodesInd[0] = g->cells[i].nodesInd[(j + 1) % 4];
				g->faces[iFace].nodesInd[1] = g->cells[i].nodesInd[(j + 2) % 4];
				g->faces[iFace].nodesInd[2] = g->cells[i].nodesInd[(j + 3) % 4];

				g->faces[iFace].cCount = 1;
				g->faces[iFace].c.x = (g->nodes[g->faces[iFace].nodesInd[0]].x + g->nodes[g->faces[iFace].nodesInd[1]].x + g->nodes[g->faces[iFace].nodesInd[2]].x) / 3.0;
				g->faces[iFace].c.y = (g->nodes[g->faces[iFace].nodesInd[0]].y + g->nodes[g->faces[iFace].nodesInd[1]].y + g->nodes[g->faces[iFace].nodesInd[2]].y) / 3.0;
				g->faces[iFace].c.z = (g->nodes[g->faces[iFace].nodesInd[0]].z + g->nodes[g->faces[iFace].nodesInd[1]].z + g->nodes[g->faces[iFace].nodesInd[2]].z) / 3.0;

				double x1 = g->nodes[g->faces[iFace].nodesInd[0]].x; double y1 = g->nodes[g->faces[iFace].nodesInd[0]].y; double z1 = g->nodes[g->faces[iFace].nodesInd[0]].z;
				double x2 = g->nodes[g->faces[iFace].nodesInd[1]].x; double y2 = g->nodes[g->faces[iFace].nodesInd[1]].y; double z2 = g->nodes[g->faces[iFace].nodesInd[1]].z;
				double x3 = g->nodes[g->faces[iFace].nodesInd[2]].x; double y3 = g->nodes[g->faces[iFace].nodesInd[2]].y; double z3 = g->nodes[g->faces[iFace].nodesInd[2]].z;

				double a1 = x1 - x3; double a2 = x2 - x3; double a3 = x3;
				double b1 = y1 - y3; double b2 = y2 - y3; double b3 = y3;
				double c1 = z1 - z3; double c2 = z2 - z3; double c3 = z3;

				g->faces[iFace].n.x = b1 * c2 - c1 * b2;
				g->faces[iFace].n.y = -(a1 * c2 - c1 * a2);
				g->faces[iFace].n.z = a1 * b2 - b1 * a2;
				g->faces[iFace].n.x /= sqrt((b1 * c2 - c1 * b2) * (b1 * c2 - c1 * b2) + (a1 * c2 - c1 * a2) * (a1 * c2 - c1 * a2) + (a1 * b2 - b1 * a2) * (a1 * b2 - b1 * a2));
				g->faces[iFace].n.y /= sqrt((b1 * c2 - c1 * b2) * (b1 * c2 - c1 * b2) + (a1 * c2 - c1 * a2) * (a1 * c2 - c1 * a2) + (a1 * b2 - b1 * a2) * (a1 * b2 - b1 * a2));
				g->faces[iFace].n.z /= sqrt((b1 * c2 - c1 * b2) * (b1 * c2 - c1 * b2) + (a1 * c2 - c1 * a2) * (a1 * c2 - c1 * a2) + (a1 * b2 - b1 * a2) * (a1 * b2 - b1 * a2));

				double x4 = g->nodes[g->cells[i].nodesInd[j]].x; double y4 = g->nodes[g->cells[i].nodesInd[j]].y; double z4 = g->nodes[g->cells[i].nodesInd[j]].z;
				if ((x4 - x3) * g->faces[iFace].n.x + (y4 - y3) * g->faces[iFace].n.y + (z4 - z3) * g->faces[iFace].n.z > 0)
				{
					g->faces[iFace].n.x *= -1.;
					g->faces[iFace].n.y *= -1.;
					g->faces[iFace].n.z *= -1.;
				}

				double a = sqrt(_sqr_(x2 - x1) + _sqr_(y2 - y1) + _sqr_(z2 - z1));
				double b = sqrt(_sqr_(x3 - x1) + _sqr_(y3 - y1) + _sqr_(z3 - z1));
				double c = sqrt(_sqr_(x3 - x2) + _sqr_(y3 - y2) + _sqr_(z3 - z2));
				double hp = 0.5 * (a + b + c);

				g->faces[iFace].S = sqrt(hp * (hp - a) * (hp - b) * (hp - c));

				g->faces[iFace].c1 = i;
				g->cells[i].facesInd[cfi[i]] = iFace;
				cfi[i]++;

				if (p > -1)
				{

					g->faces[iFace].c2 = p;
					g->cells[p].facesInd[cfi[p]] = iFace;
					cfi[p]++;
					//g->faces[iFace].cnl2 = fabs(g->faces[iFace].n.x*(g->cells[g->faces[iFace].c2].c.x - g->faces[iFace].c[0].x) + g->faces[iFace].n.y*(g->cells[g->faces[iFace].c2].c.y - g->faces[iFace].c[0].y));
					g->faces[iFace].type = Face::TYPE_INNER;
					sprintf(g->faces[iFace].typeName, "");
				}
				if (p == -2)
				{
					g->faces[iFace].c2 = -1;
					//g->faces[iFace].cnl2 = 0;
					g->faces[iFace].type = -1;

					v.clear();
					v.push_back(iFace);
					v.push_back(g->faces[iFace].nodesInd[0]);
					v.push_back(g->faces[iFace].nodesInd[1]);
					v.push_back(g->faces[iFace].nodesInd[2]);
					bnd_faces.push_back(v);
				}
				iFace++;
			}
		}
	}

	for (index_list::iterator it = faces.begin(); it != faces.end(); it++)
	{
		int iFace = findFace(--(*it)[4], --(*it)[5], --(*it)[6]);
		if (iFace > -1) {
			sprintf(g->faces[iFace].typeName, "%s", patches[(*it)[2] - 1].c_str());
		}
		else {
			throw Exception((char*)"Boundary face #%d not defined in msh file.", Exception::TYPE_MESH_GMSH_NOT_DEFINED_BND_FACE);
		}
	}

	// Рассчитываем характерный размер ячейки
	for (int i = 0; i < g->cCount; i++)
	{
		Point A = g->nodes[g->cells[i].nodesInd[0]];
		Point B = g->nodes[g->cells[i].nodesInd[1]];
		Point C = g->nodes[g->cells[i].nodesInd[2]];
		Point D = g->nodes[g->cells[i].nodesInd[3]];
		g->cells[i].V = fabs(mixed_prod(A - D, B - D, C - D)) / 6.;
		double HSize = 1.0e10;
		Point ptCellC = g->cells[i].c;
		for (int j = 0; j < g->cells[i].fCount; j++) {
			Point ptFaceC = g->faces[g->cells[i].facesInd[j]].c;	// получили центр грани
			double R = _calc_dist_2_(ptFaceC, ptCellC);
			if (R < HSize) HSize = R;
		}
		g->cells[i].H = 2.0 * sqrt(HSize);
	}

	for (int i = 0; i < g->cCount; i++)
	{
		delete[] neigh[i];
	}
	delete[] neigh;
	delete[] cfi;
}


