#include "MeshReaderSalomeUnv.h"
#include <fstream>
#include <algorithm>
#include <cstring>
#include <ctime>


int MeshReaderSalomeUnv::find_index(int val)
{
	return elements[val].ind;
}


int MeshReaderSalomeUnv::find_edge(int n1, int n2)
{
	for (int i = 0; i < edges.size(); i++) {
		if (((edges[i][3] == n1) && (edges[i][4] == n2)) || ((edges[i][3] == n2) && (edges[i][4] == n1)) ){
			return i;
		}
	}
	return -1;
}


void MeshReaderSalomeUnv::read(Grid * g)
{
	// чтение файла в контейнеры
	ifstream fin(fileName);
	if (fin.is_open()) {

		string_list sl;
		while (!fin.eof()) {
			read_block(&sl, fin);
			parse_block(sl, g);
		}
		fin.close();
	}
	else {
		throw Exception((char*)"File '%s' opening error.", Exception::FILE_OPENING_ERROR);
	}

	/* 
	 *  Упаковка данных в класс Grid
	 */

	int i;

	log((char*)"Building mesh structure:\n");

	// nodes
	log((char*)"\t- nodes;\n");
	g->nCount = points.size();
	g->nodes = new Point[g->nCount];
	i = 0;
	for (vector<Point>::iterator it = points.begin(); it != points.end(); it++, i++) {
		Point & p = g->nodes[i];
		p.x = it->x;
		p.y = it->y;
	}

	log((char*)"\t- cells;\n");
	g->cCount = cells.size();
	g->cells = new Cell[g->cCount];
	map<int, ind_set> node_cells;
	i = 0;
	for (index_list::iterator it = cells.begin(); it != cells.end(); it++, i++) {
		indexes & ind = *it;
		node_cells[ind[3]].insert(i);
		node_cells[ind[4]].insert(i);
		node_cells[ind[5]].insert(i);
	}

	i = 0;
	int** neigh = new int*[g->cCount];
	for (index_list::iterator it = cells.begin(); it != cells.end(); it++, i++) {
		neigh[i] = new int[3];
		indexes & ind = *it;
		
		elements[ind[0]] = element(i, element::TYPE_CELL);

		Cell & c = g->cells[i];
		c.nCount = 3;
		c.nodesInd = new int[g->cells[i].nCount];
		c.nodesInd[0] = ind[3];
		c.nodesInd[1] = ind[4];
		c.nodesInd[2] = ind[5];

		g->cells[i].type = 0;
		g->cells[i].c.x = (g->nodes[g->cells[i].nodesInd[0]].x + g->nodes[g->cells[i].nodesInd[1]].x + g->nodes[g->cells[i].nodesInd[2]].x) / 3.0;
		g->cells[i].c.y = (g->nodes[g->cells[i].nodesInd[0]].y + g->nodes[g->cells[i].nodesInd[1]].y + g->nodes[g->cells[i].nodesInd[2]].y) / 3.0;
		g->cells[i].HX = _max_(fabs(g->nodes[g->cells[i].nodesInd[0]].x - g->nodes[g->cells[i].nodesInd[1]].x),
			fabs(g->nodes[g->cells[i].nodesInd[1]].x - g->nodes[g->cells[i].nodesInd[2]].x),
			fabs(g->nodes[g->cells[i].nodesInd[0]].x - g->nodes[g->cells[i].nodesInd[2]].x));
		g->cells[i].HY = _max_(fabs(g->nodes[g->cells[i].nodesInd[0]].y - g->nodes[g->cells[i].nodesInd[1]].y),
			fabs(g->nodes[g->cells[i].nodesInd[1]].y - g->nodes[g->cells[i].nodesInd[2]].y),
			fabs(g->nodes[g->cells[i].nodesInd[0]].y - g->nodes[g->cells[i].nodesInd[2]].y));
		g->cells[i].eCount = 3;
		g->cells[i].edgesInd = new int[g->cells[i].eCount];

		for (int k = 0; k < 3; k++) {
			map<int, ind_set>::iterator out_it;
			indexes res(10);
			indexes::iterator it_res = set_intersection(
				node_cells[c.nodesInd[k % 3]].begin(), node_cells[c.nodesInd[k % 3]].end(),
				node_cells[c.nodesInd[(k + 1) % 3]].begin(), node_cells[c.nodesInd[(k + 1) % 3]].end(),
				res.begin());
			res.resize(it_res-res.begin());
			c.neigh[k] = -2;
			for (indexes::iterator rit = res.begin(); rit != res.end(); rit++) {
				if (*rit != i) {
					c.neigh[k] = *rit;
				}
			}
			neigh[i][k] = c.neigh[k];
		}
	}

	log((char*)"\t- edges;\n");
	g->eCount = 0;
	for (int i = 0; i < g->cCount; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			int p = neigh[i][j];
			if (p > -1)
			{
				for (int k = 0; k < 3; k++)
				{ // убираем у соседа номер этой ячейки, чтобы грань не повторялась
					if (neigh[p][k] == i) neigh[p][k] = -1;
				}
				g->eCount++;
			}
			if (p == -2) g->eCount++;
		}
	}
	g->edges = new Edge[g->eCount];

	int iEdge = 0;
	int * cfi = new int[g->cCount];
	for (int i = 0; i < g->cCount; i++)
	{
		cfi[i] = 0;
	}
	// ::memset(cfi, 0, cCount*sizeof(int));
	for (int i = 0; i < g->cCount; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			int p = neigh[i][j];
			if (p != -1)
			{
				g->edges[iEdge].n1 = g->cells[i].nodesInd[(j + 0) % 3];
				g->edges[iEdge].n2 = g->cells[i].nodesInd[(j + 1) % 3];
				//g->edges[iEdge].cCount = 1;
				//g->edges[iEdge].c = new Point[g->edges[iEdge].cCount];
				//g->edges[iEdge].c[0].x = (g->nodes[g->edges[iEdge].n1].x + g->nodes[g->edges[iEdge].n2].x) / 2.0;
				//g->edges[iEdge].c[0].y = (g->nodes[g->edges[iEdge].n1].y + g->nodes[g->edges[iEdge].n2].y) / 2.0;
				g->edges[iEdge].cCount = 3;
				g->edges[iEdge].c = new Point[g->edges[iEdge].cCount];
				double _sqrt3 = 1.0 / sqrt(3.0);
				// центр ребра
				g->edges[iEdge].c[0].x = (g->nodes[g->edges[iEdge].n1].x + g->nodes[g->edges[iEdge].n2].x) / 2.0;
				g->edges[iEdge].c[0].y = (g->nodes[g->edges[iEdge].n1].y + g->nodes[g->edges[iEdge].n2].y) / 2.0;
				// первая точка Гаусса
				g->edges[iEdge].c[1].x = (g->nodes[g->edges[iEdge].n1].x + g->nodes[g->edges[iEdge].n2].x) / 2.0 - _sqrt3*(g->nodes[g->edges[iEdge].n2].x - g->nodes[g->edges[iEdge].n1].x) / 2.0;
				g->edges[iEdge].c[1].y = (g->nodes[g->edges[iEdge].n1].y + g->nodes[g->edges[iEdge].n2].y) / 2.0 - _sqrt3*(g->nodes[g->edges[iEdge].n2].y - g->nodes[g->edges[iEdge].n1].y) / 2.0;
				// вторая точка Гаусса
				g->edges[iEdge].c[2].x = (g->nodes[g->edges[iEdge].n1].x + g->nodes[g->edges[iEdge].n2].x) / 2.0 + _sqrt3*(g->nodes[g->edges[iEdge].n2].x - g->nodes[g->edges[iEdge].n1].x) / 2.0;
				g->edges[iEdge].c[2].y = (g->nodes[g->edges[iEdge].n1].y + g->nodes[g->edges[iEdge].n2].y) / 2.0 + _sqrt3*(g->nodes[g->edges[iEdge].n2].y - g->nodes[g->edges[iEdge].n1].y) / 2.0;
				g->edges[iEdge].n.x = g->nodes[g->edges[iEdge].n2].y - g->nodes[g->edges[iEdge].n1].y;
				g->edges[iEdge].n.y = g->nodes[g->edges[iEdge].n1].x - g->nodes[g->edges[iEdge].n2].x;
				g->edges[iEdge].l = sqrt(g->edges[iEdge].n.x*g->edges[iEdge].n.x + g->edges[iEdge].n.y*g->edges[iEdge].n.y);
				g->edges[iEdge].n.x /= g->edges[iEdge].l;
				g->edges[iEdge].n.y /= g->edges[iEdge].l;
				g->edges[iEdge].c1 = i;
				g->cells[i].edgesInd[cfi[i]] = iEdge;
				cfi[i]++;
				g->edges[iEdge].cnl1 = fabs(g->edges[iEdge].n.x*(g->edges[iEdge].c[0].x - g->cells[g->edges[iEdge].c1].c.x) + g->edges[iEdge].n.y*(g->edges[iEdge].c[0].y - g->cells[g->edges[iEdge].c1].c.y));

				// коррекция направлений нормалей
				Vector vc;
				vc.x = g->cells[i].c.x - g->edges[iEdge].c[0].x;
				vc.y = g->cells[i].c.y - g->edges[iEdge].c[0].y;
				if (scalar_prod(vc, g->edges[iEdge].n) > 0) {
					g->edges[iEdge].n.x *= -1;
					g->edges[iEdge].n.y *= -1;
				}

				if (p > -1)
				{

					g->edges[iEdge].c2 = p;
					g->cells[p].edgesInd[cfi[p]] = iEdge;
					cfi[p]++;
					g->edges[iEdge].cnl2 = fabs(g->edges[iEdge].n.x*(g->cells[g->edges[iEdge].c2].c.x - g->edges[iEdge].c[0].x) + g->edges[iEdge].n.y*(g->cells[g->edges[iEdge].c2].c.y - g->edges[iEdge].c[0].y));
					g->edges[iEdge].type = Edge::TYPE_INNER;
				}
				if (p == -2)
				{
					g->edges[iEdge].c2 = -1;
					g->edges[iEdge].cnl2 = 0;
					g->edges[iEdge].type = Edge::TYPE_NAMED;
					
					int jEdge = find_edge(g->edges[iEdge].n1, g->edges[iEdge].n2);
					if (jEdge != -1) {
						elements[edges[jEdge][0]] = element(iEdge, element::TYPE_EDGE);
					}
					else {
						throw Exception((char*)"Boundary edge #%d not defined in UNV file.", Exception::TYPE_MESH_UNV_NOT_DEFINED_BND_EDGE);
					}
				}


				iEdge++;
			}
		}
	}

	log((char*)"\t- bounds;\n");
	for (bnd_map::iterator it = bounds.begin(); it != bounds.end(); it++) {
		char name[64];
		strcpy(name, it->first.c_str());
		for (indexes::iterator ind = it->second.begin(); ind != it->second.end(); ind++) {
			int iEdge, iCell;
			switch (elements[*ind].type) {
			case element::TYPE_EDGE:
				iEdge = elements[*ind].ind;
				if (iEdge >= 0) {
					strcpy(g->edges[iEdge].typeName, name);
					g->edges[iEdge].type = Edge::TYPE_NAMED;
				}
				else {
					throw Exception((char*)"Wrong edge's index when boundary parsing...", -1);
				}
				break;
			case element::TYPE_CELL:
				iCell = elements[*ind].ind;
				if (iCell >= 0) {
					strcpy(g->cells[iCell].typeName, name);
					g->cells[iCell].type = -1;
				}
				else {
					throw Exception((char*)"Wrong edge's index when boundary parsing...", -1);
				}
				break;
			}
		}
	}

	for (int i = 0; i < g->cCount; i++)	{
		double a = g->edges[g->cells[i].edgesInd[0]].l;
		double b = g->edges[g->cells[i].edgesInd[1]].l;
		double c = g->edges[g->cells[i].edgesInd[2]].l;
		double p = (a + b + c) / 2.0;
		g->cells[i].S = sqrt(p*(p - a)*(p - b)*(p - c));

		delete[] neigh[i];
	}

	delete[] neigh;
	delete[] cfi;
	log((char*)"  complete...\n");

}


void MeshReaderSalomeUnv::read_block(string_list * sl, ifstream& fin)
{
	sl->clear();
	int minusCnt = 0;
	char s[255];
	string str;
	while (minusCnt < 2 && !fin.eof()) {
		fin.getline(s, 256);
		str.assign(s);
		int pos = str.find("-1");
		if (pos != string::npos) {
			if (pos == str.length() - 2) {
				minusCnt++;
				continue;
			}
		}

		sl->push_back(string(s));
	}
}


void MeshReaderSalomeUnv::parse_block_164(string_list sl, Grid * g)
{
	log((char*)"UNV reader: Block 167 skipped...");
	
	long time_start, time_end;
	time_start = clock();

	time_end = clock();
	log((char*)"\tTime: %d\n", time_end - time_start);
}


void MeshReaderSalomeUnv::parse_block_2420(string_list sl, Grid * g)
{
	log((char*)"UNV reader: Block 2420 skipped...");

	long time_start, time_end;
	time_start = clock();

	time_end = clock();
	log((char*)"\tTime: %d\n", time_end - time_start);
}


void MeshReaderSalomeUnv::parse_block_2411(string_list sl, Grid * g)
{
	log((char*)"UNV reader: Parsing block 2411...");

	long time_start, time_end;
	time_start = clock();

	points.clear();
	string_list::iterator it = sl.begin();
	++it;
	while (it != sl.end()) {
		Point p;
		double tmpd;
		int tmpi;
		sscanf(it->c_str(), "%d %d %d %d", &tmpi, &tmpi, &tmpi, &tmpi); it++;
		sscanf(it->c_str(), "%lf %lf %lf", &p.x, &p.y, &tmpd); it++;
		points.push_back(p);
	}

	time_end = clock();
	log((char*)"\tTime: %d\n", time_end-time_start);
}


void MeshReaderSalomeUnv::parse_block_2412(string_list sl, Grid * g)
{
	log((char*)"UNV reader: Parsing block 2412...");

	long time_start, time_end;
	time_start = clock();

	int maxInd = 0;
	edges.clear();
	cells.clear();
	string_list::iterator it = sl.begin();
	++it;
	while (it != sl.end()) {
		indexes p;
		int tmp[64];
		sscanf(it->c_str(), "%d %d %d %d %d %d", &tmp[0], &tmp[1], &tmp[2], &tmp[3], &tmp[4], &tmp[5]); it++;
		tmp[0]--;
		if (tmp[0] > maxInd) maxInd = tmp[0];
		p.push_back(tmp[0]);
		p.push_back(tmp[1]);
		p.push_back(tmp[5]);
		switch (tmp[1]) {
		case 11: // ребро
			sscanf(it->c_str(), "%d %d %d", &tmp[6], &tmp[7], &tmp[8]); it++;
			sscanf(it->c_str(), "%d %d", &tmp[9], &tmp[10]); it++;
			tmp[9]--;
			tmp[10]--;
			p.push_back(tmp[9]);
			p.push_back(tmp[10]);
			edges.push_back(p);
			break;

		case 41: // треугольник
			sscanf(it->c_str(), "%d %d %d", &tmp[6], &tmp[7], &tmp[8]); it++;
			tmp[6]--;
			tmp[7]--;
			tmp[8]--;
			p.push_back(tmp[6]);
			p.push_back(tmp[7]);
			p.push_back(tmp[8]);
			cells.push_back(p);
			break;
		
		default:
			char msg[64];
			sprintf(msg, "Unknown element type '%d'.", tmp[1]);
			throw Exception(msg, Exception::TYPE_MESH_UNV_UNKNOWN_ELEMENT);
		}
	}
	elements.resize(maxInd+1);

	time_end = clock();
	log((char*)"\tTime: %d\n", time_end - time_start);
}


void MeshReaderSalomeUnv::parse_block_2467(string_list sl, Grid * g)
{
	log((char*)"UNV reader: Parsing block 2467...");

	long time_start, time_end;
	time_start = clock();

	bounds.clear();
	string_list::iterator it = sl.begin();
	++it;
	while (it != sl.end()) {
		indexes p;
		int tmp[8];
		char bnd_name[128];
		sscanf(it->c_str(), "%d %d %d %d %d %d %d %d", &tmp[0], &tmp[1], &tmp[2], &tmp[3], &tmp[4], &tmp[5], &tmp[6], &tmp[7]); it++;
		int n = tmp[7];
		sscanf(it->c_str(), "%s", bnd_name); it++; //TODO &bnd_name
		p.clear();
		for (int i = 0; i < n/2; i++) {
			sscanf(it->c_str(), "%d %d %d %d %d %d %d %d", &tmp[0], &tmp[1], &tmp[2], &tmp[3], &tmp[4], &tmp[5], &tmp[6], &tmp[7]); it++;
			p.push_back(--tmp[1]);
			p.push_back(--tmp[5]);
		}
		if (n % 2 == 1) {
			sscanf(it->c_str(), "%d %d %d %d", &tmp[0], &tmp[1], &tmp[2], &tmp[3]); it++;
			p.push_back(--tmp[1]);
		}
		bounds[string(bnd_name)] = p;
	}

	time_end = clock();
	log((char*)"\tTime: %d\n", time_end - time_start);
}


void MeshReaderSalomeUnv::parse_block(string_list sl, Grid * g)
{
	int type = atoi(sl[0].c_str());
	switch (type) {
	case 164:
		parse_block_164(sl, g);
		break;
	case 2411:
		parse_block_2411(sl, g);
		break;
	case 2412:
		parse_block_2412(sl, g);
		break;
	case 2420:
		parse_block_2420(sl, g);
		break;
	case 2467:
		parse_block_2467(sl, g);
		break;
	}
}
