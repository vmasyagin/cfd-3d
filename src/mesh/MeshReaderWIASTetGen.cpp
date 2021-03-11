#include "MeshReaderWIASTetGen.h"

void MeshReaderWIASTetGen::read(Grid* g)
{
	char str[50];
	FILE *fp;
	int tmp;

	// читаем данные об УЗЛАХ
	sprintf(str, "%s.node", fileName);
	fp = fopen(str, "r");
	fscanf(fp, "%d %d %d %d", &g->nCount, &tmp, &tmp, &tmp);
	g->nodes = new Point[g->nCount];
	for (int i = 0; i < g->nCount; i++)
	{
		fscanf(fp, "%d %lf %lf %lf", &tmp, &(g->nodes[i].x), &(g->nodes[i].y), &(g->nodes[i].z));
	}
	fclose(fp);

	// читаем данные о ЯЧЕЙКАХ
	sprintf(str, "%s.ele", fileName);
	fp = fopen(str, "r");
	fscanf(fp, "%d %d %d", &g->cCount, &tmp, &tmp);
	g->cells = new Cell[g->cCount];
	for (int i = 0; i < g->cCount; i++)
	{
		g->cells[i].nCount = 4;
		g->cells[i].nodesInd = new int[g->cells[i].nCount];
		fscanf(fp, "%d %d %d %d %d %d", &tmp, &(g->cells[i].nodesInd[0]), &(g->cells[i].nodesInd[1]), &(g->cells[i].nodesInd[2]), &(g->cells[i].nodesInd[3]), &(g->cells[i].type));
		g->cells[i].nodesInd[0]--;
		g->cells[i].nodesInd[1]--;
		g->cells[i].nodesInd[2]--;
		g->cells[i].nodesInd[3]--;
		sprintf(g->cells[i].typeName, "%d", g->cells[i].type);
		g->cells[i].c.x = (g->nodes[g->cells[i].nodesInd[0]].x + g->nodes[g->cells[i].nodesInd[1]].x + g->nodes[g->cells[i].nodesInd[2]].x + g->nodes[g->cells[i].nodesInd[3]].x) / 4.0;
		g->cells[i].c.y = (g->nodes[g->cells[i].nodesInd[0]].y + g->nodes[g->cells[i].nodesInd[1]].y + g->nodes[g->cells[i].nodesInd[2]].y + g->nodes[g->cells[i].nodesInd[3]].y) / 4.0;
		g->cells[i].HX = _max_(_max_(g->nodes[g->cells[i].nodesInd[0]].x, g->nodes[g->cells[i].nodesInd[1]].x), _max_(g->nodes[g->cells[i].nodesInd[2]].x, g->nodes[g->cells[i].nodesInd[3]].x)) -
			_min_(_min_(g->nodes[g->cells[i].nodesInd[0]].x, g->nodes[g->cells[i].nodesInd[1]].x), _min_(g->nodes[g->cells[i].nodesInd[2]].x, g->nodes[g->cells[i].nodesInd[3]].x));
		g->cells[i].HY = _max_(_max_(g->nodes[g->cells[i].nodesInd[0]].y, g->nodes[g->cells[i].nodesInd[1]].y), _max_(g->nodes[g->cells[i].nodesInd[2]].y, g->nodes[g->cells[i].nodesInd[3]].y)) -
			_min_(_min_(g->nodes[g->cells[i].nodesInd[0]].y, g->nodes[g->cells[i].nodesInd[1]].y), _min_(g->nodes[g->cells[i].nodesInd[2]].y, g->nodes[g->cells[i].nodesInd[3]].y));
		g->cells[i].HZ = _max_(_max_(g->nodes[g->cells[i].nodesInd[0]].z, g->nodes[g->cells[i].nodesInd[1]].z), _max_(g->nodes[g->cells[i].nodesInd[2]].z, g->nodes[g->cells[i].nodesInd[3]].z)) -
			_min_(_min_(g->nodes[g->cells[i].nodesInd[0]].z, g->nodes[g->cells[i].nodesInd[1]].z), _min_(g->nodes[g->cells[i].nodesInd[2]].z, g->nodes[g->cells[i].nodesInd[3]].z));
		g->cells[i].fCount = 4;
		g->cells[i].facesInd = new int[g->cells[i].fCount];
	}
	fclose(fp);

	// формируем данные о ГРАНЯХ
	sprintf(str, "%s.neigh", fileName);
	fp = fopen(str, "r");
	fscanf(fp, "%d %d", &tmp, &tmp);
	int** neigh;
	neigh = new int*[g->cCount];
	for (int i = 0; i < g->cCount; i++)
	{
		neigh[i] = new int[4];
		fscanf(fp, "%d %d %d %d %d", &tmp, &(neigh[i][0]), &(neigh[i][1]), &(neigh[i][2]), &(neigh[i][3]));
		neigh[i][0]--;
		neigh[i][1]--;
		neigh[i][2]--;
		neigh[i][3]--;
		g->cells[i].neigh[0] = neigh[i][0];
		g->cells[i].neigh[1] = neigh[i][1];
		g->cells[i].neigh[2] = neigh[i][2];
		g->cells[i].neigh[3] = neigh[i][3];
	}
	fclose(fp);
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
	int * cfi = new int[g->cCount];
	for (int i = 0; i < g->cCount; i++)
	{
		cfi[i] = 0;
	}
	// ::memset(cfi, 0, cCount*sizeof(int));
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
				//g->faces[iFace].c = new Point[g->faces[iFace].cCount];
				g->faces[iFace].c.x = (g->nodes[g->faces[iFace].nodesInd[0]].x + g->nodes[g->faces[iFace].nodesInd[1]].x + g->nodes[g->faces[iFace].nodesInd[2]].x) / 3.0;
				g->faces[iFace].c.y = (g->nodes[g->faces[iFace].nodesInd[0]].y + g->nodes[g->faces[iFace].nodesInd[1]].y + g->nodes[g->faces[iFace].nodesInd[2]].y) / 3.0;
				g->faces[iFace].c.z = (g->nodes[g->faces[iFace].nodesInd[0]].z + g->nodes[g->faces[iFace].nodesInd[1]].z + g->nodes[g->faces[iFace].nodesInd[2]].z) / 3.0;
				
				double x1 = g->nodes[g->faces[iFace].nodesInd[0]].x; double y1 = g->nodes[g->faces[iFace].nodesInd[0]].y; double z1 = g->nodes[g->faces[iFace].nodesInd[0]].z;
				double x2 = g->nodes[g->faces[iFace].nodesInd[1]].x; double y2 = g->nodes[g->faces[iFace].nodesInd[1]].y; double z2 = g->nodes[g->faces[iFace].nodesInd[1]].z;
				double x3 = g->nodes[g->faces[iFace].nodesInd[2]].x; double y3 = g->nodes[g->faces[iFace].nodesInd[2]].y; double z3 = g->nodes[g->faces[iFace].nodesInd[2]].z;

				double a1 = x1 - x3; double a2 = x2 - x3; double a3 = x3;
				double b1 = y1 - y3; double b2 = y2 - y3; double b3 = y3;
				double c1 = z1 - z3; double c2 = z2 - z3; double c3 = z3;
				
				g->faces[iFace].n.x = b1*c2 - c1*b2;
				g->faces[iFace].n.y = -(a1*c2 - c1*a2);
				g->faces[iFace].n.z = a1*b2 - b1*a2;
				g->faces[iFace].n.x /= sqrt((b1*c2 - c1*b2)*(b1*c2 - c1*b2) + (a1*c2 - c1*a2)*(a1*c2 - c1*a2) + (a1*b2 - b1*a2)*(a1*b2 - b1*a2));
				g->faces[iFace].n.y /= sqrt((b1*c2 - c1*b2)*(b1*c2 - c1*b2) + (a1*c2 - c1*a2)*(a1*c2 - c1*a2) + (a1*b2 - b1*a2)*(a1*b2 - b1*a2));
				g->faces[iFace].n.z /= sqrt((b1*c2 - c1*b2)*(b1*c2 - c1*b2) + (a1*c2 - c1*a2)*(a1*c2 - c1*a2) + (a1*b2 - b1*a2)*(a1*b2 - b1*a2));

				double x4 = g->nodes[g->cells[i].nodesInd[j]].x; double y4 = g->nodes[g->cells[i].nodesInd[j]].y; double z4 = g->nodes[g->cells[i].nodesInd[j]].z;
				if ((x4 - x3)*g->faces[iFace].n.x + (y4 - y3)*g->faces[iFace].n.y + (z4 - z3)*g->faces[iFace].n.z > 0)
				{
					g->faces[iFace].n.x *= -1;
					g->faces[iFace].n.y *= -1;
					g->faces[iFace].n.z *= -1;
				}
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
				}
				iFace++;
			}
		}
	}

	// чтение данных о граничных гранях
	sprintf(str, "%s.face", fileName);
	fp = fopen(str, "r");
	int bndCount;
	fscanf(fp, "%d %d", &bndCount, &tmp);
	for (int i = 0; i < bndCount; i++)
	{
		int n, n1, n2, n3, type;
		fscanf(fp, "%d %d %d %d %d", &n, &n1, &n2, &n3, &type);
		n1--;
		n2--;
		n3--;
		int iFace = g->findFace(n1, n2, n3);
		if (iFace >= 0) {
			g->faces[iFace].type = type;
			if (g->faces[iFace].type != 0) {
				sprintf(g->faces[iFace].typeName, "%d", g->faces[iFace].type);
			}
		}
	}
	fclose(fp);
		
	//// формируем данные о РЕБРАХ
	//std::map<std::pair<int, int>, std::set<int>> pair_faces;
	//for (int i = 0; i < g->cCount; i++)
	//{
	//	Cell& cell = g->cells[i];
	//	for (int j = 0; j < 4; j++)
	//	{
	//		Face& face = g->faces[cell.facesInd[j]];
	//		for (int k = 0; k < 3; k++) {
	//			int n1 = _min_(face.nodesInd[k], face.nodesInd[(k + 1) % 3]);
	//			int n2 = _max_(face.nodesInd[k], face.nodesInd[(k + 1) % 3]);
	//			pair_faces[std::pair<int, int>(n1, n2)].insert(cell.facesInd[j]);
	//		}
	//	}
	//}
	//g->eCount = pair_faces.size();
	//g->edges = new Edge[g->eCount];

	//int iEdge = 0;
	//for (std::map<std::pair<int, int>, std::set<int>>::iterator it = pair_faces.begin(); it != pair_faces.end(); ++it) {
	//	g->edges[iEdge].n1 = it->first.first;
	//	g->edges[iEdge].n2 = it->first.second;
	//	for (std::set<int>::iterator it_set = it->second.begin(); it_set != it->second.end(); ++it_set) {
	//		for (int j = 0; j < 3; j++) {
	//			int n1 = _min_(g->faces[*it_set].nodesInd[j], g->faces[*it_set].nodesInd[(j + 1) % 3]);
	//			int n2 = _max_(g->faces[*it_set].nodesInd[j], g->faces[*it_set].nodesInd[(j + 1) % 3]);
	//			if (n1 == g->edges[iEdge].n1 && n2 == g->edges[iEdge].n2) {
	//				g->faces[*it_set].edgesInd[(j + 2) % 3] = iEdge;
	//			}
	//		}
	//	}
	//	iEdge++;
	//}
	/*for (int i = 0; i < g->cCount; i++)
	{
		double a = g->edges[g->cells[i].edgesInd[0]].l;
		double b = g->edges[g->cells[i].edgesInd[1]].l;
		double c = g->edges[g->cells[i].edgesInd[2]].l;
		double p = (a + b + c) / 2.0;
		g->cells[i].S = sqrt(p*(p - a)*(p - b)*(p - c));
	}*/

	for (int i = 0; i < g->cCount; i++)
	{
		delete[] neigh[i];
	}
	delete[] neigh;
	delete[] cfi;


}


