#pragma once
#include "MeshReader.h"
#include <map>
#include <set>

class MeshReaderGmsh : 
	public MeshReader
{
public:
	MeshReaderGmsh(char* fName): fileName(fName) {}
	~MeshReaderGmsh() { delete[] fileName; }

	typedef std::vector<std::string> string_list;
	typedef std::vector<int> indexes;
	typedef std::vector<indexes> index_list;
	typedef std::map<std::string, indexes> bnd_map;
	typedef std::set<int> ind_set;
	
	virtual void read(Grid*);

private:
	char* fileName;
	index_list faces;
	index_list cells;
	index_list bnd_faces;

	int findFace(int n1, int n2, int n3);
	std::vector<int> getIntersection(index_list &sets);
};

