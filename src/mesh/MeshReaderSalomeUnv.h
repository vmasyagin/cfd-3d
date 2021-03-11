#pragma once
#include "MeshReader.h"
#include <vector>
#include <string>
#include <map>
#include <set>

using namespace std;



class MeshReaderSalomeUnv :
	public MeshReader
{
public:
	struct element {
		int ind;
		int type;

		static const int TYPE_EDGE = 11;
		static const int TYPE_CELL = 41;

		element() : ind(0), type(0){}
		element(int i, int t) : ind(i), type(t) {}
	};

	typedef vector<string> string_list;
	typedef vector<int> indexes;
	typedef vector<indexes> index_list;
	typedef map<string, indexes> bnd_map;
	typedef set<int> ind_set;
	typedef vector<element> ele_map;

private:
	char* fileName;
	index_list edges, cells;
	bnd_map bounds;
	vector<Point> points;
	ele_map elements;

	void read_block(string_list * sl, ifstream& fin);
	void parse_block(string_list sl, Grid * g);
	void parse_block_164(string_list sl, Grid * g);
	void parse_block_2411(string_list sl, Grid * g);
	void parse_block_2412(string_list sl, Grid * g);
	void parse_block_2467(string_list sl, Grid * g);
	void parse_block_2420(string_list sl, Grid * g);

	void edge_add(int i1, int i2);
	//bool edge_exist(int i1, int i2);
	int find_index(int val);
	int find_edge(int n1, int n2);
public:
	MeshReaderSalomeUnv(char* fName) : fileName(fName) {}
	~MeshReaderSalomeUnv() { delete[] fileName; }

	virtual void read(Grid*);
};

