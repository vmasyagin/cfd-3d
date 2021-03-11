#pragma once
#include "MeshReader.h"
class MeshReaderWIASTetGen :
	public MeshReader
{
private:
	char* fileName;
public:
	MeshReaderWIASTetGen(char* fName): fileName(fName) {}
	~MeshReaderWIASTetGen() { delete[] fileName; }
	
	virtual void read(Grid*);
};

