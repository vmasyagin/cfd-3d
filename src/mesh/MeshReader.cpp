#include "MeshReader.h"
#include "MeshReaderWIASTetGen.h"
//#include "MeshReaderSalomeUnv.h"
#include "MeshReaderGmsh.h"



MeshReader* MeshReader::create(int type, char* fileName)
{
	switch (type) {
	case TYPE_WIAS_TET:
		return new MeshReaderWIASTetGen(fileName);
		break;
	//case TYPE_SALOME:
	//	return new MeshReaderSalomeUnv(fileName);
	//	break;
	case TYPE_GMSH:
		return new MeshReaderGmsh(fileName);
		break;	
	}
}

int MeshReader::getType(char* name) {
	if (strcmp(name, "wias_tetgen") == 0) return MeshReader::TYPE_WIAS_TET;
	//else if (strcmp(name, "salome_unv") == 0) return MeshReader::TYPE_SALOME;
	else if (strcmp(name, "gmsh") == 0) return MeshReader::TYPE_GMSH;
	throw Exception((char*)"Wrong mesh type.", Exception::TYPE_MESH_WRONG_NAME);
}
