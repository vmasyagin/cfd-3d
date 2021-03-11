#include <stdlib.h>
#include <stdio.h>
#include <solver.h>
#include "global.h"
#include "mpi.h"
#include <float.h>
#include <iostream>
#include <algorithm>
#include <fstream>

#include "SolverZeidel.h"

int main(int argc, char** argv)
{

#ifdef _DEBUG
	_controlfp(~(_MCW_EM & (~_EM_INEXACT) & (~_EM_UNDERFLOW)), _MCW_EM);
#endif

	Parallel::init(&argc, &argv);

	hLog = fopen("task.log", "w"); // открываем файл для записи лога; вместо printf(...) необходимо использовать log(...)
	
	Method * method = Solver::initMethod( "task.xml" ); 
	
	Solver :: runMethod		( method );
	Solver :: destroyMethod	( method );

	fclose(hLog);

	Parallel::done();

	return 0;
}