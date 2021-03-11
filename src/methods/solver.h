#ifndef _SOLVER_H_
#define _SOLVER_H_

#include <method.h>


class Solver
{
public:
	static Method *	initMethod		( char* fileName );
	static void		runMethod		( Method* m );
	static void		destroyMethod	( Method* m );
};

#endif