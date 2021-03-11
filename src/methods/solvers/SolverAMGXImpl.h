#pragma once

#include "SolverAMGX.h"


class SolverAMGXImpl : public SolverAMGX
{
public:

	SolverAMGXImpl(char* name) {
		solver_name = name;
	}
	virtual int solve(double eps, int& maxIter);
	virtual char* getName() { return solver_name; }
};

