#pragma once

#include "fem_dg.h"

class Limiter
{
public:
	virtual ~Limiter();
	virtual void run() = 0;
	virtual const char* getName() = 0;
	static Limiter* create(const char* limiterName, FEM_DG* solver);
};

