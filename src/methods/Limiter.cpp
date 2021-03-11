#include "Limiter.h"
#include <cstring>
#include "LimiterCockburn.h"
#include "LimiterBJ.h"
#include "LimiterEntropy.h"

Limiter::~Limiter()
{
}

Limiter* Limiter::create(const char* limiterName, FEM_DG* solver)
{
	if (strcmp(limiterName, "COCKBURN") == 0) {
		return new LimiterCockburn(solver, &(solver->grid), solver->ro, solver->ru, solver->rv, solver->rw, solver->re, solver->BASE_FUNC_COUNT);
	}
	else if (strcmp(limiterName, "BARTH_JESPERSEN") == 0) {
		return new LimiterBJ(solver, &(solver->grid), solver->ro, solver->ru, solver->rv, solver->rw, solver->re, solver->BASE_FUNC_COUNT);
	}
	else if (strcmp(limiterName, "ENTROPY") == 0) {
		return new LimiterEntropy(solver, &(solver->grid), solver->ro, solver->ru, solver->rv, solver->rw, solver->re, solver->BASE_FUNC_COUNT, &(solver->limParam));
	}
	else {
		return NULL;
	}
}
