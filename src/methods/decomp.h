#include "method.h"


class Decomp : public Method
{
public:
	virtual void init(char * xmlFileName);
	virtual void run();
	virtual void done();
protected:
	Grid * grids;
	int procCount;
};