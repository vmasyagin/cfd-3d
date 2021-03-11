#ifndef _METHOD_H_
#define _METHOD_H_

#include "grid.h"
#include "global.h"
class Method
{
public:
	virtual void init(char * xmlFileName) = 0;
	virtual void run() = 0;
	virtual void done() = 0;

	void exchange(double* field)
	{
		for (int p = 0; p < Parallel::procCount; p++) {
			if (p < Parallel::procId) {
				if (grid.recvCount[p] > 0) {
					Parallel::recv(p, 0, grid.recvCount[p], &field[grid.recvShift[p]]);
				}
				int n = grid.sendInd[p].size();
				if (n > 0) {
					for (int i = 0; i < n; i++) {
						dBuf[i] = field[grid.sendInd[p][i]];
					}
					Parallel::send(p, 1, n, dBuf);
				}
			}
			else if (p > Parallel::procId) {
				int n = grid.sendInd[p].size();
				if (n > 0) {
					for (int i = 0; i < n; i++) {
						dBuf[i] = field[grid.sendInd[p][i]];
					}
					Parallel::send(p, 0, n, dBuf);
				}
				if (grid.recvCount[p] > 0) {
					Parallel::recv(p, 1, grid.recvCount[p], &field[grid.recvShift[p]]);
				}
			}
		}
	}

	void exchange(int* field)
	{
		char fName[64];
		sprintf(fName, "exchange.%04d.txt", Parallel::procId);
		FILE * fp = fopen(fName, "w");
		fprintf(fp, "COUNT:    %d:\n", grid.cCount);
		fprintf(fp, "COUNT_EX: %d:\n", grid.cCountEx);
		for (int p = 0; p < Parallel::procCount; p++) {
			if (p < Parallel::procId) {
				if (grid.recvCount[p] > 0) {
					Parallel::recv(p, 0, grid.recvCount[p], &field[grid.recvShift[p]]);
					fprintf(fp, "RECV from %d:\n", p);
					for (int i = 0; i < grid.recvCount[p]; i++) {
						fprintf(fp, "  %d: %d\n", grid.recvShift[p] + i, field[grid.recvShift[p] + i]);
					}
					fprintf(fp, "--------------------------\n");
				}
				int n = grid.sendInd[p].size();
				if (n > 0) {
					for (int i = 0; i < n; i++) {
						iBuf[i] = field[grid.sendInd[p][i]];
					}
					Parallel::send(p, 1, n, iBuf);
					fprintf(fp, "SENT to %d:\n", p);
					for (int i = 0; i < n; i++) {
						fprintf(fp, "  %d: %d\n", grid.sendInd[p][i], field[grid.sendInd[p][i]]);
					}
					fprintf(fp, "--------------------------\n");
				}
			}
			else if (p > Parallel::procId) {
				int n = grid.sendInd[p].size();
				if (n > 0) {
					for (int i = 0; i < n; i++) {
						iBuf[i] = field[grid.sendInd[p][i]];
					}
					Parallel::send(p, 0, n, iBuf);
					fprintf(fp, "SENT to %d:\n", p);
					for (int i = 0; i < n; i++) {
						fprintf(fp, "  %d: %d\n", grid.sendInd[p][i], field[grid.sendInd[p][i]]);
					}
					fprintf(fp, "--------------------------\n");
				}
				if (grid.recvCount[p] > 0) {
					Parallel::recv(p, 1, grid.recvCount[p], &field[grid.recvShift[p]]);
					fprintf(fp, "RECV from %d:\n", p);
					for (int i = 0; i < grid.recvCount[p]; i++) {
						fprintf(fp, "  %d: %d\n", grid.recvShift[p] + i, field[grid.recvShift[p] + i]);
					}
					fprintf(fp, "--------------------------\n");
				}
			}
		}
		fclose(fp);
	}

	void exchange(VECTOR* field)
	{
		for (int p = 0; p < Parallel::procCount; p++) {
			if (p < Parallel::procId) {
				if (grid.recvCount[p] > 0) {
					Parallel::recv(p, 0, grid.recvCount[p], &field[grid.recvShift[p]]);
				}
				int n = grid.sendInd[p].size();
				if (n > 0) {
					for (int i = 0; i < n; i++) {
						vBuf[i] = field[grid.sendInd[p][i]];
					}
					Parallel::send(p, 1, n, vBuf);
				}
			}
			else if (Parallel::procId < p) {
				int n = grid.sendInd[p].size();
				if (n > 0) {
					for (int i = 0; i < n; i++) {
						vBuf[i] = field[grid.sendInd[p][i]];
					}
					Parallel::send(p, 0, n, vBuf);
				}
				if (grid.recvCount[p] > 0) {
					Parallel::recv(p, 1, grid.recvCount[p], &field[grid.recvShift[p]]);
				}
			}
		}
	}

protected:
	Grid grid;

	double * dBuf;
	int    * iBuf;
	VECTOR * vBuf;
};

#endif