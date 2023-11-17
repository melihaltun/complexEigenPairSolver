
/**
* @file backSubstitution.cpp
* @author Melih Altun @2015-2017
**/

#include "backSubstitution.h"
#define EPS 2.221E-12

// Finds x in linear system: U x = y, where x, and y are vectors and U is a N x N upper triangular matrix
void backwardsSubstitution(complex U[], complex y[], unsigned int N, complex x[])
{
	int MSD, i, j;
	double absMax_U, tolr;

	absMax_U = absMaxElement(U, N, N);

	MSD = (int)ceil(log10(absMax_U));
	tolr = EPS * pow(10.0, (double)(MSD + 1));

	for (i = N-1; i >= 0; i--) {
		x[i] = equals(y[i]);
		for (j = i + 1; j < N; j++)
			x[i] = equals(complexSubtraction(x[i], complexMultiplication(U[lin_index(i, j, N)], x[j])));  //x(i) = x(i)-U(i,j)*x(j);
		if (magnitude(U[lin_index(i, i, N)]) < tolr) {
			x[i].real = 1; x[i].imag = 0;
		}
		else
			x[i] = complexDivision(x[i], U[lin_index(i, i, N)]);
	}
}
