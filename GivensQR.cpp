
/**
* @file GivensQR.cpp
* @author Melih Altun @2015-2017
**/

#include "GivensQR.h"
#define EPS 2.221E-12

// plane rotation using Givens' method
void complexGivensRotation(complex a, complex b, complex *c, complex *s)
{
	double tolr = EPS, r;
	if (magnitude(b) < tolr) {
		*c = { 1.0, 0.0 };
		*s = { 0.0, 0.0 };
	}
	else if (magnitude(a) < tolr) {
		*c = { 0.0, 0.0 };
		*s = { -1.0, 0.0 };
	}
	else {
		r = sqrt(magnitude(a)*magnitude(a) + magnitude(b)*magnitude(b));
		*c = complexDivision(complexConjugate(a), { r, 0.0 });
		*s = complexDivision(complexConjugate(b), { r, 0.0 });
	}
}


// QR decomposition => [Q,R] = qr(A), where: Q*R = A
void qrGivensComplex(complex A[], unsigned int N, unsigned int M, complex Q[], complex R[])
{
	unsigned int i, j;

	complex c, s;
	complex *G, *Gt, *Rnew, *Qnew;
	G = (complex*)malloc(N*N*sizeof(complex));
	Gt = (complex*)malloc(N*N * sizeof(complex));
	Rnew = (complex*)malloc(N*M * sizeof(complex));
	Qnew = (complex*)malloc(N*N * sizeof(complex));
	
	complexIdentity(Q, N);
	copyComplexMatrix(A, N, M, R);

	for (j = 0; j < M; j++) {
		for (i = N - 1; i > j; i--) {
			complexIdentity(G, N);
			complexGivensRotation(R[lin_index((i - 1), j, M)], R[lin_index(i, j, M)], &c, &s);
			G[lin_index((i - 1), (i - 1), N)] = c;
			G[lin_index((i - 1), i, N)] = s;
			G[lin_index(i, (i - 1), N)] = negate(complexConjugate(s));
			G[lin_index(i, i, N)] = complexConjugate(c);

			//R = G*R
			//Q = Q*G'
			hermitian(G, N, N, Gt);
			complexMatrixMultiplication(G, R, N, N, M, Rnew);
			complexMatrixMultiplication(Q, Gt, N, N, N, Qnew);
			copyComplexMatrix(Rnew, N, M, R);
			copyComplexMatrix(Qnew, N, N, Q);
		}
	}
	delete[] G;
	delete[] Gt;
	delete[] Rnew;
	delete[] Qnew;
}


//QR algortihm for finding eigen values and eigen vectors
void qrIteration(complex A[], unsigned int N, complex Q[], complex R[])
{
	unsigned int k,  maxIter = 1000;
	double n1, n2, tolr = EPS;
	complex *M, *V, *D, *Vnew;

	M = (complex*)malloc(N*N * sizeof(complex));
	V = (complex*)malloc(N*N * sizeof(complex));
	D = (complex*)malloc(N * sizeof(complex));
	Vnew = (complex*)malloc(N*N * sizeof(complex));

	memcpy(M, A, N*N * sizeof(complex));
	complexIdentity(V, N);
	
	complexDiagonal(M, N, D);
	n1 = complexVectorNorm(D,N);

	for (k = 0; k < maxIter; k++) {
		//[Q,R] = qr(M)
		qrGivensComplex(M, N, N, Q, R);

		//V = V*Q;
		complexMatrixMultiplication(V, Q, N, N, N, Vnew);
		copyComplexMatrix(Vnew, N, N, V);
		//M = R*Q;
		complexMatrixMultiplication(R, Q, N, N, N, M);

		complexDiagonal(M, N, D);
		n2 = complexVectorNorm(D, N);

		if (fabs(n1 - n2) < tolr)
			break;

		n1 = n2;
	}
	copyComplexMatrix(V, N, N, Q);
	copyComplexMatrix(M, N, N, R);

	//clear memory
	delete[] M;
	delete[] V;
	delete[] Vnew;
	delete[] D;
}
