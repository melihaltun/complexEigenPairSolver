
/**
* @file Hessenberg.cpp
* @author Melih Altun @2015-2017
**/

#include "Hessenberg.h"
#define EPS 2.221E-12

//Hessenberg reduction of matrix A. A = U*H where H is an upper Hessenberg matrix.
void hessenberg(complex A[], unsigned int N, complex H[], complex U[])
{
	double absMax_A, tolr, eps = 2.221E-16, norm_x_alpha_e1;
	int MSD, k, j, i;
	complex *x, *e1, *u, *u_conj, *x_alpha_e1, *P, *u_conj_u, *Ph, *Ph_H, *Ph_H_P, *U_P;
	complex ro, alpha;


	copyComplexMatrix(A, N, N, H);  //H=A;
	absMax_A = absMaxElement(A, N, N);


	MSD = (int)ceil(log10(absMax_A));
	tolr = EPS * pow(10.0, (double)(MSD + 1));

	complexIdentity(U, N);

	x = (complex*)malloc((N-1) * sizeof(complex));
	e1 = (complex*)malloc((N - 1) * sizeof(complex));
	u = (complex*)malloc((N - 1) * sizeof(complex));
	u_conj = (complex*)malloc((N - 1) * sizeof(complex));
	x_alpha_e1 = (complex*)malloc((N - 1) * sizeof(complex));
	P = (complex*)malloc(N * N * sizeof(complex));
	u_conj_u = (complex*)malloc((N - 1) * (N - 1) * sizeof(complex));
	Ph = (complex*)malloc(N * N * sizeof(complex));
	Ph_H = (complex*)malloc(N * N * sizeof(complex));
	Ph_H_P = (complex*)malloc(N * N * sizeof(complex));
	U_P = (complex*)malloc(N * N * sizeof(complex));

	for (k = 0; k < N - 2; k++) {
		for (j = k + 1; j < N; j++)
			x[j - k - 1] = equals(H[lin_index(j, k, N)]);

		if (magnitude(x[0]) < tolr)
			ro = equals({ -1,0 });
		else
			ro = negate(complexExp(angle(x[0])));   //ro = -(exp(j*atan2(imag(x(1)),real(x(1)))))

		alpha = complexMultiplication(ro, { complexVectorNorm(x, N - 1 - k), 0 });

		memset(e1, 0, (N - k - 1) * sizeof(complex));
		e1[0] = equals({ 1, 0 });

		// u = (x - alpha*e1) / norm(x - alpha*e1)
		complexScalarMultiplication(e1, alpha, N - k - 1, 1, x_alpha_e1);
		complexVectorSubtraction(x, x_alpha_e1, N - k - 1, x_alpha_e1);
		norm_x_alpha_e1 = complexVectorNorm(x_alpha_e1, N - k - 1);
		complexScalarMultiplication(x_alpha_e1, { 1 / norm_x_alpha_e1 , 0 }, N-k-1, 1, u);

		complexIdentity(P, N);  //P = I

		// P(k+1:N,k+1:N) = P(k+1:N,k+1:N) - 2*u*(conj(u))';
		complexMatrixConjungate(u, N - k - 1, 1, u_conj);
		complexMatrixMultiplication(u, u_conj, N - k - 1, 1, N - k - 1, u_conj_u);
		complexScalarMultiplication(u_conj_u, { 2, 0 }, N - k - 1, N - k - 1, u_conj_u);
		for (i = k + 1; i < N; i++)
			for (j = k + 1; j < N; j++)
				P[lin_index(i, j, N)] = complexSubtraction(P[lin_index(i, j, N)], u_conj_u[lin_index(i-k-1, j-k-1, N-k-1)]);

		// P = P'*H*P;
		hermitian(P, N, N, Ph);
		complexMatrixMultiplication(Ph, H, N, N, N, Ph_H);
		complexMatrixMultiplication(Ph_H, P, N, N, N, Ph_H_P);
		memcpy(H, Ph_H_P, N*N * sizeof(complex));

		// U = U*P;
		complexMatrixMultiplication(U, P, N, N, N, U_P);
		memcpy(U, U_P, N*N * sizeof(complex));
	}

	//clear memory
	delete[] x;
	delete[] e1;
	delete[] u;
	delete[] x_alpha_e1;
	delete[] P;
	delete[] u_conj;
	delete[] u_conj_u;
	delete[] Ph;
	delete[] Ph_H;
	delete[] Ph_H_P;
	delete[] U_P;
}
