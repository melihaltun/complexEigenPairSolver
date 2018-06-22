
/**
* @file complexMatrixOperations.cpp
* @author Melih Altun @2015-2017
**/

#include "complexMatrixOperations.h"

// out[] = in[]
void copyComplexMatrix(complex in[], unsigned int N, unsigned int M, complex out[])
{
	int i, j;
	for (i = 0; i < N; i++)
		for (j = 0; j < M; j++)
			out[lin_index(i, j, M)] = equals(in[lin_index(i, j, M)]);

}


// out[] = (complex)in[]
void double2complex(double in[], unsigned int M, unsigned int N, complex out[])
{
	int i, j;
	for (i = 0; i < N; i++)
		for (j = 0; j < M; j++) {
			out[lin_index(i, j, M)].real = in[lin_index(i, j, M)];
			out[lin_index(i, j, M)].imag = 0.0;
		}
}


// out[] = conj(tran(in[]))
void hermitian(complex in[], unsigned int N, unsigned int M, complex out[])
{
	unsigned int i, j;

		for (i = 0; i<N; i++)
			for (j = 0; j < M; j++)
				out[lin_index(j, i, N)] = complexConjugate(in[lin_index(i, j, M)]);
}


// out[] = tran(in[])
void complexTranspose(complex in[], unsigned int N, unsigned int M, complex out[])
{
	unsigned int i, j;

	for (i = 0; i<N; i++)
		for (j = 0; j < M; j++)
			out[lin_index(j, i, N)] = equals(in[lin_index(i, j, M)]);
}


// out[] = I
void complexIdentity(complex matrix[], int size)
{
	int i;
	memset(matrix, 0, size*size * sizeof(complex));
	for (i = 0; i < size; i++)
		matrix[lin_index(i, i, size)].real = 1;
}


// x = max(abs(in[]))
double absMaxElement(complex in[], unsigned int N, unsigned int M)
{
	double absMax_in, absElement;
	int i;

	absMax_in = magnitude(in[0]);
	for (i = 1; i < M*N; i++) {
		absElement = magnitude(in[i]);
		if (absMax_in < absElement)
			absMax_in = absElement;
	}

	return absMax_in;
}


// x = norm(V[]);
double complexVectorNorm(complex V[], unsigned int N)
{
	double sum = 0;
	int i;
	
	for (i = 0; i < N; i++)
		sum += pow(magnitude(V[i]), 2);

	return sqrt(sum);
}


// U = V1 - V2
void complexVectorSubtraction(complex V1[], complex V2[], unsigned int N, complex U[])
{
	int i;

	for (i = 0; i < N; i++)
		U[i] = complexSubtraction(V1[i],V2[i]);
}


// U = a V;
void complexScalarMultiplication(complex V[], complex alpha, unsigned int N, unsigned int M, complex U[])
{
	int i;

	for (i = 0; i < N*M; i++)
		U[i] = complexMultiplication(alpha, V[i]);
}


// Mout = M1 x M2
void complexMatrixMultiplication(complex M1[], complex M2[], unsigned int N, unsigned int K, unsigned int M, complex Mout[])
{
	int i, j, k;

	memset(Mout, 0, N*M * sizeof(complex));

	for (i = 0; i<N; i++)
		for (j = 0; j<M; j++)
			for (k = 0; k<K; k++)
				Mout[(i*M) + j] = complexAddition(Mout[(i*M) + j], complexMultiplication(M1[(i*K) + k], M2[(k*M) + j]));
}


// Mout = M1 - M2
void complexMatrixSubtraction(complex M1[], complex M2[], unsigned int N, unsigned int M, complex Mout[])
{
	int i, j;

	memset(Mout, 0, N*M * sizeof(complex));

	for (i = 0; i < N; i++)
		for (j = 0; j < M; j++)
			Mout[(i*M) + j] = complexSubtraction(M1[(i*M) + j], M2[(i*M) + j]);
}


//out[] = conj(in[])
void complexMatrixConjungate(complex in[], unsigned int N, unsigned int M, complex out[])
{
	int i;

	for (i = 0; i < N*M; i++)
		out[i] = complexConjugate(in[i]);
}


//D = diag(M[])
void complexDiagonal(complex M[], unsigned int N, complex D[])
{
	int i;

	for (i = 0; i < N; i++)
		D[i] = M[lin_index(i,i,N)];
}


//M[] = diag(D)
void complexDiagonal2Vector(complex D[], unsigned int N, complex M[])
{
	int i;

	memset(M, 0, N*N * sizeof(complex));
	for (i = 0; i < N; i++)
		M[lin_index(i, i, N)] = D[i];
}
