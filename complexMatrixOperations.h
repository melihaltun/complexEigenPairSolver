
/**
* @file complexMatrixOperations.h
* @author Melih Altun @2015-2017
**/

#if !defined (_C_MATRIX_OPS_H_)
#define _C_MATRIX_OPS_H_

#include "complexNumbers.h"
#include <string.h>
#include <stdlib.h>

#define lin_index(i, j, numCol)  ( ((i)*(numCol))+(j) )   //2D to 1D array

void copyComplexMatrix(complex in[], unsigned int N, unsigned int M, complex out[]);

void double2complex (double in[], unsigned int N, unsigned int M, complex out[]);

void hermitian(complex in[], unsigned int N, unsigned int M, complex out[]);

void complexTranspose(complex in[], unsigned int N, unsigned int M, complex out[]);

void complexIdentity(complex matrix[], int size);

double absMaxElement(complex in[], unsigned int N, unsigned int M);

double complexVectorNorm(complex V[], unsigned int N);

void complexVectorSubtraction(complex V1[], complex V2[], unsigned int N, complex U[]);

void complexScalarMultiplication(complex V[], complex alpha, unsigned int N, unsigned int M, complex U[]);

void complexMatrixMultiplication(complex M1[], complex M2[], unsigned int N, unsigned int K, unsigned int M, complex Mout[]);

void complexMatrixSubtraction(complex M1[], complex M2[], unsigned int N, unsigned int M, complex Mout[]);

void complexMatrixConjungate(complex in[], unsigned int N, unsigned int M, complex out[]);

void complexDiagonal(complex M[], unsigned int N, complex D[]);

void complexDiagonal2Vector(complex D[], unsigned int N, complex M[]);

#endif
