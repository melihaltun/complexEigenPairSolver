
/**
* @file GivensQR.h
* @author Melih Altun @2015-2017
**/

#if !defined (_QR_GIVENS_H_)
#define _QR_GIVENS_H_

#include "complexMatrixOperations.h"

void complexGivensRotation(complex a, complex b, complex *c, complex *s);

void qrGivensComplex(complex A[], unsigned int N, unsigned int M, complex Q[], complex R[]);

void qrIteration(complex A[], unsigned int N, complex Q[], complex R[]);

#endif
