
/**
* @file eigenPairs.h
* @author Melih Altun @2015-2017
**/

#if !defined (_EIGENPAIRS_H_)
#define _EIGENPAIRS_H_

#include "Hessenberg.h"
#include "implicitQR.h"
#include "GivensQR.h"
#include "backSubstitution.h"

void eig(complex A[], unsigned int N, complex V[], complex D[]);

#endif
