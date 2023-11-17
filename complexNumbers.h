
/**
* @file complexNumbers.h
* @author Melih Altun @2015-2017
**/

#if !defined (_COMPLEX_NUMS_H_)
#define _COMPLEX_NUMS_H_
#include <math.h>

typedef struct complex_ {
	double real;
	double imag;
}complex;

complex complexAddition(complex z1, complex z2);
complex complexSubtraction(complex z1, complex z2);
complex complexMultiplication(complex z1, complex z2);
complex complexDivision(complex z1, complex z2);
complex complexConjugate(complex z1);
complex negate(complex z1);
complex equals(complex z1);
double magnitude(complex z1);
double angle(complex z1);
complex complexExp(double alpha);
void z2polar(complex z, double *r, double *theta);
complex polar2z(double r, double theta);
complex complexSqrt(complex z1);

#endif