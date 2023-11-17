
/**
* @file complexNumbers.cpp
* @author Melih Altun @2015-2017
**/

#include "complexNumbers.h"


// z = z1 + z2
complex complexAddition(complex z1, complex z2)
{
	complex z;
	z.real = z1.real + z2.real;
	z.imag = z1.imag + z2.imag;
	return z;
}

// z = z1 - z2
complex complexSubtraction(complex z1, complex z2)
{
	complex z;
	z.real = z1.real - z2.real;
	z.imag = z1.imag - z2.imag;
	return z;
}

// z = z1 x z2
complex complexMultiplication(complex z1, complex z2)
{
	complex z;
	z.real = (z1.real * z2.real) - (z1.imag * z2.imag);
	z.imag = (z1.real * z2.imag) + (z1.imag * z2.real);
	return z;
}

// z = z1 / z2
complex complexDivision(complex z1, complex z2)
{
	complex z, conj_z2, z2_conj_z2;

	conj_z2 = complexConjugate(z2);
	z2_conj_z2 = complexMultiplication(z2, conj_z2);
	z = complexMultiplication(z1, conj_z2);
	z.real /= z2_conj_z2.real;
	z.imag /= z2_conj_z2.real;
	return z;
}

// z = z1*
complex complexConjugate(complex z1)
{
	complex z;
	z.real = z1.real;
	z.imag = -z1.imag;
	return z;
}

// z = -z1
complex negate(complex z1)
{
	complex z;
	z.real = -z1.real;
	z.imag = -z1.imag;
	return z;
}

// z = z1
complex equals(complex z1)
{
	complex z;
	z.real = z1.real;
	z.imag = z1.imag;
	return z;
}

// m = |z|
double magnitude(complex z1)
{
	return sqrt(z1.real*z1.real + z1.imag*z1.imag);
}

// theta = <z
double angle(complex z1)
{
	return atan2(z1.imag, z1.real);
}

// z = exp(j alpha);
complex complexExp(double alpha)
{
	complex z;
	z.real = cos(alpha);
	z.imag = sin(alpha);
	return z;
}

// r = |m|, theta = <z
void z2polar(complex z, double *r, double *theta)
{
	*r = magnitude(z);
	*theta = angle(z);
}

// z = r exp(j theta)
complex polar2z(double r, double theta)
{
	complex z;
	z.real = r*cos(theta);
	z.imag = r*sin(theta);
	return z;
}

// z2 = sqrt(z1)
complex complexSqrt(complex z1) {
	double r, theta;
	complex z2;
	r = sqrt(z1.real * z1.real + z1.imag * z1.imag);
	theta = angle(z1);
	z2.real = sqrt(r) * cos(theta / 2);
	z2.imag = sqrt(r) * sin(theta / 2);
	return z2;
}
