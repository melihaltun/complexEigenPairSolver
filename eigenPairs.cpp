
/**
* @file eigenPairs.cpp
* @author Melih Altun @2015-2017
**/

#include"eigenPairs.h"


//V, D = eig(A)
void solve2by2eig(complex A[], complex V[], complex D[])
{
	double a, b, c, delta, tolr = 2.3E-16, r, theta, nV1, nV2, m1, m2, a1, a2;
	complex x, y, z, x2, V1[2], V2[2];

	memset(D, 0, 2 * sizeof(complex));

	a = 1;
	b = -A[lin_index(0, 0, 2)].real - A[lin_index(1, 1, 2)].real;
	c = A[lin_index(0, 0, 2)].real * A[lin_index(1, 1, 2)].real - A[lin_index(1, 0, 2)].real * A[lin_index(0, 1, 2)].real;

	// use discriminant to solve for eigenvalues
	delta = b*b - 4 * a*c;

	if (delta < 0) {
		D[0].real = -b / (2 * a);
		D[0].imag = sqrt(-delta) / (2 * a);
		D[1].real = -b / (2 * a);
		D[1].imag = -sqrt(-delta) / (2 * a);
	}
	else {
		D[0].real = (-b + sqrt(delta)) / (2 * a);
		D[1].real = (-b - sqrt(delta)) / (2 * a);
	}

	// find eigen vectors
	if (magnitude(A[lin_index(0, 1, 2)]) > tolr)
		//x = (D(1, 1) - A(1, 1)) / A(1, 2);
		x = complexDivision(complexSubtraction(D[0], A[lin_index(0, 0, 2)]), A[lin_index(0, 1, 2)]);
	else
		//x = D(1, 1) - A(1, 1);
		x = complexSubtraction(D[0], A[lin_index(0, 0, 2)]);

	//z = 1+x^2;
	x2 = complexMultiplication(x, x);
	z = complexAddition({ 1, 0 }, x2);

	if (magnitude(z) > 10 * tolr) {
		//y = sqrt(z);
		z2polar(z, &r, &theta);
		y = polar2z(sqrt(r), theta / 2);
	}
	else
		y = { 1, 0 };

	//V1 = [1/y; x/y];
	//V1 = V1/norm(V1);
	V1[0] = complexDivision({ 1,0 }, y);
	V1[1] = complexDivision(x, y);

	nV1 = complexVectorNorm(V1, 2);

	if (nV1 > tolr) {
		V1[0] = complexDivision(V1[0], { nV1, 0.0 });
		V1[1] = complexDivision(V1[1], { nV1, 0.0 });
	}

	if (fabs(V1[0].imag) > tolr) {
		m1 = magnitude(V1[0]);
		m2 = magnitude(V1[1]);
		a1 = angle(V1[0]);
		a2 = angle(V1[1]);

		V1[0] = { m1, 0 };
		V1[1] = { m2*cos(a2 - a1),  m2*sin(a2 - a1) };
	}

	// V1 found, now find V2

	if (magnitude(A[lin_index(0, 1, 2)]) > tolr)
		//x = (D(1, 1) - A(1, 1)) / A(1, 2);
		x = complexDivision(complexSubtraction(D[1], A[lin_index(0, 0, 2)]), A[lin_index(0, 1, 2)]);
	else
		//x = D(1, 1) - A(1, 1);
		x = complexSubtraction(D[1], A[lin_index(0, 0, 2)]);

	//z = 1+x^2;
	x2 = complexMultiplication(x, x);
	z = complexAddition({ 1, 0 }, x2);

	if (magnitude(z) > 10 * tolr) {
		//y = sqrt(z);
		z2polar(z, &r, &theta);
		y = polar2z(sqrt(r), theta / 2);
	}
	else
		y = { 1, 0 };

	//V2 = [1/y; x/y];
	//V2 = V2/norm(V2);
	V2[0] = complexDivision({ 1,0 }, y);
	V2[1] = complexDivision(x, y);

	nV2 = complexVectorNorm(V2, 2);

	if (nV2 > tolr) {
		V2[0] = complexDivision(V2[0], { nV2, 0.0 });
		V2[1] = complexDivision(V2[1], { nV2, 0.0 });
	}

	if (fabs(V2[0].imag) > tolr) {
		m1 = magnitude(V2[0]);
		m2 = magnitude(V2[1]);
		a1 = angle(V2[0]);
		a2 = angle(V2[1]);

		V2[0] = { m1, 0 };
		V2[1] = { m2*cos(a2 - a1),  m2*sin(a2 - a1) };
	}

	V[0] = V1[0];
	V[1] = V1[1];
	V[2] = V2[0];
	V[3] = V2[1];
}


// V, D = eig(A), where A is an N x N coomplex matrix
void eig(complex A[], unsigned int N, complex V[], complex D[])
{
	unsigned int i, j, k, complexEigenPairCount = 0;
	bool offDiagonals = false;
	double tolr = 2.3E-16, nv, mag1, mag2, ang1, ang2;
	complex R2[2 * 2], V2[2 * 2], D2[2];

	complex *H, *U, *Q, *R, *R_, *Q_, *W, *RR, *DD, *rRR, *qRR, *y, *v, *Wnew;
	unsigned int *complexPairs;
	//complex H[36], U[36], Q[36], R[36], R_[36], Q_[36], W[36], RR[36], DD[36], rRR[36], qRR[36], y[6], v[6], Wnew[36];
	//complex H[64], U[64], Q[64], R[64], R_[64], Q_[64], W[64], RR[64], DD[64], rRR[64], qRR[64], y[8], v[8], Wnew[64];
	//unsigned int complexPairs[6];
	//unsigned int complexPairs[8];

	H = (complex*)malloc(N*N * sizeof(complex));
	U = (complex*)malloc(N*N * sizeof(complex));
	Q = (complex*)malloc(N*N * sizeof(complex));
	R = (complex*)malloc(N*N * sizeof(complex));
	R_ = (complex*)malloc(N*N * sizeof(complex));
	Q_ = (complex*)malloc(N*N * sizeof(complex));
	W = (complex*)malloc(N*N * sizeof(complex));
	RR = (complex*)malloc(N*N * sizeof(complex));
	DD = (complex*)malloc(N*N * sizeof(complex));
	rRR = (complex*)malloc(N*N * sizeof(complex));
	qRR = (complex*)malloc(N*N * sizeof(complex));
	y = (complex*)calloc(N, sizeof(complex));
	v = (complex*)malloc(N * sizeof(complex));
	Wnew = (complex*)malloc(N*N * sizeof(complex));;
	complexPairs = (unsigned int*)calloc(N, sizeof(unsigned int));

	//Apply Hessenberg reduction
	hessenberg(A, N, H, U);

	//Apply Francis' implicit double shift QR algorithm
	doubleShift(H, N, Q, R);

	//get eigen values
	for (i = 0; i < N; i++)
		D[i] = R[lin_index(i, i, N)];

	// solve for 2x2 eigen values where off-diagonals exist
	for (i = 0; i < N - 1; i++) {
		if (magnitude(R[lin_index(i + 1, i, N)]) > tolr) {
			offDiagonals = true;

			R2[0] = R[lin_index(i, i, N)];
			R2[1] = R[lin_index(i, i + 1, N)];
			R2[2] = R[lin_index(i + 1, i, N)];
			R2[3] = R[lin_index(i + 1, i + 1, N)];
			solve2by2eig(R2, V2, D2);

			D[i] = D2[0];
			D[i + 1] = D2[1];
		}
	}

	//find which eigen values are complex
	i = 0;
	while (i < N) {
		if (fabs(D[i].imag) > tolr) {
			complexEigenPairCount++;
			complexPairs[i] = complexEigenPairCount;
			i++;
			complexPairs[i] = complexEigenPairCount;
		}
		i++;
	}

	//if no complex eigen values exist, use QR iteration to reduce to tridiagonal
	if (offDiagonals && complexEigenPairCount == 0) {
		qrIteration(R, N, Q_, R_);
		copyComplexMatrix(R_, N, N, R);
	}

	// solve upper tridiagonal system with backwards subsitution
	for (i = 0; i < N; i++) {
		complexIdentity(DD, N);
		complexScalarMultiplication(DD, D[i], N, N, DD);
		complexMatrixSubtraction(R, DD, N, N, RR);

		qrGivensComplex(RR, N, N, qRR, rRR);

		backwardsSubstitution(rRR, y, N, v);
		nv = complexVectorNorm(v,N);
		if (nv < tolr)
			nv = 1;

		complexScalarMultiplication(v, { 1/nv,0.0 }, N, 1, v);
		for (j = 0; j < N; j++)
			W[lin_index(j, i, N)] = v[j];
	}

	//calculate eigen vectors
	if (offDiagonals && complexEigenPairCount==0) {
		complexMatrixMultiplication(Q_, W, N, N, N, Wnew);
		copyComplexMatrix(Wnew, N, N, W);
	}
	complexMatrixMultiplication(Q, W, N, N, N, Wnew);
	complexMatrixMultiplication(U, Wnew, N, N, N, V);

	// rotate complex eigen vectors
	for (i = 0; i < N; i++) {
		if (complexPairs[i]) {
			if (i > 0 && complexPairs[i - 1] == complexPairs[i])
				k = i - 1;
			else
				k = i;
			mag1 = magnitude(V[lin_index(k, i, N)]);
			ang1 = angle(V[lin_index(k, i, N)]);
			for (j = 0; j < N; j++) {
				if (j == k)
					V[lin_index(k, i, N)] = { mag1, 0.0 };
				else {
					mag2 = magnitude(V[lin_index(j, i,N)]);
					ang2 = angle(V[lin_index(j, i, N)]);
					V[lin_index(j, i, N)] = { mag2*cos(ang2 - ang1), mag2*sin(ang2 - ang1) };
				}
			}
		}
	}

	//clear memory
	delete[] H;
	delete[] U;
	delete[] Q;
	delete[] R;
	delete[] Q_;
	delete[] R_;
	delete[] W;
	delete[] DD;
	delete[] RR;
	delete[] qRR;
	delete[] rRR;
	delete[] y;
	delete[] v;
	delete[] Wnew;
	delete[] complexPairs;
}
