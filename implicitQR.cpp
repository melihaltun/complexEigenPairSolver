
/**
* @file implicitQR.cpp
* @author Melih Altun @2015-2017
**/


#include "implicitQR.h"
#define EPS 2.221E-12

// for the math behind this implementation see:
// http://people.inf.ethz.ch/arbenz/ewp/Lnotes/2010/chapter3.pdf

// calculates Householder reflector vectors
void householder(complex x[], unsigned int n, complex u[])
{
	double mx, xc, sigma, s;
	complex *x2;
	unsigned int i;

	x2 = (complex*)malloc(n * sizeof(complex));

	mx = magnitude(x[0]);
	for (i = 1; i < n; i++) {
		xc = magnitude(x[i]);
		if (mx < xc)
			mx = xc;
	}

	for (i = 0; i < n; i++)
		x2[i] = complexDivision(x[i], { mx, 0.0 });

	s = copysign(1.0, x2[0].real);

	sigma = s * complexVectorNorm(x2, n);

	u[0] = complexAddition(x2[0], { sigma, 0.0 });
	for (i = 1; i < n; i++)
		u[i] = x2[i];

	delete[] x2;
}


// inner loop of the Francis' algorithm
void doubleShiftInnerLoop(complex H[], unsigned int N, unsigned int p, complex NH[], complex Q[])
{
	complex t, d, xyz[3], u3[3], c3[3*3], I3[3*3], u3_u3t[3*3], ut_u, u3t[3], I2[2*2], u2[2], u2t[2], u2_u2t[2*2], c2[2*2];
	complex *Pk, *Pkt,  *Qnew, *Hnew, *H_Pk;
	unsigned int k, j, i, MSD;
	double maxM, absM, tolr;

	Pk = (complex*)malloc(N*N * sizeof(complex));
	Pkt = (complex*)malloc(N*N * sizeof(complex));
	Qnew = (complex*)malloc(N*N * sizeof(complex));
	Hnew = (complex*)malloc(N*N * sizeof(complex));
	H_Pk = (complex*)malloc(N*N * sizeof(complex));

	complexIdentity(Q, N);
	memcpy(NH, H, N*N * sizeof(complex));

	//t = H(p - 1, p - 1) + H(p, p);  %trace
	t = complexAddition(NH[lin_index(p - 2, p - 2, N)], NH[lin_index(p - 1, p - 1, N)]);
	//d = H(p - 1, p - 1) * H(p, p) - H(p, p - 1) * H(p - 1, p);   % det
	d = complexSubtraction(complexMultiplication(NH[lin_index(p - 2, p - 2, N)], NH[lin_index(p - 1, p - 1, N)]), complexMultiplication(NH[lin_index(p - 1, p - 2, N)], NH[lin_index(p - 2, p - 1, N)]));
	//x = H(1, 1) ^ 2 + H(1, 2)*H(2, 1) - t* H(1, 1) + d;
	xyz[0] = complexMultiplication(NH[lin_index(0, 0, N)], NH[lin_index(0, 0, N)]);
	xyz[0] = complexAddition(xyz[0], complexMultiplication(NH[lin_index(0, 1, N)], NH[lin_index(1, 0, N)]));
	xyz[0] = complexSubtraction(xyz[0], complexMultiplication(t, NH[lin_index(0, 0, N)]));
	xyz[0] = complexAddition(xyz[0], d);
	//y = H(2, 1) * (H(1, 1) + H(2, 2) - t);
	xyz[1] = complexSubtraction(complexAddition(NH[lin_index(0, 0, N)], NH[lin_index(1, 1, N)]), t);
	xyz[1] = complexMultiplication(NH[lin_index(1,0,N)], xyz[1]);
	//z = H(2, 1) * H(3, 2);
	xyz[2] = complexMultiplication(NH[lin_index(1, 0, N)], NH[lin_index(2, 1, N)]);

	for (k = 0; k <= p - 3; k++) {
		householder(xyz, 3, u3);

		//c = eye(length(u)) - 2 * (u * u') / (u' * u);
		complexIdentity(I3, 3);
		hermitian(u3, 3, 1, u3t);
		complexMatrixMultiplication(u3, u3t, 3, 1, 3, u3_u3t);
		complexMatrixMultiplication(u3t, u3, 1, 3, 1, &ut_u);
		for (j = 0; j < 9; j++)
			u3_u3t[j] = complexDivision(u3_u3t[j], ut_u);
		complexScalarMultiplication(u3_u3t, { 2.0, 0.0 }, 3, 3, u3_u3t);
		complexMatrixSubtraction(I3, u3_u3t, 3, 3, c3);

		complexIdentity(Pk, N);

		//Pk(k + 1:k + 3, k + 1 : k + 3) = c;
		for (i = 0; i < 3; i++)
			for (j = 0; j < 3; j++)
				Pk[lin_index(k + i, k + j, N)] = c3[lin_index(i, j, 3)];

		//Q = Q*Pk
		complexMatrixMultiplication(Q, Pk, N, N, N, Qnew);
		copyComplexMatrix(Qnew, N, N, Q);

		//clear rounding errors in Q
		maxM = magnitude(Q[0]);
		for (i = 1; i < N*N; i++) {
			absM = magnitude(Q[i]);
			if (maxM < absM)
				maxM = absM;
		}

		MSD = ceil(log10(maxM));
		tolr = EPS * pow(10.0, (double)MSD);

		for (i = 1; i < N*N; i++) {
			if (tolr > fabs(Q[i].imag)) {
				Q[i].real = copysign(magnitude(Q[i]), Q[i].real);
				Q[i].imag = 0.0;

			}
		}

		//H = Pk'*H*Pk
		hermitian(Pk, N, N, Pkt);
		complexMatrixMultiplication(NH, Pk, N, N, N, H_Pk);
		complexMatrixMultiplication(Pkt, H_Pk, N, N, N, Hnew);
		copyComplexMatrix(Hnew, N, N, NH);

		//clear rounding errors in H
		maxM = magnitude(NH[0]);

		for (i = 1; i < N*N; i++) {
			absM = magnitude(NH[i]);
			if (maxM < absM)
				maxM = absM;
		}

		MSD = ceil(log10(maxM));
		tolr = EPS * pow(10.0, (double)MSD);

		for (i = 1; i < N*N; i++) {
			if (tolr > fabs(NH[i].imag)) {
				NH[i].real = copysign(magnitude(NH[i]), NH[i].real);
				NH[i].imag = 0.0;
			}
		}

		xyz[0] = NH[lin_index(k + 1, k, N)];
		xyz[1] = NH[lin_index(k + 2, k, N)];
		if (k < p-3)
			xyz[2] = NH[lin_index(k + 3, k, N)];
	}

	householder(xyz, 2, u2);  // get reflector vector

	//c = eye(length(u)) - 2 * (u * u') / (u' * u);
	complexIdentity(I2, 2);
	hermitian(u2, 2, 1, u2t);
	complexMatrixMultiplication(u2, u2t, 2, 1, 2, u2_u2t);
	complexMatrixMultiplication(u2t, u2, 1, 2, 1, &ut_u);
	for (j = 0; j < 4; j++)
		u2_u2t[j] = complexDivision(u2_u2t[j], ut_u);
	complexScalarMultiplication(u2_u2t, { 2.0, 0.0 }, 2, 2, u2_u2t);
	complexMatrixSubtraction(I2, u2_u2t, 2, 2, c2);

	complexIdentity(Pk, N);
	//Pk(p-1:p,p-1:p) = c;
	for (i = 0; i < 2; i++)
		for (j = 0; j < 2; j++)
			Pk[lin_index(k + i, k + j, N)] = c2[lin_index(i, j, 2)];

	//Q = Q*Pk
	complexMatrixMultiplication(Q, Pk, N, N, N, Qnew);
	copyComplexMatrix(Qnew, N, N, Q);

	//clear rounding errors
	maxM = magnitude(Q[0]);
	for (i = 1; i < N*N; i++) {
		absM = magnitude(Q[i]);
		if (maxM < absM)
			maxM = absM;
	}

	MSD = ceil(log10(maxM));
	tolr = EPS * pow(10.0, (double)MSD);

	for (i = 1; i < N*N; i++) {
		absM = magnitude(Q[i]);
		if (tolr > absM)
			Q[i] = { 0.0, 0.0 };
		else if (tolr > fabs(Q[i].imag)) {
			Q[i].real = copysign(magnitude(Q[i]), Q[i].real);
			Q[i].imag = 0.0;
		}
	}

	//H = Pk'*H*Pk
	hermitian(Pk, N, N, Pkt);
	complexMatrixMultiplication(NH, Pk, N, N, N, H_Pk);
	complexMatrixMultiplication(Pkt, H_Pk, N, N, N, Hnew);
	copyComplexMatrix(Hnew, N, N, NH);

	maxM = magnitude(NH[0]);

	//clear rounding errors
	for (i = 1; i < N*N; i++) {
		absM = magnitude(NH[i]);
		if (maxM < absM)
			maxM = absM;
	}

	MSD = ceil(log10(maxM));
	tolr = EPS * pow(10.0, (double)MSD);

	for (i = 1; i < N*N; i++) {
		absM = magnitude(NH[i]);
		if (tolr > absM)
			NH[i] = { 0.0, 0.0 };
		else if (tolr > fabs(NH[i].imag)) {
			NH[i].real = copysign(magnitude(NH[i]), NH[i].real);
			NH[i].imag = 0.0;
		}
	}

	// clear dynamic memory
	delete[] Pk;
	delete[] Pkt;
	delete[] Qnew;
	delete[] Hnew;
	delete[] H_Pk;
}


// outer loop of the implicit double shift QR algorithm
void doubleShift(complex H[], unsigned int N, complex Q[], complex R[])
{
	unsigned int p = N, maxIter = 1000, k = 0;
	double tolr = EPS;
	complex *HH, *NH, *QQ, *Qnew;
	
	HH = (complex*)malloc(N * N * sizeof(complex));
	NH = (complex*)malloc(N * N * sizeof(complex));
	QQ = (complex*)malloc(N * N * sizeof(complex));
	Qnew = (complex*)malloc(N * N * sizeof(complex));

	copyComplexMatrix(H, N, N, R);  //R = H

	complexIdentity(Q, N);  // Q = eye(N)

	while (p > 2) {
		k++;
		copyComplexMatrix(R, N, N, HH);  //HH = R

		doubleShiftInnerLoop(HH, N, p, NH, QQ);

		copyComplexMatrix(NH, N, N, R);  // R = NH

		complexMatrixMultiplication(Q, QQ, N, N, N, Qnew);
		copyComplexMatrix(Qnew, N, N, Q);   // Q = Q * QQ

		// reduce matrix if certain off-diagonals become zero
		if (magnitude(R[lin_index(p - 2, p - 3, N)]) < tolr * (magnitude(R[lin_index(p - 2, p - 2, N)]) + magnitude(R[lin_index(p - 3, p - 3, N)]))) {
			R[lin_index(p - 2, p - 3, N)] = equals({ 0,0 });
			p -= 2;
			continue;
		}

		if (magnitude(R[lin_index(p - 1, p - 2, N)]) < tolr * (magnitude(R[lin_index(p - 1, p - 1, N)]) + magnitude(R[lin_index(p - 2, p - 2, N)]))) {
			R[lin_index(p - 1, p - 2, N)] = equals({ 0,0 });
			p--;
		}

		if (k > maxIter)
			break;
	}

	// clear dynamic memory
	delete[] HH;
	delete[] QQ;
	delete[] NH;
	delete[] Qnew;
}
