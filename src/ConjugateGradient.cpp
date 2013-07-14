#include<R.h>
#include<Rdefines.h>
#include<iostream>
#include <stdlib.h>

// Nota Bene: R stores matrices columnwise

// dimensions: p covariates, m points
// transpar: 2x(1+2+p)
#define transpar(i,j) transpar[i+2*(j)]
#define transparnew(i,j) transparnew[i+2*(j)]
#define transparold(i,j) transparold[i+2*(j)]
// covariates: is a matrix of size m x p
#define covariates(i,j) covariates[i+m*(j)]
// A: 2 x 2 x m-1
#define A(i,j,k) A[i+2*(j)+4*(k)]
// dgammA: 2 x 2 x m-1
#define dgammA(i,j,k) dgammA[i+2*(j)+4*(k)]
// gamma: m x 2
#define gammA(i,j) gammA[i+m*(j)]
// phi: 2x(1+2+p)
#define phi(i,j) phi[i+2*(j)]

extern "C" {

void C_piiA(int* m_ptr, int* p_ptr,
	double* covariates,
	int* distancesincluded,
	double* transpar,
	int* v,
	double* pii,
	double* A)
	{
	unsigned int p = *p_ptr, m = *m_ptr;
	unsigned int q = 1+2+p;
	double p0, p1;
	p0 = transpar(0,0);
	p1 = transpar(1,0);
	for (unsigned int j=0; j<p; j++){
		p0 += transpar(0,3+j) * covariates(0,j);
		p1 += transpar(1,3+j) * covariates(0,j);
	}
	pii[0] = exp(p0) / (exp(p0)+exp(p1));
	pii[1] = 1 - pii[0];
	
	double num11 = exp(transpar(0,1));
	double num21 = exp(transpar(0,2));
	for (unsigned int k=0; k<m-1; k++){
		double denom12 = 0.;
		for (unsigned int j=0; j<p; j++){
			denom12 += transpar(1,3+j) * covariates(k+1,j);
		}
		double denom22 = denom12;
		if (*distancesincluded){
			denom22 -= transpar(1,3) * covariates(k+1,0); // undo plus
			denom22 -= transpar(1,3) * covariates(k+1,0); // and minus
		}
		double num12 = exp(transpar(1,1) + denom12);
		double num22 = exp(transpar(1,2) + denom22);
		A(0,0,k) = num11/(num11 + num12);
		A(0,1,k) = num12/(num11 + num12);
		A(1,0,k) = num21/(num21 + num22);
		A(1,1,k) = num22/(num21 + num22);
}
}

void C_LineSearch(int* m_ptr, int* p_ptr,
	double* covariates,
	int* distancesincluded,
	double* dgammA,
	double* gammA,
	double* transpar,
	double* phi,
	int* iterCG,
	double* ptol,
	int* v,
	double* nu,
	double* transparnew,
	double* pii,
	double* A)
{
	unsigned int p = *p_ptr, m = *m_ptr;
	int i, j, n, k;
	int niter = 0;
	double difference = 1.0;
	double dQ, dQ2, dQ_tmp, dQ2_tmp1, dQ2_tmp2, tmp;
	
	while (difference > *ptol && niter < *iterCG)
	{
		niter++;
		C_piiA( m_ptr, p_ptr, covariates, distancesincluded, transparnew, v, pii, A);
		
		// dQ
		dQ = 0.0;
		for (j=0; j < 2; j++)
		{
			dQ_tmp = 0.0;
			for (k=3; k < p+3; k++)
			{
				if (*distancesincluded && k == 3)
					dQ_tmp += -phi(j,k) * covariates(0,k-3);
				else
					dQ_tmp += phi(j,k) * covariates(0,k-3);
			}
			dQ += (phi(j,0) + dQ_tmp) * (gammA(0,j) - pii[j]);
		}
		
		for (i=0; i < 2; i++)
		{
			for (j=0; j < 2; j++)
			{
				for (n = 1; n < m; n++)
				{
					dQ_tmp = 0.0;
					for (k=3; k < p+3; k++)
					{
						if (*distancesincluded && i == 1 && k == 3)
							dQ_tmp += -phi(j,k) * covariates(n, k-3);
						else
							dQ_tmp += phi(j,k) * covariates(n, k-3);
					}
					dQ += (phi(j,i+1) + dQ_tmp) * (dgammA(i,j,n) - gammA(n-1,i) * A(i,j,n));
				}
			}
		}
//		Rprintf("dQ : %f\n", dQ);
		
		// dQ2
		dQ2 = 0.0;
		for (j=0; j < 2; j++)
		{
			dQ2_tmp1 = 0.0;
			for (k=3; k < p+3; k++)
			{
				if (*distancesincluded && k == 3)
					dQ2_tmp1 += -phi(j,k) * covariates(0,k-3);
				else
					dQ2_tmp1 += phi(j,k) * covariates(0,k-3);
			}
			dQ2 += (phi(j,0) + dQ2_tmp1) * (phi(j,0) + dQ2_tmp1) * pii[j] * (1 - pii[j]);
		}
		dQ2 = -dQ2;
		
		for (i=0; i < 2; i++)
		{
			for (j=0; j < 2; j++)
			{
				dQ2_tmp2 = 0.0;
				for (n = 1; n < m; n++)
				{
					dQ2_tmp1 = 0.0;
					for (k=3; k < p+3; k++)
					{
						if (*distancesincluded && i == 1 && k == 3)
							dQ2_tmp1 += -phi(j,k) * covariates(n, k-3);
						else
							dQ2_tmp1 += phi(j,k) * covariates(n, k-3);
					}
					dQ2_tmp2 += (phi(j,i+1) + dQ2_tmp1) * (phi(j,i+1) + dQ2_tmp1) * A(i,j,n) * (1 - A(i,j,n)) * gammA(n-1,i);
				}
				dQ2 -= dQ2_tmp2;
			}
		}
//		Rprintf("dQ2 : %f\n", dQ2);
//		dQ2 = - dQ2 - dQ2_tmp2;
		
		if (dQ2 == 0 && v)
		{
			Rprintf("Error in line search\n");
			error("Error in line search\n");
			difference = *ptol;
			niter = *iterCG;
			*nu = 0;
			for (j=0; j < 2; j++)
			{
				for (k=0; k < 3+p; k++)
				{
					transparnew(j,k) = transpar(j,k);
				}
			}
		}
		else
		{
			*nu = *nu - dQ / dQ2;
			difference = abs(dQ/ dQ2);
			for (j=0; j < 2; j++)
			{
				for (k=0; k < 3+p; k++)
				{
					transparnew(j,k) = transpar(j,k) + *nu * phi(j,k);
				}
			}
		}
	}
}

void C_ComputeGradient(int* m_ptr, int* p_ptr,
	double* covariates,
	int* distancesincluded,
	double* dgammA,
	double* gammA,
	double* transpar,
	int* v,
	double* pii,
	double* A,
	double* gradient)
{
	unsigned int p = *p_ptr, m = *m_ptr;
	C_piiA(m_ptr, p_ptr, covariates, distancesincluded, transpar, v, pii, A);
	gradient[0] = gammA(0,1) - pii[1];
	for (unsigned int k=0; k<m-1; k++){
		gradient[1] += dgammA(0,1,k) - gammA(k,0)*A(0,1,k);
		gradient[2] += dgammA(1,1,k) - gammA(k,1)*A(1,1,k);
	}
	double* tmp = new double[m-1];
	for (unsigned int k=0; k<m-1; k++){
		tmp[k] = gammA(k+1,1) - (gammA(k,0)*A(0,1,k) + gammA(k,1)*A(1,1,k));
	}
	unsigned int first = 0;
	if (*distancesincluded){
		first = 1; // gradient[3] not concerned
	}
	for (unsigned int j=first; j<p; j++){
		gradient[3+j] =(gammA(0,1) - pii[1])*covariates(0,j);
		for (unsigned int k=0; k<m-1; k++){
			gradient[3+j] += tmp[k]*covariates(k+1,j);
		}
	}
	delete[] tmp;
	if (*distancesincluded){
		gradient[3] =(gammA(0,1) - pii[1])*covariates(0,0);
		for (unsigned int k=0; k<m-1; k++){
			gradient[3] += dgammA(0,1,k) * covariates(k+1,0)
				- dgammA(1,1,k) * covariates(k+1,0)
				- gammA(k,0)*A(0,1,k)*covariates(k+1,0)
				+ gammA(k,1)*A(1,1,k)*covariates(k+1,0);
		}
	}
}

void C_ComputeCG(
	int* m_ptr, int* p_ptr, double* covariates, int* distancesincluded,
	double* dgammA, double* gammA,
	double* transpar,
	double* transparnew,
	double* transparold,
	double* phi,
	double* nu,
	int* iterCG,
	double* ptol,
	double* gradientnew,
	double* gradientold,
	int* v, double* pii, double* A
	)
{
	unsigned int p = *p_ptr, m = *m_ptr;
	int difference = 1;
	int niter = 0;
	int i, j, k;
	double PR;
	C_ComputeGradient(m_ptr, p_ptr, covariates, distancesincluded, dgammA, gammA, transpar, v, pii, A, gradientold);
	
	while (difference > *ptol && niter < *iterCG)
	{
		for (j=0; j < 2; j++)
		{
			for (k=0; k < 3+p; k++)
			{
				transparold(j,k) = transpar(j,k);
			}
		}
		niter++;
		C_LineSearch(m_ptr, p_ptr, covariates, distancesincluded, dgammA, gammA, transpar, phi, iterCG, ptol, v, nu, transparnew, pii, A);
		
		for (j=0; j < 2; j++)
		{
			for (k=0; k < 3+p; k++)
			{
				transpar(j,k) = transparnew(j,k);
			}
		}
		C_ComputeGradient(m_ptr, p_ptr, covariates, distancesincluded, dgammA, gammA, transpar, v, pii, A, gradientnew);
		
		PR = 0.0;
		for (k = 0; k < 3+p; k++)
		{
			PR += (( gradientnew[k] - gradientold[k] ) * gradientnew[k]) / ( gradientold[k] * gradientold[k]);
		}
		if (PR < 0)
		{
			PR = 0.0;
		}
		
		for (j=0; j < 2; j++)
		{
			for (k=0; k < 3+p; k++)
			{
				if ( j == 0 )
					phi(j,k) = 0;
				else
					phi(j,k) = -gradientnew[k] + PR * phi(j,k);
			}
		}
		
		for (j = 0; j < 2; j++)
		{
			if (*distancesincluded)
				transpar(j,3) = abs(transpar(j,3));
		}
		if (*distancesincluded)
			transpar(2,4) = abs(transpar(2,4));
		
		difference = 0.0;
		for (k = 1; k < 3+p; k++)
		{
			if (difference < abs(transpar(2,k) - transparold(2,k)))
				difference = abs(transpar(2,k) - transparold(2,k));
		}
	}
	C_piiA( m_ptr, p_ptr, covariates, distancesincluded, transpar, v, pii, A);
}

void Call_piiA(SEXP A)
{
	if (!isMatrix(A))
		error("a matrix is required");
	double *pA;
	PROTECT(A = AS_NUMERIC(A));
	pA = NUMERIC_POINTER(A);
	pA[0] = 77;
	UNPROTECT(1);
}

}
