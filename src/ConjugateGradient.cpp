#include<R.h>
#include<iostream>
extern "C" {

// Nota Bene: R stores matrices columnwise:
// transpar: 2x(1+2+p)
#define transpar(i,j) transpar[i+2*(j)]
// covariates: is a matrix of size m x p
#define covariates(i,j) covariates[i+m*(j)]
// A: 2 x 2 x m-1
#define A(i,j,k) A[i+2*(j)+4*(k)]

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
      double num12 = exp(transpar(1,1) + denom12);
      double num22 = exp(transpar(1,2) + denom22);
      A(0,0,k) = num11/(num11 + num12);
      A(0,1,k) = num12/(num11 + num12);
      A(1,0,k) = num21/(num21 + num22);
      A(1,1,k) = num22/(num21 + num22);
  }
}

}
