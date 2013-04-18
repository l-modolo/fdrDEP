#include<R.h>
#include<Rdefines.h>
#include<iostream>

// Nota Bene: R stores matrices columnwise:
// transpar: 2x(1+2+p)
#define transpar(i,j) transpar[i+2*(j)]
// covariates: is a matrix of size m x p
#define covariates(i,j) covariates[i+m*(j)]
// A: 2 x 2 x m-1
#define A(i,j,k) A[i+2*(j)+4*(k)]
// dgammA: 2 x 2 x m-1
#define dgammA(i,j,k) dgammA[i+2*(j)+4*(k)]
// gamma: m x 2
#define gammA(i,j) gammA[i+m*(j)]

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
        gradient[3+j] =  (gammA(0,1) - pii[1])*covariates(0,j);
        for (unsigned int k=0; k<m-1; k++){
            gradient[3+j] += tmp[k]*covariates(k+1,j);
        }
    }
    delete[] tmp;
    if (*distancesincluded){
        gradient[3] =  (gammA(0,1) - pii[1])*covariates(0,0);
        for (unsigned int k=0; k<m-1; k++){
            gradient[3] += dgammA(0,1,k) * covariates(k+1,0)
                - dgammA(1,1,k) * covariates(k+1,0)
                - gammA(k,0)*A(0,1,k)*covariates(k+1,0)
                + gammA(k,1)*A(1,1,k)*covariates(k+1,0);
        }
    }
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
