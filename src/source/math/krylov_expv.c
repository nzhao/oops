#include <stdio.h>
#include <stdlib.h>
#include <complex>

#include "expokit.h"

#ifdef USE_MKL
  #define MKL_Complex16 double _Complex
  #include "mkl.h"
#else
  #include "cblas.h"
#endif

#ifdef USE_OMP
  #include "omp.h"
#endif

int krylov_zgexpv(const int n, const double _Complex *a, const double _Complex *v, const int tn, const double *ta, double _Complex *w_seq, const int klim, const int m, const double tol,  const int itrace)
{
  int             iflag = 0;
  
  double          anorm = 0.0;
  int             lwsp = 0, liwsp = 0, w_seq_len = 0;
  double _Complex *wsp  = NULL;
  int             *iwsp = NULL;
  
  lwsp      = n * (m + 2) + 5 * (m + 2) * (m + 2) + 7;
  liwsp     = m + 2;
  w_seq_len = n * tn;
  
  wsp = (double _Complex *) malloc(sizeof(double _Complex) * lwsp);
  if (wsp == NULL) return -1;
  iwsp = (int *) malloc(sizeof(int) * liwsp);
  if (iwsp == NULL) {
    free(wsp);
    return -1;
  }
  
  krylov_zgexpv_(&n, a, v, &tn, ta, &m, &tol, &anorm, wsp, &lwsp, iwsp, &liwsp, &itrace, &iflag, w_seq, &w_seq_len);
  
  free(iwsp);
  free(wsp);
  
  return iflag;
}

int krylov_zcooexpv(const int n, const int nz, const int *ia, const int *ja, const double _Complex *a, const double _Complex *v, const int tn, const double *ta, double _Complex *w_seq, const int klim, const int m, const double tol,  const int itrace)
{
  int             iflag = 0;
  
  double          anorm = 0.0;
  int             lwsp = 0, liwsp = 0, w_seq_len = 0;
  double _Complex *wsp  = NULL;
  int             *iwsp = NULL;
  
  lwsp      = n * (m + 2) + 5 * (m + 2) * (m + 2) + 7;
  liwsp     = m + 2;
  w_seq_len = n * tn;
  
  wsp = (double _Complex *) malloc(sizeof(double _Complex) * lwsp);
  if (wsp == NULL) return -1;
  iwsp = (int *) malloc(sizeof(int) * liwsp);
  if (iwsp == NULL){
    free(wsp);
    return -1;
  }
  
  krylov_zcooexpv_(&n, &nz, ia, ja, a, v, &tn, ta, &m, &tol, &anorm, wsp, &lwsp, iwsp, &liwsp, &itrace, &iflag, w_seq, &w_seq_len);
  
  free(iwsp);
  free(wsp);
  
  return iflag;
}

