#include <stdlib.h>

#include "ham.h"

//  get the index of c_{n beta} in coeff_lst[:];
unsigned int getSingleIdx(const unsigned int nspin, const unsigned int n, \
  const unsigned int beta){
  
  unsigned int  len;
  
  if (n > nspin)
    return (3 * nspin + 9 * nspin * (nspin - 1) / 2);
  
  len = nspin - 1 - n;
  
  return (3 * n + 9 * (nspin + len) * n / 2 + len * 3 * beta + beta);
}

//  get the index of c_{m alpha n beta} in coeff_lst[:];
unsigned int getPairIdx(const unsigned int nspin, const unsigned int n, \
  const unsigned int beta, const unsigned int m, const unsigned int alpha){
  
  if (n > nspin || m > nspin)
    return (3 * nspin + 9 * nspin * (nspin - 1) / 2);
  
  return (getSingleIdx(nspin, n, beta) + 1 + 3 * (m - 1 - n) + alpha);
}

ham hamGen(int nspin, double *coeff_lst){
  ham       h;
  int       i, j, nnz;
  double    a, ax, ay, axx, axy, ayx, ayy;
  
  h.nzia = (unsigned int*)malloc(sizeof(unsigned int) * 2);
  h.nzja = (unsigned int*)malloc(sizeof(unsigned int) * nspin);
  h.nza = (double*)malloc(sizeof(double) * nspin);
  
  h.nxyia = (unsigned int*)malloc(sizeof(unsigned int) * 2);
  h.nxyja = (unsigned int*)malloc(sizeof(unsigned int) * nspin);
  h.nxyax = (double*)malloc(sizeof(double) * nspin);
  h.nxyay = (double*)malloc(sizeof(double) * nspin);
  
  h.nxymzia = (unsigned int*)malloc(sizeof(unsigned int) * (nspin + 1));
  h.nxymzja = (unsigned int*)malloc(sizeof(unsigned int) * nspin * nspin);
  h.nxymzax = (double*)malloc(sizeof(double) * nspin * nspin);
  h.nxymzay = (double*)malloc(sizeof(double) * nspin * nspin);
  
  i = nspin * (nspin - 1) / 2;
  h.nzmzia = (unsigned int*)malloc(sizeof(unsigned int) * (nspin + 1));
  h.nzmzja = (unsigned int*)malloc(sizeof(unsigned int) * i);
  h.nzmza = (double*)malloc(sizeof(double) * i);
  
  h.nxymxyia = (unsigned int*)malloc(sizeof(unsigned int) * (nspin + 1));
  h.nxymxyja = (unsigned int*)malloc(sizeof(unsigned int) * i);
  h.nxymxyaxx = (double*)malloc(sizeof(double) * i);
  h.nxymxyaxy = (double*)malloc(sizeof(double) * i);
  h.nxymxyayx = (double*)malloc(sizeof(double) * i);
  h.nxymxyayy = (double*)malloc(sizeof(double) * i);
  
  //  convert coeff_lst to crs format;
  
  nnz = 0;
  h.nzia[0] = nnz;
  for (j = 0; j < nspin; j++){
    a = coeff_lst[getSingleIdx(nspin, j, 2)];
    if (a){
      h.nzja[nnz] = j;
      h.nza[nnz] = a;
      nnz++;
    }
  }
  h.nzia[1] = nnz;
  
  nnz = 0;
  h.nxyia[0] = nnz;
  for (j = 0; j < nspin; j++){
    ax = coeff_lst[getSingleIdx(nspin, j, 0)];
    ay = coeff_lst[getSingleIdx(nspin, j, 1)];
    if (ax || ay){
      h.nxyja[nnz] = j;
      h.nxyax[nnz] = ax;
      h.nxyay[nnz] = ay;
      nnz++;
    }
  }
  h.nxyia[1] = nnz;
  
  nnz = 0;
  for (i = 0; i < nspin; i++){
    h.nzmzia[i] = nnz;
    for (j = 0; j < i; j++){
      a = coeff_lst[getPairIdx(nspin, j, 2, i, 2)];
      if (a){
        h.nzmzja[nnz] = j;
        h.nzmza[nnz] = a;
        nnz++;
      }
    }
  }
  h.nzmzia[nspin] = nnz;
  
  nnz = 0;
  for (i = 0; i < nspin; i++){
    h.nxymzia[i] = nnz;
    for (j = 0; j < nspin; j++ ){
      if (i == j) continue;
      if (i < j){
        ax = coeff_lst[getPairIdx(nspin, i, 0, j, 2)];
        ay = coeff_lst[getPairIdx(nspin, i, 1, j, 2)];
      }
      else {
        ax = coeff_lst[getPairIdx(nspin, j, 2, i, 0)];
        ay = coeff_lst[getPairIdx(nspin, j, 2, i, 1)];
      }
      if (ax || ay){
        h.nxymzja[nnz] = j;
        h.nxymzax[nnz] = ax;
        h.nxymzay[nnz] = ay;
        nnz++;
      }
    }
  }
  h.nxymzia[nspin] = nnz;
  
  nnz = 0;
  for (i = 0; i < nspin; i++){
    h.nxymxyia[i] = nnz;
    for (j = 0; j < i; j++){
      axx = coeff_lst[getPairIdx(nspin, j, 0, i, 0)];
      axy = coeff_lst[getPairIdx(nspin, j, 1, i, 0)];
      ayx = coeff_lst[getPairIdx(nspin, j, 0, i, 1)];
      ayy = coeff_lst[getPairIdx(nspin, j, 1, i, 1)];
      if (axx || axy || ayx || ayy){
        h.nxymxyja[nnz] = j;
        h.nxymxyaxx[nnz] = axx;
        h.nxymxyaxy[nnz] = axy;
        h.nxymxyayx[nnz] = ayx;
        h.nxymxyayy[nnz] = ayy;
        nnz++;
      }
    }
  }
  h.nxymxyia[nspin] = nnz;
  
  return h;
}

void hamFree(ham h){
  
  free(h.nxymxyayy);
  free(h.nxymxyayx);
  free(h.nxymxyaxy);
  free(h.nxymxyaxx);
  free(h.nxymxyja);
  free(h.nxymxyia);
  free(h.nxymzay);
  free(h.nxymzax);
  free(h.nxymzja);
  free(h.nxymzia);
  free(h.nzmza);
  free(h.nzmzja);
  free(h.nzmzia);
  free(h.nxyay);
  free(h.nxyax);
  free(h.nxyja);
  free(h.nxyia);
  free(h.nza);
  free(h.nzja);
  free(h.nzia);
  
  return;
}
