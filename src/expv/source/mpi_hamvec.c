//==============================================================================
//  u = H * v + u;
//  H consist of the single and the pair terms of a many qubit system.
//==============================================================================

#include <complex.h>
#include <mpi.h>
#include <omp.h>

void mpi_hamvec(const unsigned int nspin, const double *coeff_lst, \
  const double _Complex *v, double _Complex *u, MPI_Request *reqs, \
  MPI_Status *status, double _Complex *cache, const unsigned int nchunk){
  
  int                       nprocs, rank, dest;
  unsigned int              np, rk;
  unsigned long             lproc;
  unsigned int              pn, qn, in, pm, qm, im;
  unsigned long             idxnx, idxny, idxnz, idxmx, idxmy, idxmz, idx;
  double                    cz;
  double _Complex           cxy;
  int                       chunk;
  //unsigned long             nchunk;
  
  unsigned int              n, m, k;
  unsigned long             i;
  
  double _Complex           *w, *wcache, *wex;
  
  unsigned int              nd;
  
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  np = nprocs;
  rk = rank;
  
  i = 0;
  while (np >>= 1){
    i++;
  }
  np = i;
  
  nd = nspin - np;
  
  lproc = 1 << nd;
  
  //chunk = 1048576;
  //nchunk = (lproc - 1) / chunk + 1;
  chunk = lproc / nchunk;
  
  w = cache;
  wcache = cache + lproc;
  
  //============================================================================
  //  idx = ip * 2^{q + 1} + im * 2^{q} + iq;
  //============================================================================
  
  //============================================================================
  //  determine the coefficient for each vector element;
  //  calculate:
  //    1.  sigma_{nz},                 0 <= n < nspin;
  //    2.  sigma_{mz} * sigma_{nz},    0 <= n < m < nspin;
  //  no need for mpi;
  //============================================================================
  #pragma omp parallel for private(cz, n, pn, in, idxnz, m, im, idxmz)
  for (i = 0; i < lproc; i++){
    cz = 0.0;
    
    for (n = 0; n < nspin; n++){
      
      //  for sigma_{nz};
      pn = nspin - 1 - n;
      in = (n < nd) ? ((i >> n) & 1) : ((rk >> (n - nd)) & 1);
      
      //  sigma_{nz}, withou mpi;
      idxnz = 3 * n + 9 * (((nspin + pn) * n) >> 1) + ((3 * pn + 1) << 1);
      cz += (in > 0) ? (-coeff_lst[idxnz]) : (+coeff_lst[idxnz]);
      
      for (m = n + 1; m < nspin; m++){
        
        //  for sigma_{mz};
        im = (m < nd) ? ((i >> m) & 1) : ((rk >> (m - nd)) & 1);
        
        //  sigma_{mz} * sigma_{nz}, withou mpi;
        idxmz = idxnz + 3 * (m - n);
        cz += ((in ^ im) > 0) ? (-coeff_lst[idxmz]) : (+coeff_lst[idxmz]);
      }// loop for m;
    }// loop for n;
    
    //  sigma_{nz} and sigma_{mz} * sigma_{nz};
    u[i] += cz * v[i];
  }// loop for i;
  
  //============================================================================
  //  determine the coefficient for each vector element;
  //  calculate:
  //    1.  sigma_{nxy},                0 <= n < nd;
  //    2.  sigma_{mz} * sigma_{nxy},   0 <= n < m < nd;
  //    3.  sigma_{nxy} * sigma_{mz},   0 <= m < n < nd;
  //    4.  sigma_{mxy} * sigma_{nxy},  0 <= n < m < nd;
  //  all operators in M_{D}, no need for mpi;
  //============================================================================
  #pragma omp parallel for private(cxy, n, in, idxnx, idxny, m, pm, im, \
    idxmx, idxmy, idx)
  for (i = 0; i < lproc; i++){
    
    //  for sigma_{nxy} in M_{D};
    for (n = 0; n < nd; n++){
      cxy = 0.0;
      
      //  for simga_{nxy};
      in = (i >> n) & 1;
      
      //  sigma_{nx}, withou mpi;
      idxnx = n * (3 + 9 * nspin) - 9 * (((n + 1) * n) >> 1);
      cxy += coeff_lst[idxnx];
      
      //  sigma_{ny}, without mpi;
      idxny = idxnx + 3 * (nspin - n) - 2;
      cxy += (in > 0) ? (+I * coeff_lst[idxny]) : (-I * coeff_lst[idxny]);
      
      //  for sigma_{mz} * sigma_{nxy}, m > n and m < n;
      for (m = 0; m < nspin; m++){
        if (m == n) continue;
        
        pm = nspin - 1 - m;
        im = (m < nd) ? ((i >> m) & 1) : ((rk >> (m - nd)) & 1);
        
        if (m > n){
          //  for sigma_{mz} * sigma_{nxy}, m > n;
          idxmx = idxnx + 3 * (m - n);
          idxmy = idxny + 3 * (m - n);
        }
        else {
          //  for sigma_{nxy} * sigma_{mz}, m < n;
          idxmx = 3 * n + ((9 * (nspin + pm) * m) >> 1) \
                + ((3 * pm + 1) << 1) - 2;
          idxmy = idxmx + 1;
        }
        
        //  sigma_{mz} * sigma_{nx}, without mpi;
        cxy += (im > 0) ? (-coeff_lst[idxmx]) : (+coeff_lst[idxmx]);
        //  sigma_{mz} * sigma_{ny}, without mpi;
        cxy += ((in ^ im) > 0) ? (+I * coeff_lst[idxmy]) \
                               : (-I * coeff_lst[idxmy]);
      }// loop for m;
      
      idx = (in > 0) ? (i - (1 << n)) : (i + (1 << n));
      
      //  sigma_{nxy} and sigma_{mz} * sigma_{nxy};
      u[i] += cxy * v[idx];
      
      //  for sigma_{mxy} * sigma_{nxy};
      for (m = n + 1; m < nd; m++){
        cxy = 0.0;
        
        im = (idx >> m) & 1;
        
        //  sigma_{mx} * sigma_{nx}, without mpi;
        idxmx = idxnx + 3 * (m - n) - 2;
        cxy += coeff_lst[idxmx];
        //  sigma_{my} * sigma_{nx}, without mpi;
        cxy += (im > 0) ? (+I * coeff_lst[idxmx + 1]) \
                        : (-I * coeff_lst[idxmx + 1]);
        
        //  sigma_{mx} * sigma_{ny}, without mpi;
        idxmy = idxny + 3 * (m - n) - 2;
        cxy += (in > 0) ? (+I * coeff_lst[idxmy]) : (-I * coeff_lst[idxmy]);
        //  sigma_{my} * sigma_{ny}, without mpi;
        cxy += ((in ^ im) > 0) ? (+coeff_lst[idxmy + 1]) \
                               : (-coeff_lst[idxmy + 1]);
        
        //  sigma_{mxy} * sigma_{nxy};
        u[i] += (im > 0) ? (cxy * v[idx - (1 << m)]) \
                         : (cxy * v[idx + (1 << m)]);
      }// loop for m;
      
    }// loop for n;
  }// loop for i;
  
  //============================================================================
  //  determine the coefficient for each composition;
  //  calculate:
  //    1.  sigma_{nxy},                nd <= n < nspin;
  //    2.  sigma_{mz} * sigma_{nxy},   nd <= n < nspin, m > n;
  //    3.  sigma_{nxy} * sigma_{mz},   nd <= n < nspin, m < n;
  //    4.  sigma_{nxy} * sigma_{mxy},  nd <= n < nspin, 
  //                                    and 0 <= m < nd;
  //  operators sigma_{nxy} in M_{N}, need for mpi; sigma_{mxy} in M_{D};
  //============================================================================
  //  prefetch vector for the first iteration;
  if (np > 0){
    n = nd;
    in = rk & 1;
    dest = (in > 0) ? (rank - 1) : (rank + 1);
    
    idxnx = n * (3 + 9 * nspin) - 9 * (((n + 1) * n) >> 1);
    //  start mpi sendrecv;
    for (k = 0; k < nchunk; k++){
      MPI_Irecv((void*)(wcache + chunk * k), 2 * chunk, MPI_DOUBLE, dest, \
        idxnx * nchunk + k, MPI_COMM_WORLD, reqs + 2 * k);
      MPI_Isend((void*)(v + chunk * k), 2 * chunk, MPI_DOUBLE, dest, \
        idxnx * nchunk + k, MPI_COMM_WORLD, reqs + 2 * k + 1);
    }
  }
  for (n = nd; n < nspin; n++){
    //  complete mpi sendrecv;
    MPI_Waitall(2 * nchunk, reqs, status);
    wex = w;
    w = wcache;
    wcache = wex;
    
    //  prefetch vector for the next iteration;
    if (n < nspin - 1){
      pn = nspin - 1 - (n + 1);
      qn = n + 1 - nd;
      in = (rk >> qn) & 1;
      dest = (in > 0) ? (rank - (1 << qn)) : (rank + (1 << qn));
      idxnx = 3 * (n + 1) + 9 * (((nspin + pn) * (n + 1)) >> 1);
      //  start mpi sendrecv;
      for (k = 0; k < nchunk; k++){
        MPI_Irecv((void*)(wcache + chunk * k), 2 * chunk, MPI_DOUBLE, dest, \
          idxnx * nchunk + k, MPI_COMM_WORLD, reqs + 2 * k);
        MPI_Isend((void*)(v + chunk * k), 2 * chunk, MPI_DOUBLE, dest, \
          idxnx * nchunk + k, MPI_COMM_WORLD, reqs + 2 * k + 1);
      }
    }
    
    //  for sigma_{nxy};
    pn = nspin - 1 - n;
    qn = n - nd;
    in = (rk >> qn) & 1;
    dest = (in > 0) ? (rank - (1 << qn)) : (rank + (1 << qn));
    
    idxnx = 3 * n + 9 * (((nspin + pn) * n) >> 1);
    idxny = idxnx + 3 * pn + 1;
    
    // determine the coefficient for each vector element;
    #pragma omp parallel for private(cxy, m, pm, qm, im, \
      idxmx, idxmy)
    for (i = 0; i < lproc; i++){
      cxy = 0.0;
      
      //  sigma_{nx}, nx in mpi;
      cxy += coeff_lst[idxnx];
      //  sigma_{ny}, ny in mpi;
      cxy += (in > 0) ? (I * coeff_lst[idxny]) : (-I * coeff_lst[idxny]);
    
      //  for sigma_{mz} * sigma_{nxy}, m > n and m < n;
      for (m = 0; m < nspin; m++){
        if (m == n) continue;
        
        pm = nspin - 1 - m;
        im = (m < nd) ? ((i >> m) & 1) : ((rk >> (m - nd)) & 1);
        
        if (m > n){
          //  for sigma_{mz} * sigma_{nxy}, m > n;
          idxmx = idxnx + 3 * (m - n);
          idxmy = idxny + 3 * (m - n);
        }
        else {
          //  for sigma_{nxy} * sigma_{mz}, m < n;
          idxmx = 3 * n + ((9 * (nspin + pm) * m) >> 1) \
                + ((3 * pm + 1) << 1) - 2;
          idxmy = idxmx + 1;
        }
        
        //  sigma_{mz} * sigma_{nx}, nx in mpi;
        cxy += (im > 0) ? (-coeff_lst[idxmx]) : (+coeff_lst[idxmx]);
        //  sigma_{mz} * sigma_{ny}, ny in mpi;
        cxy += ((in ^ im) > 0) ? (+I * coeff_lst[idxmy]) \
                               : (-I * coeff_lst[idxmy]);
      }// loop for m;
      
      //  sigma_{mz} * sigma_{nxy};
      u[i] += cxy * w[i];
      
      //  for sigma_{nxy} * sigma_{mxy}, nxy in mpi;
      for (m = 0; m < nd; m++){
        cxy = 0.0;
        
        im = (i >> m) & 1;
        
        idxmx = 9 * m * nspin - 9 * (((m + 1) * m) >> 1) + 3 * n - 2;
        idxmy = idxmx + 3 * (nspin - 1 - m) + 1;
        
        //  sigma_{mx} * sigma_{nx}, nx in mpi;
        cxy += coeff_lst[idxmx];
        //  sigma_{mx} * sigma_{ny}, ny in mpi;
        cxy += (in > 0) ? (+I * coeff_lst[idxmx + 1]) \
                        : (-I * coeff_lst[idxmx + 1]);
        
        //  sigma_{my} * sigma_{nx}, nx in mpi;
        cxy += (im > 0) ? (+I * coeff_lst[idxmy]) \
                        : (-I * coeff_lst[idxmy]);
        //  sigma_{my} * sigma_{ny}, ny in mpi;
        cxy += ((im ^ in) > 0) ? (+coeff_lst[idxmy + 1]) \
                               : (-coeff_lst[idxmy + 1]);
        
        //  sigma_{nxy} * sigma_{mxy};
        u[i] += (im > 0) ? (cxy * w[i - (1 << m)]) : (cxy * w[i + (1 << m)]);
      }// loop for m;
    }// loop for i;
  }// loop for n;
  
  //============================================================================
  //  determine the coefficient for each vector element;
  //  calculate:
  //    1.  sigma_{mxy} * sigma_{nxy},  nd <= n < m < nspin;
  //  sigma_{mxy} and sigma_{nxy} in M_{N}, need for mpi;
  //============================================================================
  //  sigma_{nxy}, sigma_{mxy} in M_{N};
  //  prefetch vector for the first iteration;
  if (np > 1){
    n = nd;
    pn = nspin - 1 - n;
    qn = n - nd;
    in = (rk >> qn) & 1;
    idx = (in > 0) ? (rank - (1 << qn)) : (rank + (1 << qn));
    idxnx = 3 * n + 9 * (((nspin + pn) * n) >> 1);
    idxny = idxnx + 3 * pn + 1;
    
    m = nd + 1;
    pm = nspin - 1 - m;
    qm = m - nd;
    im = (rk >> qm) & 1;
    dest = (im > 0) ? (idx - (1 << qm)) : (idx + (1 << qm));
    idxmx = idxny + 3 * (m - n) - 2;
    //  start mpi sendrecv;
    for (k = 0; k < nchunk; k++){
      MPI_Irecv((void *)(wcache + chunk * k), 2 * chunk, MPI_DOUBLE, dest, \
        idxmx * nchunk + k, MPI_COMM_WORLD, reqs + 2 * k);
      MPI_Isend((void *)(v + chunk * k), 2 * chunk, MPI_DOUBLE, dest, \
        idxmx * nchunk + k, MPI_COMM_WORLD, reqs + 2 * k + 1);
    }
  }
  for (n = nd; n < nspin - 1; n++){
    
    for (m = n + 1; m < nspin; m++){
      
      //  complete mpi sendrecv;
      MPI_Waitall(2 * nchunk, reqs, status);
      wex = w;
      w = wcache;
      wcache = wex;
      
      //  prefetch vector for the next iteration;
      if (m < nspin - 1){
        pn = nspin - 1 - n;
        qn = n - nd;
        in = (rk >> qn) & 1;
        idx = (in > 0) ? (rank - (1 << qn)) : (rank + (1 << qn));
        idxnx = 3 * n + 9 * (((nspin + pn) * n) >> 1);
        idxny = idxnx + 3 * pn + 1;
        
        pm = nspin - 1 - (m + 1);
        qm = (m + 1) - nd;
        im = (rk >> qm) & 1;
        dest = (im > 0) ? (idx - (1 << qm)) : (idx + (1 << qm));
        idxmx = idxny + 3 * ((m + 1) - n) - 2;
        //  start mpi sendrecv;
        for (k = 0; k < nchunk; k++){
          MPI_Irecv((void *)(wcache + chunk * k), 2 * chunk, MPI_DOUBLE, dest, \
            idxmx * nchunk + k, MPI_COMM_WORLD, reqs + 2 * k);
          MPI_Isend((void *)(v + chunk * k), 2 * chunk, MPI_DOUBLE, dest, \
            idxmx * nchunk + k, MPI_COMM_WORLD, reqs + 2 * k + 1);
        }
      }
      else if (n < nspin - 2){
        pn = nspin - 1 - (n + 1);
        qn = (n + 1) - nd;
        in = (rk >> qn) & 1;
        idx = (in > 0) ? (rank - (1 << qn)) : (rank + (1 << qn));
        idxnx = 3 * (n + 1) + 9 * (((nspin + pn) * (n + 1)) >> 1);
        idxny = idxnx + 3 * pn + 1;
        
        pm = nspin - 1 - (n + 2);
        qm = (n + 2) - nd;
        im = (rk >> qm) & 1;
        dest = (im > 0) ? (idx - (1 << qm)) : (idx + (1 << qm));
        idxmx = idxny + 3 * ((n + 2) - n) - 2;
        //  start mpi sendrecv;
        for (k = 0; k < nchunk; k++){
          MPI_Irecv((void *)(wcache + chunk * k), 2 * chunk, MPI_DOUBLE, dest, \
            idxmx * nchunk + k, MPI_COMM_WORLD, reqs + 2 * k);
          MPI_Isend((void *)(v + chunk * k), 2 * chunk, MPI_DOUBLE, dest, \
            idxmx * nchunk + k, MPI_COMM_WORLD, reqs + 2 * k + 1);
        }
      }
      
      pn = nspin - 1 - n;
      qn = n - nd;
      in = (rk >> qn) & 1;
      idx = (in > 0) ? (rank - (1 << qn)) : (rank + (1 << qn));
      idxnx = 3 * n + 9 * (((nspin + pn) * n) >> 1);
      idxny = idxnx + 3 * pn + 1;
      
      cxy = 0.0;
      
      pm = nspin - 1 - m;
      qm = m - nd;
      im = (rk >> qm) & 1;
      
      dest = (im > 0) ? (idx - (1 << qm)) : (idx + (1 << qm));
      
      //  sigma_{mx} * sigma_{nx};
      idxmx = idxnx + 3 * (m - n) - 2;
      cxy += coeff_lst[idxmx];
      //  sigma_{my} * sigma_{nx};
      cxy += (im > 0) ? (I * coeff_lst[idxmx + 1]) : (-I * coeff_lst[idxmx + 1]);
      
      //  sigma_{mx} * sigma_{ny};
      idxmy = idxny + 3 * (m - n) - 2;
      cxy += (in > 0) ? (I * coeff_lst[idxmy]) : (-I * coeff_lst[idxmy]);
      //  sigma_{my} * sigma_{ny};
      cxy += ((in ^ im) > 0) ? (coeff_lst[idxmy + 1]) : (-coeff_lst[idxmy + 1]);
      
      #pragma omp parallel for
      //  set each vector element;
      for (i = 0; i < lproc; i++)
        u[i] += cxy * w[i];
      
    }// loop of m;
  }// loop of n;
  
  return;
}

#define MPI_HAMVEC      mpi_hamvec_

#if defined(__cplusplus)
extern "C" {
#endif /* __cplusplus */
void MPI_HAMVEC(long* nspin_ptr, double *coeff_lst, \
  double _Complex *v, double _Complex *u, MPI_Request *reqs, \
  MPI_Status *status, double _Complex *cache, long *nchunk_ptr);
#if defined(__cplusplus)
}
#endif /* __cplusplus */

void MPI_HAMVEC(long* nspin_ptr, double *coeff_lst, \
  double _Complex *v, double _Complex *u, MPI_Request *reqs, \
  MPI_Status *status, double _Complex *cache, long *nchunk_ptr){
  
  unsigned int nspin = *nspin_ptr;
  unsigned int nchunk = *nchunk_ptr;
  mpi_hamvec(nspin, coeff_lst, v, u, reqs, status, cache, nchunk);
}
