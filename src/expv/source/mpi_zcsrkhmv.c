//==============================================================================
//  u = H * v + u;
//  H consist of the single and the pair terms of a many qubit system.
//==============================================================================

#include <complex.h>
#include <mpi.h>
#include <omp.h>

#include "ham.h"

void mpi_zcsrkhmv(const unsigned int nspin, const ham h, \
  const double _Complex coeff, const double _Complex *v, double _Complex *u, \
  const unsigned int nchunk, MPI_Request *reqs, \
  MPI_Status *status, double _Complex *cache){
  
  int                       nprocs, rank, dest, chunk;
  unsigned int              rk, np, nd;
  unsigned long             lproc, idx;
  unsigned int              qn, in, qm, im;
  double _Complex           c;
  
  unsigned int              n, m, k, j, nnz, jidx, nnzc, nnzcex;
  unsigned int              nidx, nidxex, nex, nnzex, qnex, inex, idxex;
  unsigned int              jex, jidxex, mex;
  unsigned long             l;
  
  double _Complex           *w, *wcache, *wex;
  
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  np = nprocs;
  rk = rank;
  
  l = 0;
  while (np >>= 1){
    l++;
  }
  np = l;
  
  nd = nspin - np;
  
  lproc = 1 << nd;
  
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
  //    2.  sigma_{nz} * sigma_{mz},    0 <= m < n < nspin;
  //============================================================================
  #pragma omp parallel for private(c, j, n, in, nnz, jidx, m, im)
  for (l = 0; l < lproc; l++){
    c = 0.0;
    
    //  for sigma_{nz}, 0 <= n < nspin;
    for (j = 0; j < h.nzia[1]; j++){
      n = h.nzja[j];
      in = (n < nd) ? ((l >> n) & 1) : ((rk >> (n - nd)) & 1);
      c += (in) ? (-h.nza[j]) : (+h.nza[j]);
    }
    
    //  for simga_{nz} * sigma_{mz}, 0 <= m < n < nspin;
    for (n = 0; n < nspin; n++){
      nnz = h.nzmzia[n + 1] - h.nzmzia[n];
      if (nnz){
        in = (n < nd) ? ((l >> n) & 1) : ((rk >> (n - nd)) & 1);
        for (j = 0; j < nnz; j++){
          jidx = h.nzmzia[n] + j;
          m = h.nzmzja[jidx];
          im = (m < nd) ? ((l >> m) & 1) : ((rk >> (m - nd)) & 1);
          c += (in ^ im) ? (-h.nzmza[jidx]) : (+h.nzmza[jidx]);
        }
      }
    }
    
    c *= coeff;
    u[l] += c * v[l];
  }
  
  //============================================================================
  //  determine the coefficient for each vector element;
  //  calculate:
  //    1.  sigma_{nxy},                0 <= n < nd;
  //    3.  sigma_{nxy} * sigma_{mz},   0 <= n < nd, 0 <= m < nspin;
  //    4.  sigma_{nxy} * sigma_{mxy},  0 <= m < n < nd;
  //============================================================================
  #pragma omp parallel for private(n, c, in, idx, j, nnz, jidx, m, im)
  for (l = 0; l < lproc; l++){
    for (n = 0; n < nd; n++){
      c = 0.0;
      in = (l >> n) & 1;
      idx = (in) ? (l - (1 << n)) : (l + (1 << n));
      
      //  for sigma_{nxy}, 0 <= n < nd;
      for (j = 0; j < h.nxyia[1]; j++){
        if (n == h.nxyja[j]){
          c += h.nxyax[j];
          c += (in) ? (+I * h.nxyay[j]) : (-I * h.nxyay[j]);
          break;
        }
      }
      
      //  for sigma_{nxy} * sigma_{mz}, 0 <= n < nd, 0 <= m < nspin;
      nnz = h.nxymzia[n + 1] - h.nxymzia[n];
      if (nnz){
        for (j = 0; j < nnz; j++){
          jidx = h.nxymzia[n] + j;
          m = h.nxymzja[jidx];
          im = (m < nd) ? ((l >> m) & 1) : ((rk >> (m - nd)) & 1);
          c += (im) ? (-h.nxymzax[jidx]) : (+h.nxymzax[jidx]);
          c += (in ^ im) ? (+I * h.nxymzay[jidx]) : (-I * h.nxymzay[jidx]);
        }
      }
      
      c *= coeff;
      u[l] += c * v[idx];
    }
  }
  #pragma omp parallel for private(n, nnz, in, idx, j, jidx, m, im, c)
  for (l = 0; l < lproc; l++){
    for (n = 0; n < nd; n++){
      
      //  for sigma_{nxy} * sigma_{mxy}, 0 <= m < n < nd;
      nnz = h.nxymxyia[n + 1] - h.nxymxyia[n];
      if (nnz){
        in = (l >> n) & 1;
        idx = (in) ? (l - (1 << n)) : (l + (1 << n));
        for (j = 0; j < nnz; j++){
          jidx = h.nxymxyia[n] + j;
          m = h.nxymxyja[jidx];
          im = (idx >> m) & 1;
          c = h.nxymxyaxx[jidx];
          c += (im) ? (+I * h.nxymxyaxy[jidx]) : (-I * h.nxymxyaxy[jidx]);
          c += (in) ? (+I * h.nxymxyayx[jidx]) : (-I * h.nxymxyayx[jidx]);
          c += (in ^ im) ? (+h.nxymxyayy[jidx]) : (-h.nxymxyayy[jidx]);
          
          c *= coeff;
          u[l] += (im) ? (c * v[idx - (1 << m)]) : (c * v[idx + (1 << m)]);
        }
      }
    }
  }
  
  //============================================================================
  //  determine the coefficient for each composition;
  //  calculate:
  //    1.  sigma_{nxy},                nd <= n < nspin;
  //    2.  sigma_{nxy} * sigma_{mz},   nd <= n < nspin, 0 <= m < nspin;
  //    4.  sigma_{nxy} * sigma_{mxy},  nd <= n < nspin, 0 <= m < nd;
  //============================================================================
  if (np > 0){
    
    //  prefetch vector for the first iteration;
    nidx = nspin;
    for (n = nd; n < nspin; n++){
      nnzc = 0;
      
      //  for sigma_{nxy}, nd <= n < nspin;
      for (j = 0; j < h.nxyia[1]; j++){
        if (n == h.nxyja[j]){
          nnzc++;
          break;
        }
        if (n < h.nxyja[j]) break;
      }
      
      //  for sigma_{nxy} * sigma_{mz}, nd <= n < nspin, 0 <= m < nspin;
      nnzc += h.nxymzia[n + 1] - h.nxymzia[n];
      
      //  for sigma_{nxy} * sigma_{mxy}, nd <= n < nspin, 0 <= m < nd;
      nnz = h.nxymxyia[n + 1] - h.nxymxyia[n];
      if (nnz && (h.nxymxyja[h.nxymxyia[n]] < nd)) nnzc++;
      
      if (nnzc){
        //  MPI sendrecv;
        qn = n - nd;
        in = (rk >> qn) & 1;
        dest = (in) ? (rank - (1 << qn)) : (rank + (1 << qn));
        for (k = 0; k < nchunk; k++){
          MPI_Irecv((void*)(wcache + chunk * k), 2 * chunk, MPI_DOUBLE, dest, \
            n * nchunk + k, MPI_COMM_WORLD, reqs + 2 * k);
          MPI_Isend((void*)(v + chunk * k), 2 * chunk, MPI_DOUBLE, dest, \
            n * nchunk + k, MPI_COMM_WORLD, reqs + 2 * k + 1);
        }
        nidx = n;
        break;
      }
    }
    
    for (n = nidx; n < nspin; n++){
      if (nnzc){
        //  MPI complete;
        MPI_Waitall(2 * nchunk, reqs, status);
        wex = w;
        w = wcache;
        wcache = wex;
      }
      
      //  prefetch vector for the next iteration;
      for (nex = n + 1; nex < nspin; nex++){
        nnzcex = 0;
        
        //  for sigma_{nxy}, nd <= nex < nspin;
        for (j = 0; j < h.nxyia[1]; j++){
          if (nex == h.nxyja[j]){
            nnzcex++;
            break;
          }
          if (nex < h.nxyja[j]) break;
        }
        
        //  for sigma_{nxy} * sigma_{mz}, nd <= nex < nspin, 0 <= mex < nspin;
        nnzcex += h.nxymzia[nex + 1] - h.nxymzia[nex];
        
        //  for sigma_{nxy} * sigma_{mxy}, nd <= nex < nspin, 0 <= mex < nd;
        nnzex = h.nxymxyia[nex + 1] - h.nxymxyia[nex];
        if (nnzex && (h.nxymxyja[h.nxymxyia[nex]] < nd)) nnzcex++;
        
        if (nnzcex){
          qn = nex - nd;
          in = (rk >> qn) & 1;
          dest = (in) ? (rank - (1 << qn)) : (rank + (1 << qn));
          for (k = 0; k < nchunk; k++){
            MPI_Irecv((void*)(wcache + chunk * k), 2 * chunk, MPI_DOUBLE, dest,\
              nex * nchunk + k, MPI_COMM_WORLD, reqs + 2 * k);
            MPI_Isend((void*)(v + chunk * k), 2 * chunk, MPI_DOUBLE, dest, \
              nex * nchunk + k, MPI_COMM_WORLD, reqs + 2 * k + 1);
          }
          break;
        }
      }
      
      if (nnzc){
        qn = n - nd;
        in = (rk >> qn) & 1;
        
        #pragma omp parallel for private(c, j, jidx, m, im)
        for (l = 0; l < lproc; l++){
          c = 0.0;
          
          //  sigma_{nxy}, nd <= n < nspin;
          for (j = 0; j < h.nxyia[1]; j++){
            if (n == h.nxyja[j]){
              c += h.nxyax[j];
              c += (in) ? (+I * h.nxyay[j]) : (-I * h.nxyay[j]);
              break;
            }
            if (n < h.nxyja[j]) break;
          }
          
          //  sigma_{nxy} * sigma_{mz}, nd <= n < nspin, 0 <= m < nspin;
          for (j = 0; j < (h.nxymzia[n + 1] - h.nxymzia[n]); j++){
            jidx = h.nxymzia[n] + j;
            m = h.nxymzja[jidx];
            im = (m < nd) ? ((l >> m) & 1) : ((rk >> (m - nd)) & 1);
            c += (im) ? (-h.nxymzax[jidx]) : (+h.nxymzax[jidx]);
            c += (in ^ im) ? (+I * h.nxymzay[jidx]) : (-I * h.nxymzay[jidx]);
          }
          
          c *= coeff;
          u[l] += c * w[l];
        }
        
        #pragma omp parallel for private(j, jidx, m, im, c)
        for (l = 0; l < lproc; l++){
          //  sigma_{nxy} * sigma_{mxy}, nd <= n < nspin, 0 <= m < nd;
          for (j = 0; j < (h.nxymxyia[n + 1] - h.nxymxyia[n]); j++){
            jidx = h.nxymxyia[n] + j;
            m = h.nxymxyja[jidx];
            if (m >= nd) break;
            im = (l >> m) & 1;
            c = h.nxymxyaxx[jidx];
            c += (im) ? (+I * h.nxymxyaxy[jidx]) : (-I * h.nxymxyaxy[jidx]);
            c += (in) ? (+I * h.nxymxyayx[jidx]) : (-I * h.nxymxyayx[jidx]);
            c += (in ^ im) ? (+h.nxymxyayy[jidx]) : (-h.nxymxyayy[jidx]);
            
            c *= coeff;
            u[l] += (im) ? (c * w[l - (1 << m)]) : (c * w[l + (1 << m)]);
          }
        }
      }
      nnzc = nnzcex;
    }
  }
  
  //============================================================================
  //  determine the coefficient for each vector element;
  //  calculate:
  //    1.  sigma_{mxy} * sigma_{nxy},  nd <= m < n < nspin;
  //============================================================================
  if (np > 1){
    
    //  prefetch vector for the first iteration;
    nidx = 0;
    nidxex = 0;
    for (nex = nd; nex < nspin; nex++){
      //  sigma_{nxy} * sigma_{mxy}, nd <= m < n < nspin;
      nnzex = h.nxymxyia[nex + 1] - h.nxymxyia[nex];
      if (nnzex){
        qnex = nex - nd;
        inex = (rk >> qnex) & 1;
        idxex = (inex) ? (rank - (1 << qnex)) : (rank + (1 << qnex));
        for (jex = 0; jex < nnzex; jex++){
          jidxex = h.nxymxyia[nex] + jex;
          mex = h.nxymxyja[jidxex];
          if (mex < nd) continue;
          nidxex = nspin * nex + mex;
          qm = mex - nd;
          im = (rk >> qm) & 1;
          dest = (im) ? (idxex - (1 << qm)) : (idxex + (1 << qm));
          //  MPI sendrecv;
          for (k = 0; k < nchunk; k++){
            MPI_Irecv((void*)(wcache + chunk * k), 2 * chunk, MPI_DOUBLE, dest,\
              nidxex * nchunk + k, MPI_COMM_WORLD, reqs + 2 * k);
            MPI_Isend((void*)(v + chunk * k), 2 * chunk, MPI_DOUBLE, dest, \
              nidxex * nchunk + k, MPI_COMM_WORLD, reqs + 2 * k + 1);
          }
          break;
        }
        if (nidx < nidxex){
          nidx = nidxex;
          break;
        }
      }
    }
    
    for (n = nd; n < nspin; n++){
      //  sigma_{nxy} * sigma_{mxy}, nd <= m < n < nspin;
      nnz = h.nxymxyia[n + 1] - h.nxymxyia[n];
      if (nnz){
        qn = n - nd;
        in = (rk >> qn) & 1;
        for (j = 0; j < nnz; j++){
          jidx = h.nxymxyia[n] + j;
          m = h.nxymxyja[jidx];
          if ((nspin * n + m) < nidx) continue;
          
          //  MPI complete;
          MPI_Waitall(2 * nchunk, reqs, status);
          wex = w;
          w = wcache;
          wcache = wex;
          
          //  prefetch vector for the next iteration;
          for (nex = n; nex < nspin; nex++){
            nnzex = h.nxymxyia[nex + 1] - h.nxymxyia[nex];
            if (nnzex){
              qnex = nex - nd;
              inex = (rk >> qnex) & 1;
              idxex = (inex) ? (rank - (1 << qnex)) : (rank + (1 << qnex));
              for (jex = 0; jex < nnzex; jex++){
                jidxex = h.nxymxyia[nex] + jex;
                mex = h.nxymxyja[jidxex];
                if (mex < nd) continue;
                nidxex = nspin * nex + mex;
                if (nidxex <= nidx) continue;
                qm = mex - nd;
                im = (rk >> qm) & 1;
                dest = (im) ? (idxex - (1 << qm)) : (idxex + (1 << qm));
                //  MPI sendrecv;
                for (k = 0; k < nchunk; k++){
                  MPI_Irecv((void*)(wcache + chunk * k), 2 * chunk, MPI_DOUBLE,\
                    dest, nidxex * nchunk + k, MPI_COMM_WORLD, \
                    reqs + 2 * k);
                  MPI_Isend((void*)(v + chunk * k), 2 * chunk, MPI_DOUBLE, \
                    dest, nidxex * nchunk + k, MPI_COMM_WORLD, \
                    reqs + 2 * k + 1);
                }
                break;
              }
              if (nidx < nidxex){
                nidx = nidxex;
                break;
              }
            }
          }
          
          qm = m - nd;
          im = (rk >> qm) & 1;
          c = h.nxymxyaxx[jidx];
          c += (im) ? (+I * h.nxymxyaxy[jidx]) : (-I * h.nxymxyaxy[jidx]);
          c += (in) ? (+I * h.nxymxyayx[jidx]) : (-I * h.nxymxyayx[jidx]);
          c += (in ^ im) ? (+h.nxymxyayy[jidx]) : (-h.nxymxyayy[jidx]);
          
          c *= coeff;
          #pragma omp parallel for
          for (l = 0; l < lproc; l++)
            u[l] += c * w[l];
        }
      }
    }
  }
  
  return;
}
#define MPI_ZCSRKHMV    mpi_zcsrkhmv_

#if defined(__cplusplus)
extern "C" {
#endif /* __cplusplus */
void MPI_ZCSRKHMV(long* nspin_ptr, long *h_ptr, double _Complex *coeff_ptr,\
  double _Complex *v, double _Complex *u, long *nchunk_ptr, MPI_Request *reqs, \
  MPI_Status *status, double _Complex *cache);
#if defined(__cplusplus)
}
#endif /* __cplusplus */

void MPI_ZCSRKHMV(long* nspin_ptr, long *h_ptr, double _Complex *coeff_ptr,\
  double _Complex *v, double _Complex *u, long *nchunk_ptr, MPI_Request *reqs, \
  MPI_Status *status, double _Complex *cache){
  
  unsigned int nspin = *nspin_ptr;
  unsigned int nchunk = *nchunk_ptr;
  ham h = *((ham*)h_ptr);
  double _Complex coeff = *coeff_ptr;
  
/*  // for debug;*/
/*  printf("[zcsrkhmv] nspin  = %d\n", nspin);*/
/*  printf("[zcsrkhmv] nchunk = %d\n", nchunk);*/
/*  printf("[zcsrkhmv] h_ptr  = %p\n", h_ptr);*/
/*  printf("[zcsrkhmv] *h_ptr = %p\n", *h_ptr);*/
/*  printf("[zcsrkhmv] &h     = %p\n", &h);*/
/*  printf("[zcsrkhmv] h.nzia = %p\n", h.nzia);*/
  
  mpi_zcsrkhmv(nspin, h, coeff, v, u, nchunk, reqs, status, cache);
}
