// |w> = |w> + H * |v>;
// H = sum(H_i) + sum(H_ij);
// H_i = h_i * S_(i,xyz) or h_ij * S_(i,xyz) * S_(j,xyz);
// w_med = H_ij * w_med 

#include <iostream>
#include <complex>
#include <cmath>

#include "cuda_runtime.h"
#include "cublas_v2.h"

// global variable used for texture memory optimization;
texture< int2, 1, cudaReadModeElementType > texRef;

//======================================================================
// kron_cuda_v1, y = (Ip @ Am @ Iq) * x;
//======================================================================
__global__ void kron_cuda_v1( const unsigned int m, const size_t s, const unsigned int n, const cuDoubleComplex *A, const unsigned int mat_i_idx_idx, const cuDoubleComplex *x, cuDoubleComplex *y )
{
  cuDoubleComplex     res;
  cuDoubleComplex     mid;
  int2                a1, a2;  

  unsigned int        k, x_idx;
  
  unsigned int        sex((unsigned int)s);
  unsigned int        x_line, A_tex_idx, i;
  
  extern __shared__   cuDoubleComplex x_shd[ ];
  
  //====================================================================
  //  determine the correspondence between thread and vector;
  //  k:      element index, unique for each CUDA thread;
  //  x_idx:  element index, unique in vector x;
  //====================================================================
//  k = blockDim.y * ( blockIdx.x + gridDim.x * blockIdx.y + gridDim.x * gridDim.y * blockIdx.z ) + threadIdx.y;
  
  k = blockIdx.z * gridDim.y + blockIdx.y;
  k = k * gridDim.x + blockIdx.x;
  k = k * blockDim.y + threadIdx.y;
  
  //  idx2 = s * n * ( k / n ) + k % n;
  //  idx2 = ( s - 1 ) * n * ( k / n ) + k;
  //  x_idx = ( sex - 1 ) * n * ( k / n ) + k + n * threadIdx.x;
  x_idx = k / n;
  x_idx = x_idx * (sex - 1) + threadIdx.x;
  x_idx = x_idx * n + k;
  
  x_line = blockDim.x * threadIdx.y;
  
  //====================================================================
  // copy x to the share memory x_shd;
  //====================================================================
  //  x_shd[ threadIdx.x + blockDim.x * threadIdx.y ] = x[ idx2 + n * threadIdx.x];
  x_shd[threadIdx.x + x_line] = x[x_idx];
  __syncthreads();
  
  //====================================================================
  // matrix multiplication using the global memory with texture memory;
  //====================================================================
//  A_tex_idx = (mat_i_idx_idx + threadIdx.x) * 2;
//  a1 = tex1Dfetch(texRef, A_tex_idx);
//  a2 = tex1Dfetch(texRef, A_tex_idx + 1);
//  mid.x = __hiloint2double(a1.y, a1.x);
//  mid.y = __hiloint2double(a2.y, a2.x);
//  res = cuCmul(x[x_idx - n * threadIdx.x], mid);
//  for (i = 1; i < s; i++)
//  {
//    A_tex_idx = (mat_i_idx_idx + threadIdx.x + sex * i) * 2;
//    a1 = tex1Dfetch(texRef, A_tex_idx);
//    a2 = tex1Dfetch(texRef, A_tex_idx + 1);
//    mid.x = __hiloint2double(a1.y, a1.x);
//    mid.y = __hiloint2double(a2.y, a2.x);
//    res = cuCadd(res, cuCmul(x[x_idx - n * (threadIdx.x - i)], mid));
//  }
  
  //====================================================================
  // matrix multiplication using the share memory without texture memory;
  //====================================================================
//  A_tex_idx = mat_i_idx_idx + threadIdx.x;
//  res = cuCmul(x_shd[x_line], A[A_tex_idx]);
//  for (i = 1; i < s; i++)
//  {
//    res = cuCadd(res, cuCmul(x_shd[i + x_line], A[A_tex_idx + sex * i]));
//  }
  
  //====================================================================
  // matrix multiplication using the share memory with texture memory;
  //====================================================================
  A_tex_idx = (mat_i_idx_idx + threadIdx.x) * 2;
  a1 = tex1Dfetch(texRef, A_tex_idx);
  a2 = tex1Dfetch(texRef, A_tex_idx + 1);
  mid.x = __hiloint2double(a1.y, a1.x);
  mid.y = __hiloint2double(a2.y, a2.x);
  res = cuCmul(x_shd[x_line], mid);
  for (i = 1; i < s; i++)
  {
    A_tex_idx = (mat_i_idx_idx + threadIdx.x + sex * i) * 2;
    a1 = tex1Dfetch(texRef, A_tex_idx);
    a2 = tex1Dfetch(texRef, A_tex_idx + 1);
    mid.x = __hiloint2double(a1.y, a1.x);
    mid.y = __hiloint2double(a2.y, a2.x);
    res = cuCadd(res, cuCmul(x_shd[i + x_line], mid));
  }
  
  //====================================================================
  // add to the result;
  //====================================================================
  y[x_idx] = res;
}
//======================================================================

//======================================================================
// kron_cuda_v3, y = y + coeff * (Ip @ Am @ Iq) * x;
//======================================================================
__global__ void kron_cuda_v3( const unsigned int m, const size_t s, const unsigned int n, const cuDoubleComplex *A, const unsigned int mat_i_idx_idx, const cuDoubleComplex *x, cuDoubleComplex *y, const cuDoubleComplex *coeff )
{
  cuDoubleComplex     res;
  cuDoubleComplex     mid;
  int2                a1, a2;  

  unsigned int        k, x_idx;
  
  unsigned int        sex((unsigned int)s);
  unsigned int        x_line, A_tex_idx, i;
  
  extern __shared__   cuDoubleComplex x_shd[ ];
  
  //====================================================================
  //  determine the correspondence between thread and vector;
  //  k:      element index, unique for each CUDA thread;
  //  x_idx:  element index, unique in vector x;
  //====================================================================
  //  k = blockDim.y * ( blockIdx.x + gridDim.x * blockIdx.y + gridDim.x * gridDim.y * blockIdx.z ) + threadIdx.y;
  k = blockIdx.z * gridDim.y + blockIdx.y;
  k = k * gridDim.x + blockIdx.x;
  k = k * blockDim.y + threadIdx.y;
  
  //  idx2 = s * n * ( k / n ) + k % n;
  //  idx2 = ( s - 1 ) * n * ( k / n ) + k;
  //  x_idx = ( s - 1 ) * n * ( k / n ) + k + n * threadIdx.x;
  x_idx = k / n;
  x_idx = x_idx * (sex - 1) + threadIdx.x;
  x_idx = x_idx * n + k;
  
  x_line = blockDim.x * threadIdx.y;
  
  //====================================================================
  // copy x to the share memory x_shd;
  //====================================================================
  //  x_shd[ threadIdx.x + blockDim.x * threadIdx.y ] = x[ idx2 + n * threadIdx.x];
  x_shd[threadIdx.x + x_line] = x[x_idx];
  __syncthreads();
  
  //====================================================================
  // matrix multiplication using the global memory with texture memory;
  //====================================================================
//  A_tex_idx = (mat_i_idx_idx + threadIdx.x) * 2;
//  a1 = tex1Dfetch(texRef, A_tex_idx);
//  a2 = tex1Dfetch(texRef, A_tex_idx + 1);
//  mid.x = __hiloint2double(a1.y, a1.x);
//  mid.y = __hiloint2double(a2.y, a2.x);
//  res = cuCmul(x[x_idx - n * threadIdx.x], mid);
//  for (i = 1; i < s; i++)
//  {
//    A_tex_idx = (mat_i_idx_idx + threadIdx.x + sex * i) * 2;
//    a1 = tex1Dfetch(texRef, A_tex_idx);
//    a2 = tex1Dfetch(texRef, A_tex_idx + 1);
//    mid.x = __hiloint2double(a1.y, a1.x);
//    mid.y = __hiloint2double(a2.y, a2.x);
//    res = cuCadd(res, cuCmul(x[x_idx - n * (threadIdx.x - i)], mid));
//  }
  
  //====================================================================
  // matrix multiplication using the share memory without texture memory;
  //====================================================================
//  A_tex_idx = mat_i_idx_idx + threadIdx.x;
//  res = cuCmul(x_shd[x_line], A[A_tex_idx]);
//  for (i = 1; i < s; i++)
//  {
//    res = cuCadd(res, cuCmul(x_shd[i + x_line], A[A_tex_idx + sex * i]));
//  }
  
  //====================================================================
  // matrix multiplication using the share memory with texture memory;
  //====================================================================
  A_tex_idx = (mat_i_idx_idx + threadIdx.x) * 2;
  a1 = tex1Dfetch(texRef, A_tex_idx);
  a2 = tex1Dfetch(texRef, A_tex_idx + 1);
  mid.x = __hiloint2double(a1.y, a1.x);
  mid.y = __hiloint2double(a2.y, a2.x);
  res = cuCmul(x_shd[x_line], mid);
  for (i = 1; i < s; i++)
  {
    A_tex_idx = (mat_i_idx_idx + threadIdx.x + sex * i) * 2;
    a1 = tex1Dfetch(texRef, A_tex_idx);
    a2 = tex1Dfetch(texRef, A_tex_idx + 1);
    mid.x = __hiloint2double(a1.y, a1.x);
    mid.y = __hiloint2double(a2.y, a2.x);
    res = cuCadd(res, cuCmul(x_shd[i + x_line], mid));
  }
  
  //====================================================================
  // add to the result;
  //====================================================================
  res = cuCmul(res, *coeff);
  y[x_idx] = cuCadd(y[x_idx], res);
}
//======================================================================

__global__ void vecrzt_kernel( cuDoubleComplex *x )
{
  size_t  idx;
  
  idx = blockDim.x * ( blockIdx.x + gridDim.x * blockIdx.y + gridDim.x * gridDim.y * blockIdx.z ) + threadIdx.x;
  
  x[ idx ] = make_cuDoubleComplex(0.0,0.0);
}

#include "include/math/get_grid.h"
//======================================================================
// kron_cuda_v4, y = y + coeff * (Ip @ Am @ Iq @ Bn @ Ir) * x;
//======================================================================
__global__ void kron_cuda_v4(const unsigned int p, const unsigned int m, const unsigned int q, const size_t n, const unsigned int r, const cuDoubleComplex *A, const unsigned int A_idx, const unsigned int B_idx, const cuDoubleComplex *coeff, const cuDoubleComplex *x, cuDoubleComplex *y)
{
  extern __shared__   cuDoubleComplex x_shd[ ];
  cuDoubleComplex     res;
  
  unsigned int        k;
  unsigned int        iq, ir, ip, im, in;
  unsigned int        x_idx;
  
  unsigned int        x_line, Aidx, Bidx, A_tex_idx, B_tex_idx;
  unsigned int        nex((unsigned int)n);
  
  cuDoubleComplex     mul;
  cuDoubleComplex     mid;
  int2                a1, a2;
  
  unsigned int        i, j;
  
  //====================================================================
  //  determine the correspondence between thread and vector;
  //  q * r * p = blockDim.y * gridDim.x * gridDim.y * gridDim.z;
  //  0 <= k < qrp;
  //  
  //  get the numbering of the block;
  //  k:            line index, unique for every line of CUDA block;
  //                determine ip, iq, ir;
  //                k = (blockIdx.x + blockIdx.y * gridDim.x 
  //                    + blockIdx.z * gridDim.y * gridDim.x ) 
  //                    * blockDim.y + threadIdx.y;
  //                  = r * q * ip + r * iq + ir;
  //  threadIdx.x:  row index, unique for every column of CUDA block;
  //                determine im, in;
  //                threadIdx.x = n * im + in;
  //  t, s:         perfect shuffle permutation index;
  //                s = P_{q,m}(t), P denotes the permutation;
  //                t = m * iq + im, numbering in share memory;
  //                s = q * im + iq, numbering in global memory;
  //  x_idx:        element index, unique in vector x;
  //                x_idx = ip * q * m * n * r + iqs * m * n * r 
  //                        + ims * n * r + in * r + ir;
  //====================================================================
  
  //  k = threadIdx.y + (blockIdx.x + blockIdx.y * gridDim.x + blockIdx.z * gridDim.y * gridDim.x ) * blockDim.y;
  k = blockIdx.z * gridDim.y + blockIdx.y;
  k = k * gridDim.x + blockIdx.x;
  k = k * blockDim.y + threadIdx.y;
  
  //  ir = k % r;
  //  iq = (k / r) % q;
  //  ip = k / (r * q);
//  ip = k / (r * q);
//  ir = k - r * q * ip;
//  iq = ir / r;
//  ir -= r * iq;
  iq = r * q;
  ip = k / iq;
  ir = k - iq * ip;
  iq = ir / r;
  ir -= r * iq;
  
  // threadIdx.x = n * im + in;
  im = threadIdx.x / nex;
  in = threadIdx.x - nex * im; 
  
  //  s = q * im + iq;
  //  iqs = s / m;
  //  ims = s % m;
  //  x_idx = ip * q * m * n * r + iqs * m * n * r + ims * n * r + in * r + ir;
  x_idx = (((ip * m + im) * q + iq) * nex + in) * r + ir;
  
  x_line = blockDim.x * threadIdx.y;
  
  //====================================================================
  // copy x to the share memory x_shd;
  // x_shd is in row-major order, x is in column-major order;
  //====================================================================
  x_shd[threadIdx.x + x_line] = x[x_idx];
  __syncthreads();
  
  //====================================================================
  // y_shd[threadIdx.y, threadIdx.x]  = x_shd[threadIdx.y,:] 
  //                                  * kron(AT, BT)[:,threadIdx.x];
  // AT (or BT) is the transpose of A (or B);
  // row index of x_shd is threadIdx.y;
  // column indexes of AT and BT are "im" and "in";
  // row indexes of A and B are "im" and "in";
  //====================================================================
  Bidx = B_idx + in;
  Aidx = A_idx + im;
  
  res = make_cuDoubleComplex(0.0, 0.0);
  
  //====================================================================
  // matrix multiplication using the global memory with texture memory;
  //====================================================================
//  for (i = 0; i < m; i++)
//  {
//    mul = make_cuDoubleComplex(0.0, 0.0);
//    x_line = ((ip * m + i) * q + iq) * nex;
//    for (j = 0; j < n; j++)
//    {
//      B_tex_idx = (Bidx + nex * j) * 2;
//      // get B[in + n * j];
//      a1 = tex1Dfetch(texRef, B_tex_idx);
//      a2 = tex1Dfetch(texRef, B_tex_idx + 1);
//      mid.x = __hiloint2double(a1.y, a1.x);
//      mid.y = __hiloint2double(a2.y, a2.x);
//      mul = cuCadd(mul, cuCmul(x[(x_line + j) * r + ir], mid));
//    }
//    A_tex_idx = (Aidx + m * i) * 2;
//    // get A[im + m * i];
//    a1 = tex1Dfetch(texRef, A_tex_idx);
//    a2 = tex1Dfetch(texRef, A_tex_idx + 1);
//    mid.x = __hiloint2double(a1.y, a1.x);
//    mid.y = __hiloint2double(a2.y, a2.x);
//    res = cuCadd(res, cuCmul(mul, mid));
//  }
  
  //====================================================================
  // matrix multiplication using the share memory without texture memory;
  //====================================================================
//  for (i = 0; i < m; i++)
//  {
//    mul = make_cuDoubleComplex(0.0, 0.0);
//    ip = n * i + x_line;
//    for (j = 0; j < n; j++)
//    {
//      B_tex_idx = Bidx + nex * j;
//      mul = cuCadd(mul, cuCmul(x_shd[ip + j], A[B_tex_idx]));
//    }
//    A_tex_idx = Aidx + m * i;
//    res = cuCadd(res, cuCmul(mul, A[A_tex_idx]));
//  }
  
  //====================================================================
  // matrix multiplication using the share memory with texture memory;
  //====================================================================
  for (i = 0; i < m; i++)
  {
    mul = make_cuDoubleComplex(0.0, 0.0);
    for (j = 0; j < n; j++)
    {
      B_tex_idx = (Bidx + nex * j) * 2;
      // get B[in + n * j];
      a1 = tex1Dfetch(texRef, B_tex_idx);
      a2 = tex1Dfetch(texRef, B_tex_idx + 1);
      mid.x = __hiloint2double(a1.y, a1.x);
      mid.y = __hiloint2double(a2.y, a2.x);
      mul = cuCadd(mul, cuCmul(x_shd[n * i + j + x_line], mid));
    }
    A_tex_idx = (Aidx + m * i) * 2;
    // get A[im + m * i];
    a1 = tex1Dfetch(texRef, A_tex_idx);
    a2 = tex1Dfetch(texRef, A_tex_idx + 1);
    mid.x = __hiloint2double(a1.y, a1.x);
    mid.y = __hiloint2double(a2.y, a2.x);
    res = cuCadd(res, cuCmul(mul, mid));
  }
  
  //====================================================================
  // add to the result;
  //====================================================================
  res = cuCmul(res, *coeff);
  y[x_idx] = cuCadd(y[x_idx], res);
}
//======================================================================

void hamvec_cuda3( cublasHandle_t cublas_handle, int nspin, int nTerm, std::complex<double> *coeff_lst_zplx, size_t *nbody_lst, size_t *pos_i_idx, size_t *pos_i_lst, size_t *dim_i_lst, size_t *mat_i_idx, std::complex<double> *dev_mat_i_lst, size_t vlen, size_t *nspin_dim, size_t *nspin_m_lst, size_t *nspin_n_lst, std::complex<double> *dev_v, std::complex<double> *dev_w, std::complex<double> *dev_w_med, std::complex<double> *dev_coeff_lst_zplx, size_t maxThreadsPerBlock, size_t *maxGridSize )
{
/*
  Calculate the action of a Hamiltonian on a state vector.
  The Hamiltonian can be decomposed into many terms, where each term 
  consists of 
  
  input:
    
    cublas_handle,  handle of cublas;
    nspin,          number of bodies;
    nTerm,          number of interactions in the Hamiltonian;
    coeff_lst_zplx, list of the coefficient in each interaction;
    nbody_lst,      list of the number of bodies in each interaction;
    pos_i_idx,      list of the index of the position list of the body 
                    in each interaction;
    pos_i_lst,      list of the positions of bodies in each interaction;
    dim_i_lst,      list of the dimension of operator of each body in 
                    each interaction;
    mat_i_idx,      list of the index of the operator list of the body
                    in each interaction;
    dev_mat_i_lst,  list of the operators of bodies in each interaction; 
    vlen,           dimension of the state vector;
    nspin_dim,      list of the dimension of each body;
    nspin_m_lst,    list of the dimension of first m bodies;
    nspin_n_lst,    list of the dimension of last n bodies;
    dev_v,          input state vector;
    dev_w_med,      intermediate state vector;
    dev_coeff_lst_zplx, the same as "coeff_lst_zplx";
    maxThreadsPerBlock,   max block size for CUDA;
    maxGridSize,    grid size limit for CUDA;
  
  output:
    
    dev_w,          output state vector, after the Hamiltonian acting on
                    the input state;
*/
  
  size_t  nT, nbody, nb;
  size_t  idx, pos_i, dim_i;
  std::complex<double>  coeff;
  
  size_t  m, n;
  
  dim3    grid_dim, block_dim;
  size_t  dimex1, dimex2;
  
  // optimization by texture memory, bind the operator list;
  cudaBindTexture( 0, texRef, dev_mat_i_lst );
  
  size_t          idx1, idx2, pos1, pos2, i;
  size_t          p, q, r;
  unsigned long   nthreads;
  unsigned int    nthreads_max_per_block;
  dim3            blim(1024,1024,64);
  dim3            glim(2147483647,65535,65535);
//  unsigned long   bnlimit(256);
  dim3            grid_limit, block_limit, grid_size, block_size;
  //====================================================================
  // initialization of w;
  //====================================================================
  
  // set block size, 2D;
  block_dim.z = 1;
  block_dim.y = 1;
  dimex1 = maxThreadsPerBlock;
  dimex2 = vlen;
  while ( (dimex2 % dimex1) != 0 ) dimex1--;
  block_dim.x = dimex1;
  dimex2 /= dimex1;
  if ( dimex2 <= maxGridSize[0] * maxGridSize[1] * maxGridSize[2] )
  {
    dimex1 = dimex2 / ( maxGridSize[0] * maxGridSize[1] );
    if ( ( dimex2 % ( maxGridSize[0] * maxGridSize[1] ) ) > 0 )
      dimex1++;
    while ( (dimex2 % dimex1) != 0 ) dimex1++;
    grid_dim.z = dimex1;
    dimex2 = dimex2 / grid_dim.z;
    dimex1 = dimex2 / maxGridSize[0];
    if ( ( dimex2 % maxGridSize[0] ) > 0)
      dimex1++;
    while ( (dimex2 % dimex1) != 0 ) dimex1++;
    grid_dim.y = dimex1;
    grid_dim.x = dimex2 / grid_dim.y;
  }
  else
  {
    std::cout << "block number exceeds limit." << std::endl;
    return;
  }
  
  // reset the output vector;
  vecrzt_kernel<<< grid_dim, block_dim >>>( (cuDoubleComplex*)dev_w );
  
  //====================================================================
  // action of the Hamiltonian on the input state;
  //====================================================================
  
//  // obtain grid;
//  m = dim_i_lst[0];
//  n = dim_i_lst[0];
//  grid_limit = glim;
//  block_limit = blim; 
//  nthreads = vlen / (m * n);
//  nthreads_max_per_block = maxThreadsPerBlock / (m * n);
//  block_limit.x = 1;
//  block_limit.z = 1;
//  get_grid(nthreads, nthreads_max_per_block, grid_limit, block_limit, &grid_size, &block_size );
//  block_size.x = m * n;
  
  // loop over each term of the Hamiltonian;
  for ( nT = 0; nT < nTerm; nT++ )
  {
    coeff = coeff_lst_zplx[ nT ];
    nbody = nbody_lst[ nT ];
    
    if ( abs(coeff) == 0 ) continue;
    
    // loop over each body in each term of the Hamiltonian;
    for ( nb = 0; nb < nbody; nb++ )
    {
      idx = pos_i_idx[ nT ] + nbody - 1 - nb;
      pos_i = pos_i_lst[ idx ];
      dim_i = dim_i_lst[ idx ];
      
      m = nspin_m_lst[ pos_i ];
      n = nspin_n_lst[ pos_i ];
      
      // set block size, 2D;
      block_dim.z = 1;
      block_dim.x = dim_i;
      dimex1 = maxThreadsPerBlock / block_dim.x;
      dimex2 = m * n;
      while ( (dimex2 % dimex1) != 0 ) dimex1--;
      block_dim.y = dimex1;
      // set grid size, 1D-3D;
      dimex2 /= dimex1;// grid size is m * n / block_dim.y;
      if ( dimex2 <= maxGridSize[0] * maxGridSize[1] * maxGridSize[2] )
      {
        dimex1 = dimex2 / ( maxGridSize[0] * maxGridSize[1] );
        if ( ( dimex2 % ( maxGridSize[0] * maxGridSize[1] ) ) > 0 )
          dimex1++;
        while ( (dimex2 % dimex1) != 0 ) dimex1++;
        grid_dim.z = dimex1;
        dimex2 = dimex2 / grid_dim.z;
        dimex1 = dimex2 / maxGridSize[0];
        if ( ( dimex2 % maxGridSize[0] ) > 0)
          dimex1++;
        while ( (dimex2 % dimex1) != 0 ) dimex1++;
        grid_dim.y = dimex1;
        grid_dim.x = dimex2 / grid_dim.y;
      }
      else
      {
        std::cout << "block number exceeds limit." << std::endl;
        return;
      }
      
      if ( nbody == 1 )
      {
//        // for debug;
//        std::cout << "run for nbody == 1. " << std::endl;
        
        kron_cuda_v3<<< grid_dim, block_dim, block_dim.x * block_dim.y * sizeof(std::complex<double>) >>>( m, dim_i, n, (cuDoubleComplex*)dev_mat_i_lst, mat_i_idx[idx], (cuDoubleComplex*)dev_v, (cuDoubleComplex*)dev_w, (cuDoubleComplex*)&dev_coeff_lst_zplx[nT] );
      }
      else if (nbody == 2)
      {
        // new 2to1 mat-vec;
        idx2 = pos_i_idx[ nT ] + nbody - 1;
        idx1 = idx2 - 1;
        n = dim_i_lst[ idx2 ];
        m = dim_i_lst[ idx1 ];
        pos2 = pos_i_lst[ idx2 ];
        pos1 = pos_i_lst[ idx1 ];
        
        p = nspin_m_lst[ pos1 ];
        r = nspin_n_lst[ pos2 ];
        q = 1;
        for (i = ( pos1 < pos2 ? pos1 : pos2) + 1; i < ( pos1 < pos2 ? pos2 : pos1); i++)
          q *= nspin_dim[i];
        
        // obtain grid;
        grid_limit = glim;
        block_limit = blim; 
        nthreads = vlen / (m * n);
        nthreads_max_per_block = maxThreadsPerBlock / (m * n);
        block_limit.x = 1;
        block_limit.z = 1;
        get_grid(nthreads, nthreads_max_per_block, grid_limit, block_limit, &grid_size, &block_size );
        block_size.x = m * n;
        
//        // for debug;
//        std::cout << "run for nbody == 2. " << std::endl;
//        std::cout << "p, m, q, n, r = " << p << ", " << m << ", " << q << ", " << n << ", " << r << ", " << std::endl;
//        std::cout << "grid_size.xyz = " << grid_size.x << ", " << grid_size.y << ", " << grid_size.z << std::endl;
//        std::cout << "block_size.xyz = " << block_size.x << ", " << block_size.y << ", " << block_size.z << std::endl;
////        return;
//        std::cout << "A mat idx = " << mat_i_idx[idx1] << std::endl;
//        std::cout << "B mat idx = " << mat_i_idx[idx2] << std::endl;
//        std::complex<double>  *matA, *matB;
//        cudaError_t           cuda_status;
//        int i;
//        matA = new std::complex<double> [4];
//        matB = new std::complex<double> [4];
//        cuda_status = cudaMemcpy( matA, &dev_mat_i_lst[mat_i_idx[idx1]], ( 4 * sizeof( matA[0] ) ), cudaMemcpyDeviceToHost );
//        if (cuda_status != cudaSuccess)
//          std::cout << "matA memcpy failed!" << std::endl;
//        cuda_status = cudaMemcpy( matB, &dev_mat_i_lst[mat_i_idx[idx2]], ( 4 * sizeof( matB[0] ) ), cudaMemcpyDeviceToHost );
//        if (cuda_status != cudaSuccess)
//          std::cout << "matB memcpy failed!" << std::endl;
//        std::cout << "matA: " << std::endl;
//        for (i = 0; i < 4; i++)
//          std::cout << i << ", " << matA[i] << std::endl;
//        std::cout << "matB: " << std::endl;
//        for (i = 0; i < 4; i++)
//          std::cout << i << ", " << matB[i] << std::endl;
//        delete[] matA;
//        delete[] matB;
        
        // launch kernel;
        kron_cuda_v4<<< grid_size, block_size, block_size.x * block_size.y * sizeof(std::complex<double>) >>>( p, m, q, n, r, (cuDoubleComplex*)dev_mat_i_lst, mat_i_idx[idx1], mat_i_idx[idx2], (cuDoubleComplex*)&dev_coeff_lst_zplx[nT], (cuDoubleComplex*)dev_v, (cuDoubleComplex*)dev_w );
        break;
      }
      else
      {
        // for debug;
        std::cout << "nbody run! skip calculation!" << std::endl;
        return;
        
        if ( nb == 0 )
          kron_cuda_v1<<< grid_dim, block_dim, block_dim.x * block_dim.y * sizeof(std::complex<double>) >>>( m, dim_i, n, (cuDoubleComplex*)dev_mat_i_lst, mat_i_idx[idx], (cuDoubleComplex*)dev_v, (cuDoubleComplex*)dev_w_med );
        else if ( nb == nbody - 1 )
          kron_cuda_v3<<< grid_dim, block_dim, block_dim.x * block_dim.y * sizeof(std::complex<double>) >>>( m, dim_i, n, (cuDoubleComplex*)dev_mat_i_lst, mat_i_idx[idx], (cuDoubleComplex*)dev_w_med, (cuDoubleComplex*)dev_w, (cuDoubleComplex*)&dev_coeff_lst_zplx[nT] );
        else
        {
          std::cout << "error using, nb should be 1 or 2. " << std::endl;
          return;
//          kron_cuda_v2<<< grid_dim, block_dim, block_dim.x * block_dim.y * sizeof(std::complex<double>) >>>( m, dim_i, n, (cuDoubleComplex*)dev_mat_i_lst, mat_i_idx[idx], (cuDoubleComplex*)dev_w_med );
//          kron_cuda_v1<<< grid_dim, block_dim, block_dim.x * block_dim.y * sizeof(std::complex<double>) >>>( m, dim_i, n, (cuDoubleComplex*)dev_mat_i_lst, mat_i_idx[idx], (cuDoubleComplex*)dev_w_med, (cuDoubleComplex*)dev_w_med );
        }
      }
      
    }
    
  }
  
  // optimization by texture memory, release the bind;
  cudaUnbindTexture( texRef );
  
//  // for debug;
//  std::cout << "grid_size.xyz = " << grid_size.x << ", " << grid_size.y << ", " << grid_size.z << std::endl;
//  std::cout << "block_size.xyz = " << block_size.x << ", " << block_size.y << ", " << block_size.z << std::endl;
  
}

//======================================================================
// interface for FORTRAN;
//======================================================================

// declaration of function handle;
#define HAMVEC_CUDA3       hamvec_cuda3_
#define HAMVEC_CUDA3_INIT  hamvec_cuda3_init_
#define HAMVEC_CUDA3_TERM  hamvec_cuda3_term_

// delcaration of function;
#if defined(__cplusplus)
extern "C" {
#endif /* __cplusplus */

void HAMVEC_CUDA3( size_t *cublas_handle_ptr, int *nspin_ptr, int *nTerm_ptr, std::complex<double> *coeff_lst_zplx_ptr, size_t *nbody_lst_ptr, size_t *pos_i_idx_ptr, size_t *pos_i_lst_ptr, size_t *dim_i_lst_ptr, size_t *mat_i_idx_ptr, size_t *dev_mat_i_lst_ptr, size_t *ham_dim_ptr, size_t *nspin_dim_ptr, size_t *nspin_m_lst_ptr, size_t *nspin_n_lst_ptr, size_t *dev_v_ptr, size_t *dev_w_ptr, size_t *dev_w_med_ptr, size_t *dev_coeff_lst_zplx_ptr, size_t *maxThreadsPerBlock_ptr, size_t *maxGridSize_ptr );

void HAMVEC_CUDA3_INIT( size_t *cublas_handle_ptr );

void HAMVEC_CUDA3_TERM( size_t *cublas_handle_ptr );

#if defined(__cplusplus)
}
#endif /* __cplusplus */

// interface of function;
void HAMVEC_CUDA3( size_t *cublas_handle_ptr, int *nspin_ptr, int *nTerm_ptr, std::complex<double> *coeff_lst_zplx_ptr, size_t *nbody_lst_ptr, size_t *pos_i_idx_ptr, size_t *pos_i_lst_ptr, size_t *dim_i_lst_ptr, size_t *mat_i_idx_ptr, size_t *dev_mat_i_lst_ptr, size_t *ham_dim_ptr, size_t *nspin_dim_ptr, size_t *nspin_m_lst_ptr, size_t *nspin_n_lst_ptr, size_t *dev_v_ptr, size_t *dev_w_ptr, size_t *dev_w_med_ptr, size_t *dev_coeff_lst_zplx_ptr, size_t *maxThreadsPerBlock_ptr, size_t *maxGridSize_ptr )
{
  cublasHandle_t        cublas_handle   = (cublasHandle_t)*cublas_handle_ptr;
  int                   nspin           = *nspin_ptr;
  int                   nTerm           = *nTerm_ptr;
  std::complex<double>  *coeff_lst_zplx = coeff_lst_zplx_ptr;
  size_t                *nbody_lst      = nbody_lst_ptr;
  size_t                *pos_i_idx      = pos_i_idx_ptr;
  size_t                *pos_i_lst      = pos_i_lst_ptr;
  size_t                *dim_i_lst      = dim_i_lst_ptr;
  size_t                *mat_i_idx      = mat_i_idx_ptr;
  std::complex<double>  *dev_mat_i_lst  = (std::complex<double>*)(*dev_mat_i_lst_ptr);
  size_t                vlen            = *ham_dim_ptr;
  size_t                *nspin_dim      = nspin_dim_ptr;
  size_t                *nspin_m_lst    = nspin_m_lst_ptr;
  size_t                *nspin_n_lst    = nspin_n_lst_ptr;
  std::complex<double>  *dev_v          = (std::complex<double>*)(*dev_v_ptr);
  std::complex<double>  *dev_w          = (std::complex<double>*)(*dev_w_ptr);
  std::complex<double>  *dev_w_med      = (std::complex<double>*)(*dev_w_med_ptr);
  std::complex<double>  *dev_coeff_lst_zplx = (std::complex<double>*)(*dev_coeff_lst_zplx_ptr);
  size_t            maxThreadsPerBlock  = *maxThreadsPerBlock_ptr;
  size_t                *maxGridSize    = maxGridSize_ptr;
  
  hamvec_cuda3( cublas_handle, nspin, nTerm, coeff_lst_zplx, nbody_lst, pos_i_idx, pos_i_lst, dim_i_lst, mat_i_idx, dev_mat_i_lst, vlen, nspin_dim, nspin_m_lst, nspin_n_lst, dev_v, dev_w, dev_w_med, dev_coeff_lst_zplx, maxThreadsPerBlock, maxGridSize );
}

void HAMVEC_CUDA3_INIT( size_t *cublas_handle_ptr )
{
  // initialization of cublas handle for FORTRAN;
  cublasHandle_t cublas_handle;
  cublasCreate( &cublas_handle );
  *cublas_handle_ptr = (size_t)cublas_handle;
}

void HAMVEC_CUDA3_TERM( size_t *cublas_handle_ptr )
{
  // termination of cublas handle for FORTRAN;
  cublasHandle_t cublas_handle;
  cublas_handle = (cublasHandle_t)*cublas_handle_ptr;
  cublasDestroy( cublas_handle );
}

