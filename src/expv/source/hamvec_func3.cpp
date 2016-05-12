// |w> = |w> + H * |v>;
// H = sum(H_i) + sum(H_ij);
// H_i = h_i * S_(i,xyz) or h_ij * S_(i,xyz) * S_(j,xyz);
// w_med = H_ij * w_med 

#include <iostream>
#include <complex>
#include <cmath>

#include "omp.h"
#define MKL_Complex16 std::complex<double>
#include "mkl.h"

//======================================================================
// kron_func_v2.cpp
//======================================================================
void kron_func_v2( size_t m, size_t s, size_t n, std::complex<double> *A, std::complex<double> *x )
{
  std::complex<double> zero(0.0,0.0);
  std::complex<double> one(1.0,0.0);
  
  size_t inc1, inc2, inc3;
  
  inc1 = s * n;
  inc2 = n;
  inc3 = s;
  
  #pragma omp parallel
  {
  size_t k;
  size_t i, j, p, q, idx1, idx2;
  std::complex<double> res;
  std::complex<double> y[s];
  #pragma omp for
  for ( k = 0; k < m*n; k++ )
  {
    i = k / n;
    j = k % n;
//  for ( i = 0; i < m; i++ )
//  {
    idx1 = inc1 * i;
//    for ( j = 0; j < n; j++ )
//    {
      idx2 = idx1 + j;
      
      for ( p = 0; p < s; p++ )
      {
//        idx3 = s * p;
        res = zero;
        for ( q = 0; q < s; q++ )
        {
          res += x[ idx2 + inc2 * q ] * A[ p + inc3 * q ];
        }
        y[ p ] = res;
      }
      for ( p = 0; p < s; p++ )
        x[ idx2 + inc2 * p] = y[p];
      
//    }
//  }
  }
  } //  end of pragma omp parallel;
  
}
//======================================================================

//======================================================================
// kron_func_v3
//======================================================================
void kron_func_v3( size_t m, size_t s, size_t n, std::complex<double> *A, std::complex<double> *x )
{
  std::complex<double> zero(0.0,0.0);
  
  size_t k, p, q, idx2;
  std::complex<double> res;
  std::complex<double> x_shd[s];
  
//  #pragma omp parallel for private( res, x_shd, k, idx2, p, q ) firstprivate( A, x, m, n, s, zero )
  for ( k = 0; k < m*n; k++ )
  {
    idx2 = ( s - 1 ) * n * ( k / n ) + k;
    for ( q = 0; q < s; q++ ) x_shd[ q ] = x[ idx2 + n * q ];
    for ( p = 0; p < s; p++ )
    {
      res = zero;
      for ( q = 0; q < s; q++ )
      {
        res += x_shd[ q ] * A[ p + s * q ];
      }
      x[ idx2 + n * p ] = res;
    }
  }
}
//======================================================================

//======================================================================
// kron_func_v32
//======================================================================
void kron_func_v32( size_t m, size_t s, size_t n, std::complex<double> *A, std::complex<double> *x )
{
  std::complex<double> zero(0.0,0.0);
  
  size_t k, p, q, idx2;
  std::complex<double> res;
  std::complex<double> x_shd[s];
  
  for ( k = 0; k < m*n; k++ )
  {
    idx2 = ( s - 1 ) * n * ( k / n ) + k;
    for ( q = 0; q < s; q++ ) x_shd[ q ] = x[ idx2 + n * q ];
    #pragma omp parallel for private( res, p, q )
    for ( p = 0; p < s; p++ )
    {
      res = zero;
      for ( q = 0; q < s; q++ )
      {
        res += x_shd[ q ] * A[ p + s * q ];
      }
      x[ idx2 + n * p ] = res;
    }
  }
}
//======================================================================

//======================================================================
// kron_func_v4
//======================================================================
void kron_func_v4( size_t m, size_t s, size_t n, std::complex<double> *A, std::complex<double> *x, int bdy )
{
  std::complex<double> zero(0.0,0.0), one(1.0,0.0);
  
  size_t k, p, q, idx2;
  std::complex<double> res;
  
  int bdx = s, bds = s * bdy;
  size_t  idy[ bdy ];
  std::complex<double> x_shd[ bds ], y_shd[ bds ];
  
//  #pragma omp parallel for private( k, idy, x_shd, y_shd, p, q ) firstprivate( A, x, m, n, s, zero, one, bds, bdy )
  for ( k = 0; k < m*n; k += bdy )
  {
    for ( p = 0; p < bdy; p++ )
    {
      idy[ p ] = s * n * ( ( k + p ) / n ) + ( k + p ) % n;
      for ( q = 0; q < s; q++ )
        x_shd[ p + bdy * q ] = x[ idy[ p ] + n * q ];
    }
//    for ( p = 0; p < bds; p++ ) y_shd[ p ] = zero;
    cblas_zgemm( CblasColMajor, CblasNoTrans, CblasTrans, bdy, s, s, &one, x_shd, bdy, A, s, &zero, y_shd, bdy );
    for ( p = 0; p < bdy; p++ )
    {
      for ( q = 0; q < s; q++ )
        x[ idy[ p ] + n * q ] = y_shd[ p + bdy * q ];
    }
  }
}
//======================================================================

//======================================================================
// kron_func_v5
//======================================================================
void kron_func_v5( size_t m, size_t s, size_t n, std::complex<double> *A, std::complex<double> *x, int bdy )
{
  std::complex<double> zero(0.0,0.0), one(1.0,0.0);
  
  size_t k, p, q, l;
  std::complex<double> res;
  
  int bdx = s, bds = s * bdy;
  size_t  idy[ bdy ];
  std::complex<double> x_shd[ bds ], y_shd[ bds ];
  
//  #pragma omp parallel for private( k, p, q, l ) firstprivate( A, x, m, n, s, bdy, idy, x_shd, y_shd, one, zero )
  for ( k = 0; k < m*n; k += bdy )
  {
    for ( p = 0; p < bdy; p++ )
    {
      idy[ p ] = s * n * ( ( k + p ) / n ) + ( k + p ) % n;
      for ( q = 0; q < s; q++ )
        x_shd[ p + bdy * q ] = x[ idy[ p ] + n * q ];
    }
    
    // matrix-matrix multiplication;
//    cblas_zgemm( CblasColMajor, CblasNoTrans, CblasTrans, bdy, s, s, &one, x_shd, bdy, A, s, &zero, y_shd, bdy );
    for ( p = 0; p < bds; p++ ) y_shd[ p ] = zero;
    for ( q = 0; q < s; q++ )
      for ( l = 0; l < s; l++ )
        for ( p = 0; p < bdy; p++ )
          y_shd[ p + bdy * q ] += x_shd[ p + bdy * l ] * A[ s * l + q ];
    
    for ( p = 0; p < bdy; p++ )
    {
      for ( q = 0; q < s; q++ )
        x[ idy[ p ] + n * q ] = y_shd[ p + bdy * q ];
    }
  }
}
//======================================================================

//======================================================================
// kron_func_v6
//======================================================================
void kron_func_v6( size_t m, size_t s, size_t n, std::complex<double> *A, std::complex<double> *x, int bdy )
{
  std::complex<double> zero(0.0,0.0), one(1.0,0.0);
  
  size_t k, p, q, l;
  std::complex<double> res;
  
  int bdx = s, bds = s * bdy;
  size_t  idy[ bdy ];
  std::complex<double> x_shd[ bds ], y_shd[ bds ];
  
//  #pragma omp parallel for private( k, p, q, l ) firstprivate( A, x, m, n, s, bdy, idy, x_shd, y_shd, one, zero )
//  #pragma omp parallel private( k, p, q, l, x_shd, y_shd, idy )
  {
//  #pragma omp for
  for ( k = 0; k < m*n; k += bdy )
  {
    for ( p = 0; p < bdy; p++ )
    {
      idy[ p ] = s * n * ( ( k + p ) / n ) + ( k + p ) % n;
      for ( q = 0; q < s; q++ )
        x_shd[ p + bdy * q ] = x[ idy[ p ] + n * q ];
    }
    
    // matrix-matrix multiplication;
//    cblas_zgemm( CblasColMajor, CblasNoTrans, CblasTrans, bdy, s, s, &one, x_shd, bdy, A, s, &zero, y_shd, bdy );
    for ( p = 0; p < bds; p++ ) y_shd[ p ] = zero;
    for ( q = 0; q < s; q++ )
      for ( l = 0; l < s; l++ )
        for ( p = 0; p < bdy; p++ )
          y_shd[ p + bdy * q ] += x_shd[ p + bdy * l ] * A[ s * l + q ];
    
    for ( p = 0; p < bdy; p++ )
    {
      for ( q = 0; q < s; q++ )
        x[ idy[ p ] + n * q ] = y_shd[ p + bdy * q ];
    }
  }
  }// end of omp parallel;
}
//======================================================================

//======================================================================
// kron_func_v7
//======================================================================
void kron_func_v7( size_t m, size_t s, size_t n, std::complex<double> *A, std::complex<double> *x, int bdy )
{
  std::complex<double> zero(0.0,0.0), one(1.0,0.0);
  
  size_t k, p, q, l;
  std::complex<double> res;
  
  int bdx = s, bds = s * bdy;
  size_t  idy[ bdy ];
  std::complex<double> x_shd[ bds ], y_shd[ bds ];
  
  for ( k = 0; k < m*n; k += bdy )
  {
    for ( p = 0; p < bdy; p++ )
    {
      idy[ p ] = s * n * ( ( k + p ) / n ) + ( k + p ) % n;
      for ( q = 0; q < s; q++ )
        x_shd[ p + bdy * q ] = x[ idy[ p ] + n * q ];
    }
    
    // matrix-matrix multiplication;
//    cblas_zgemm( CblasColMajor, CblasNoTrans, CblasTrans, bdy, s, s, &one, x_shd, bdy, A, s, &zero, y_shd, bdy );
//    for ( p = 0; p < bds; p++ ) y_shd[ p ] = zero;
//    for ( q = 0; q < s; q++ )
//      for ( l = 0; l < s; l++ )
//        for ( p = 0; p < bdy; p++ )
//          y_shd[ p + bdy * q ] += x_shd[ p + bdy * l ] * A[ s * l + q ];
    
//    for ( q = 0; q < s; q++ )
//      for ( p = 0; p < bdy; p++ )
//        x[ idy[ p ] + n * q ] = y_shd[ p + bdy * q ];
    
//    #pragma omp parallel for private( res, p, q, l )
    for ( p = 0; p < bdy; p++ )
    {
      for ( q = 0; q < s; q++ )
      {
        res = zero;
        for ( l = 0; l < s; l++ )
          res += x_shd[ p + bdy * l ] * A[ s * l + q ];
//        *( x + ( *( idy + p ) + n * q ) ) = res;
        x[ idy[ p ] + n * q ] = res;
      }
    }
    
  }
}
//======================================================================

void hamvec_func3( int nspin, int nTerm, std::complex<double> *coeff_lst_zplx, size_t *nbody_lst, size_t *pos_i_idx, size_t *pos_i_lst, size_t *dim_i_lst, size_t *mat_i_idx, std::complex<double> *mat_i_lst, size_t vlen, size_t *nspin_dim, size_t *nspin_m_lst, size_t *nspin_n_lst, std::complex<double> *v, std::complex<double> *w, std::complex<double> *w_med )
{
/*
  Calculate the action of a Hamiltonian on a state vector.
  The Hamiltonian can be decomposed into many terms, where each term 
  consists of 
  
  input:
    
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
    mat_i_lst,      list of the operators of bodies in each interaction; 
    vlen,           dimension of the state vector;
    nspin_dim,      list of the dimension of each body;
    nspin_m_lst,    list of the dimension of first m bodies;
    nspin_n_lst,    list of the dimension of last n bodies;
    v,              input state vector;
    w_med,          intermediate state vector;
  
  output:
    
    w,              output state vector, after the Hamiltonian acting on
                    the input state;
*/
  size_t                nT, nbody, nb;
  size_t                idx, pos_i, dim_i;
  std::complex<double>  coeff;
  
  size_t                i, m, n;
  
  std::complex<double> zero(0.0,0.0);
  
  // initialization of w;
//  #pragma omp parallel
//  {
//  #pragma omp parallel for firstprivate( zero, w )
  for ( i = 0; i < vlen; i++ ) w[i] = zero;
//  }
  
  // kron_func_v4;
  int bds = 512, bdy;
  
  for ( nT = 0; nT < nTerm; nT++ )
  {
    coeff = coeff_lst_zplx[ nT ];
    nbody = nbody_lst[ nT ];
    
    cblas_zcopy( vlen, v, 1, w_med, 1 );
    for ( nb = 0; nb < nbody; nb++ )
    {
      idx = pos_i_idx[ nT ] + nbody - 1 - nb;
      pos_i = pos_i_lst[ idx ];
      dim_i = dim_i_lst[ idx ];
      
      m = nspin_m_lst[ pos_i ];
      n = nspin_n_lst[ pos_i ];
      
//      kron_func_v2( m, dim_i, n, &mat_i_lst[ mat_i_idx[ idx ] ], w_med );
      kron_func_v3( m, dim_i, n, &mat_i_lst[ mat_i_idx[ idx ] ], w_med );
//      kron_func_v32( m, dim_i, n, &mat_i_lst[ mat_i_idx[ idx ] ], w_med );
      
      // kron_func_v4;
      bdy = bds / dim_i;
//      bdy = m * n;
      while ( ( ( m * n ) % bdy ) != 0 ) bdy--;
      bds = dim_i * bdy;
//      kron_func_v4( m, dim_i, n, &mat_i_lst[ mat_i_idx[ idx ] ], w_med, bdy );
//      kron_func_v5( m, dim_i, n, &mat_i_lst[ mat_i_idx[ idx ] ], w_med, bdy );
//      kron_func_v6( m, dim_i, n, &mat_i_lst[ mat_i_idx[ idx ] ], w_med, bdy );
//      kron_func_v7( m, dim_i, n, &mat_i_lst[ mat_i_idx[ idx ] ], w_med, bdy );
    }
    cblas_zaxpy( vlen, &coeff, w_med, 1, w, 1 );
  }
}

// interface for FORTRAN;

#define HAMVEC_FUNC3       hamvec_func3_

#if defined(__cplusplus)
extern "C" {
#endif /* __cplusplus */

void HAMVEC_FUNC3( int *nspin_ptr, int *nTerm_ptr, std::complex<double> *coeff_lst_zplx_ptr, size_t *nbody_lst_ptr, size_t *pos_i_idx_ptr, size_t *pos_i_lst_ptr, size_t *dim_i_lst_ptr, size_t *mat_i_idx_ptr, std::complex<double> *mat_i_lst_ptr, size_t *ham_dim_ptr, size_t *nspin_dim_ptr, size_t *nspin_m_lst_ptr, size_t *nspin_n_lst_ptr, std::complex<double> *v_ptr, std::complex<double> *w_ptr, std::complex<double> *w_med_ptr );

#if defined(__cplusplus)
}
#endif /* __cplusplus */

void HAMVEC_FUNC3( int *nspin_ptr, int *nTerm_ptr, std::complex<double> *coeff_lst_zplx_ptr, size_t *nbody_lst_ptr, size_t *pos_i_idx_ptr, size_t *pos_i_lst_ptr, size_t *dim_i_lst_ptr, size_t *mat_i_idx_ptr, std::complex<double> *mat_i_lst_ptr, size_t *ham_dim_ptr, size_t *nspin_dim_ptr, size_t *nspin_m_lst_ptr, size_t *nspin_n_lst_ptr, std::complex<double> *v_ptr, std::complex<double> *w_ptr, std::complex<double> *w_med_ptr )
{
  int                   nspin           = *nspin_ptr;
  int                   nTerm           = *nTerm_ptr;
  std::complex<double>  *coeff_lst_zplx = coeff_lst_zplx_ptr;
  size_t                *nbody_lst      = nbody_lst_ptr;
  size_t                *pos_i_idx      = pos_i_idx_ptr;
  size_t                *pos_i_lst      = pos_i_lst_ptr;
  size_t                *dim_i_lst      = dim_i_lst_ptr;
  size_t                *mat_i_idx      = mat_i_idx_ptr;
  std::complex<double>  *mat_i_lst      = mat_i_lst_ptr;
  size_t                vlen            = *ham_dim_ptr;
  size_t                *nspin_dim      = nspin_dim_ptr;
  size_t                *nspin_m_lst    = nspin_m_lst_ptr;
  size_t                *nspin_n_lst    = nspin_n_lst_ptr;
  std::complex<double>  *v              = v_ptr;
  std::complex<double>  *w              = w_ptr;
  std::complex<double>  *w_med          = w_med_ptr;
  
  hamvec_func3( nspin, nTerm, coeff_lst_zplx, nbody_lst, pos_i_idx, pos_i_lst, dim_i_lst, mat_i_idx, mat_i_lst, vlen, nspin_dim, nspin_m_lst, nspin_n_lst, v, w, w_med );
}

