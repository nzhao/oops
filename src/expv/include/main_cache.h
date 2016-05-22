#ifdef __cplusplus
extern"C" {
#endif

// Expokit in FORTRAN using GPU;
void main_cache_( size_t *nspin, size_t *nTerm, double *coeff_lst, size_t *nbody_lst, size_t *pos_i_idx, size_t *pos_i_lst, size_t *dim_i_lst, size_t *mat_i_idx, std::complex<double> *mat_i_lst, size_t *ham_dim, size_t *nspin_dim, std::complex<double> *v, size_t *pos_idx, size_t *mat_idx, size_t *k, size_t *maxThreadsPerBlock, size_t *maxGridSize, size_t *tn, double *ta, size_t *m, double *tol, size_t *itrace, std::complex<double> *w_seq, size_t *w_seq_len );

#ifdef __cplusplus
}
#endif

