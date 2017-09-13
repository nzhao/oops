#ifdef __cplusplus
extern"C" {
#endif
// exp(-i * h * t);
void mpi_hzcsrkexpv_(long *n, long *m, double _Complex *v, double _Complex *w, \
  double *tol, double *anorm, double _Complex *wsp, long *lwsp, \
  long *iwsp, long *liwsp, long *itrace, long *iflag, long *nspin, \
  ham *hlst, double _Complex *prefactor, long *nterm, double *tlst, long *tn, \
  ham *alst, long *n_op, double _Complex *oplst, \
  MPI_Request *reqs, MPI_Status *status, double _Complex *cache, long *nchunk);
#ifdef __cplusplus
}
#endif
