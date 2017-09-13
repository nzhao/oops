#include "ham.h"

#ifdef __cplusplus
extern"C" {
#endif

void mpi_zcsrkhmv_(const unsigned int nspin, const ham h, \
  const double _Complex coeff, const double _Complex *v, double _Complex *u, \
  const unsigned int nchunk, MPI_Request *reqs, \
  MPI_Status *status, double _Complex *cache);

#ifdef __cplusplus
}
#endif
