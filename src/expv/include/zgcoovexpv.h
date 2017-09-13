#ifdef __cplusplus
extern"C" {
#endif
// exp(-i * h * t);
void zgcoovexpv_(long *n, long *m, double *t,double _Complex *v, \
  double _Complex *w, double *tol, double *anorm, \
  double _Complex *wsp, long *lwsp, long *iwsp, \
  long *liwsp, long *itrace, long *iflag);
#ifdef __cplusplus
}
#endif
