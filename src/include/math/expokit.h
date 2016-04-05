#ifdef __cplusplus
extern"C" {
#endif

// irreducible rational Pade approximation in Expokit;
void zgpadm_(int *ideg, int *m, double *t, std::complex<double> *H, int *ldh, std::complex<double> *wsp, int *lwsp, int *ipiv, int *iexph, int *ns, int *iflag);

// zgexpv in Expokit, implement blas on matvec;
void krylov_zgexpv_(const int *n, const double _Complex *a, const double _Complex *v, const int *tn, const double *ta, const int *m, const double *tol, double *anorm, double _Complex *wsp, const int *lwsp, int *iwsp, const int *liwsp, const int *itrace, int *iflag, double _Complex *w_seq, const int *w_seq_len);

// zgexpv in Expokit, implement blas on matvec;
void krylov_zcooexpv_(const int *n, const int *nz, const int *ia, const int *ja, const double _Complex *a, const double _Complex *v, const int *tn, const double *ta, const int *m, const double *tol, double *anorm, double _Complex *wsp, const int *lwsp, int *iwsp, const int *liwsp, const int *itrace, int *iflag, double _Complex *w_seq, const int *w_seq_len);

#ifdef __cplusplus
}
#endif

