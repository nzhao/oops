//#ifdef __cplusplus
//extern"C" {
//#endif

// interface to zgexpv in Expokit, dense a;
int krylov_zgexpv(const int n, const double _Complex *a, const double _Complex *v, const int tn, const double *ta, double _Complex *w_seq, const int klim, const int m, const double tol,  const int itrace);

// interface to zgexpv in Expokit, sparse coo a;
int krylov_zcooexpv(const int n, const int nz, const int *ia, const int *ja, const double _Complex *a, const double _Complex *v, const int tn, const double *ta, double _Complex *w_seq, const int klim, const int m, const double tol,  const int itrace);

//#ifdef __cplusplus
//}
//#endif

