#ifdef __cplusplus
extern"C" {
#endif

#ifndef HAM_H
#define HAM_H

typedef struct {
  unsigned int      *nzia;
  unsigned int      *nzja;
  double            *nza;
  
  unsigned int      *nzmzia;
  unsigned int      *nzmzja;
  double            *nzmza;
  
  unsigned int      *nxyia;
  unsigned int      *nxyja;
  double            *nxyax;
  double            *nxyay;
  
  unsigned int      *nxymzia;
  unsigned int      *nxymzja;
  double            *nxymzax;
  double            *nxymzay;
  
  unsigned int      *nxymxyia;
  unsigned int      *nxymxyja;
  double            *nxymxyaxx;
  double            *nxymxyaxy;
  double            *nxymxyayx;
  double            *nxymxyayy;
} ham;

#endif// HAM_H

//  get the index of c_{n beta} in coeff_lst[:];
unsigned int getSingleIdx(const unsigned int nspin, const unsigned int n, \
  const unsigned int beta);

//  get the index of c_{m alpha n beta} in coeff_lst[:];
unsigned int getPairIdx(const unsigned int nspin, const unsigned int n, \
  const unsigned int beta, const unsigned int m, const unsigned int alpha);

//  convert coeff_lst to Ham in CSR format;
ham hamGen(int nspin, double *coeff_lst);
void hamFree(ham h);

#ifdef __cplusplus
}
#endif
