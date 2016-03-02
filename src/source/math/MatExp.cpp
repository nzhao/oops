#include "include/math/MatExp.h" 
#include "include/math/expokit.h" 
#include <complex>

MatExp::MatExp(const cx_mat& m, cx_double prefactor, MatExpMethod method)
{
    _matrix = m;
    _prefactor = prefactor;
    _method = method;
}

void MatExp::run()
{
    switch (_method) {
        case ArmadilloExpMat:
            _resMatrix = armadillo_exp_mat();
            break;
        case PadeApproximation:
            _resMatrix = pade_exp_mat();
            break;
        default:
            cout << "Exp method not sopport." << endl;
            assert(0);
    }
}

cx_mat MatExp::pade_exp_mat()
{
    cx_mat res;
    
    // computes exp(t*H), irreducible rational Pade approximation;
    int       ideg(6);// 
    int       m(0);  // order of H;
    double    t(1.0); // time-scale;
    std::complex<double> *H(NULL);// matrix;
    int       ldh(0);
    std::complex<double> *wsp(NULL);
    int       lwsp(0);
    int       *ipiv(NULL);
    int       iexph(0);
    int       ns(0);
    int       iflag(0);
    
    res = _prefactor * _matrix;
    m = res.n_cols;
    H = res.memptr();
    if (H == NULL)
      return zeros<cx_mat>(m, m);
    ldh = m;
    lwsp = 4 * m * m + ideg + 1;
    
    wsp = new (std::nothrow) std::complex<double> [lwsp];
    ipiv = new (std::nothrow) int [m];
    if ((ipiv == NULL) || (wsp == NULL))
    {
      delete []wsp;
      delete []ipiv;
      return zeros<cx_mat>(m, m);
    }
    
    zgpadm_(&ideg, &m, &t, H, &ldh, wsp, &lwsp, ipiv, &iexph, &ns, &iflag);
    
    if (iflag < 0)
    {
      std::cout << "problem in ZGPADM, iflag = " << iflag << std::endl;
      delete []wsp;
      delete []ipiv;
      return zeros<cx_mat>(m, m);
    }
    res = cx_mat( &wsp[iexph-1], m, m);// zero-based numbering;
    
    delete []wsp;
    delete []ipiv;
    return res;
}

