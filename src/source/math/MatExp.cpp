#include "include/math/MatExp.h" 
#include "include/math/expokit.h" 
#include <complex>

////////////////////////////////////////////////////////////////////////////////
//{{{  MatExp
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
//}}}
////////////////////////////////////////////////////////////////////////////////




////////////////////////////////////////////////////////////////////////////////
//{{{  MatExpVector
MatExpVector::MatExpVector(const SumKronProd& skp, cx_double prefactor, const cx_vec& v, const vec& time_list)
{
    _skp = skp;
    _vector = v;
    _prefactor = prefactor;
    _time_list = time_list;

    _nTime = time_list.n_elem;
    _dim = skp.getDim();

    _klim = 1;
    _krylov_m = 30;
    _krylov_tol = 1e-7;
    _itrace = 1;
}

void MatExpVector::run()
{
    //////////////////////////////////////////////////////////////////////////////
    //parameter preparation
    DIM_LIST           spinDim = _skp.getDimList();
    vector<MULTIPLIER> coeff = _skp.getCoeffList();
    vector<size_t>     nBody = _skp.getMatNumList();
    vector<INDICES>    indices_list = _skp.getIndicesList(); 
    vector<TERM>       term_list = _skp.getTermList();
    INDICES            full_indices = join_all(indices_list);
    TERM               full_term = join_all(term_list);

    vector<size_t> dim_vector;
    for(size_t i=0; i<full_term.size(); ++i)
        dim_vector.push_back( full_term[i].n_cols );

    size_t nSpin = _skp.getKronNum();
    size_t nTerm = _skp.getKronProdSize();
    double * coeff_list = coeff.data();
    size_t * nBody_list = nBody.data();
    
    size_t * pos_offset = new size_t [nTerm+1];
    pos_offset[0]=0; partial_sum (nBody_list, nBody_list+nTerm, pos_offset+1);
    size_t total_nbody = pos_offset[nTerm];

    size_t * pos_list   = full_indices.data();
    size_t * dim_list   =  dim_vector.data();
    size_t * dim2 = new size_t [total_nbody];
    for(int i=0; i<total_nbody; ++i)
        dim2[i] = dim_list[i]*dim_list[i];
    
    size_t * mat_offset = new size_t [total_nbody+1];
    mat_offset[0]=0; partial_sum(dim2, dim2+total_nbody, mat_offset+1);
    size_t total_dim = mat_offset[total_nbody];

    complex<double> * matC = new complex<double> [total_dim];
    int count = 0;
    for(int i=0; i<full_term.size();++i)
        for(int j=0; j<dim2[i]; ++j)
        {
            matC[count] = full_term[i](j);
            count ++;
        }
    size_t nDim = _skp.getDim();
    size_t * spin_dim   = spinDim.data();

    complex<double> * vecC = _vector.memptr();
    
    size_t nt = _time_list.n_elem; 
    double * tlist =  _time_list.memptr();

    complex<double> * w_seq;
    size_t w_seq_len= nDim * nt;
    w_seq = new complex<double> [w_seq_len];
    //////////////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////////////
    // Print to screen
    {/*{{{*/
    cout << "#1. nSpin = " << nSpin << endl << endl;
    cout << "#2. nTerm = " << nTerm << endl << endl;
    cout << "#3. coeff_list = " << endl;
    for(int i=0; i<nTerm;++i)
        cout << coeff_list[i] << "\t";
    cout << endl << endl;
    cout << "#4. nBody_list = " << endl;
    for(int i=0; i<nTerm;++i)
        cout << nBody_list[i] << "\t";
    cout << endl << endl;
    cout << "#5. pos_offset = " << endl;
    for(int i=0; i<nTerm+1;++i)
        cout <<pos_offset[i] << "\t";
    cout << endl << endl;

    cout << "#6. pos_list = " << endl;
    for(int i=0; i<total_nbody; ++i)
        cout << pos_list[i] << "\t";
    cout << endl << endl;
    cout << "#7. dim_list = " << endl;
    for(int i=0; i<total_nbody; ++i)
        cout <<dim_list[i] << "\t";
    cout << endl << endl;
    cout << "#8. mat_offset = " << endl;
    for(int i=0; i<total_nbody+1; ++i)
        cout <<mat_offset[i] << "\t";
    cout << endl << endl;

    cout << "#9. matC = " << endl;
    for(int i=0; i<total_dim; ++i)
        cout << matC[i] << "\t";
    cout << endl << endl;
    cout << "#10. nDim = " << nDim << endl << endl;
    cout << "#11. sin_dim = " << endl;
    for(int i=0; i<nSpin; ++i)
        cout << spin_dim[i] << "\t";
    cout << endl << endl;
    cout << "#12. vecC = " << endl;
    for(int i=0; i<nDim; ++i) 
        cout << vecC[i] << "\t";
    cout << endl << endl;
    cout << "#13. total_nbody = " << total_nbody << endl << endl;
    cout << "#14. total_dim = " << total_dim << endl << endl;
    cout << "#15. klim = " << _klim << endl << endl;
    cout << "#16. nt = " << nt << endl << endl;
    cout << "#17. tlist = " << endl;
    for(int i=0; i<nt; ++i)
        cout << tlist[i] << "\t";
    cout << endl << endl;
    cout << "#18. m= " << _krylov_m << endl << endl;
    cout << "#19. tol = " << _krylov_tol << endl << endl;
    cout << "#20. itrace = " << _itrace << endl << endl;
    cout << "#21. w_seq = " << w_seq << endl << endl;
    cout << "#22. w_seq_len = " << w_seq_len << endl << endl;

    }/*}}}*/
    
    //main_mkl_(  &nSpin,
                //&nTerm, 
                //coeff_list, 
                //nBody_list, 
                //pos_offset, 
                //pos_list, 
                //dim_list, 
                //mat_offset, 
                //matC, 
                //&nDim, 
                //spin_dim, 
                //vecC, 
                //&total_nbody, 
                //&total_dim, 
                //&klim, 
                //&nt, 
                //tlist,
                //&m, 
                //&tol, 
                //&itrace,
                //w_seq,
                //&w_seq_len );
}
//}}}
////////////////////////////////////////////////////////////////////////////////
