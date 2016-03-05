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
MatExpVector::MatExpVector(const SumKronProd& skp, cx_double prefactor, const vec& time_list)
{
    _skp = skp;
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
    size_t            dim = _skp.getDim();
    DIM_LIST   spinDim = _skp.getDimList();  //spin_dim
    size_t      factor_num = _skp.getKronNum();  //nspin
    size_t       prod_size = _skp.getKronProdSize();  //nTerm

    vector<MULTIPLIER> coeff = _skp.getCoeffList();  //coeff_list
    vector<size_t>        nBody = _skp.getMatNumList(); //nBody_list
    vector<INDICES>  indices_list = _skp.getIndicesList(); 
    vector<TERM>        term_list = _skp.getTermList();

    INDICES full_indices = join_all(indices_list); //pos_list
    TERM       full_term = join_all(term_list);    //matC

    vector<size_t> dim_vector; //dim_list
    for(size_t i=0; i<full_term.size(); ++i)
        dim_vector.push_back( full_term[i].n_cols );

    //////////////////////////////////////////////////////////////////////////////
    //parameter preparation
    size_t nSpin = factor_num;
    size_t nTerm = prod_size;
    double * coeff_list = coeff.data();
    size_t * nBody_list = nBody.data();
    //double * coeff_list = new double [nTerm];
    //size_t * nBody_list = new size_t [nTerm];
    size_t * pos_offset = new size_t [nTerm+1];
    //for(int i=0; i<nTerm;++i)
    //{
        //coeff_list[i] = coeff[i];
        //nBody_list[i] = nBody[i];
    //}
    pos_offset[0]=0; partial_sum (nBody_list, nBody_list+nTerm, pos_offset+1);
    size_t total_nbody = pos_offset[nTerm];

    size_t * pos_list   = new size_t [total_nbody];
    size_t * dim_list   = new size_t [total_nbody];
    size_t * mat_offset = new size_t [total_nbody+1];
    size_t * dim2 = new size_t [total_nbody];
    for(int i=0; i<total_nbody; ++i)
    {
        pos_list[i] = full_indices[i];
        dim_list[i] = dim_vector[i];
        dim2[i] = dim_list[i]*dim_list[i];
    }
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
    size_t nDim = dim;
    size_t * spin_dim   = new size_t [nSpin];
    for(int i=0; i<nSpin; ++i)
        spin_dim[i] = spinDim[i];

    complex<double> * vecC = new complex<double> [ nDim ];
    for (int i=0; i<nDim; ++i )
        vecC[i] = 0.0;
    vecC[0] = 1.0;
    
    size_t nt = _time_list.n_elem; 
    double * tlist = new double [nt];
    for(int i=0; i<nt; ++i)
        tlist[i] = _time_list(i);

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
