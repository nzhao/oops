#include "include/math/MatExp.h" 
#include "include/math/expokit.h" 
#include "include/math/main_mkl.h" 
#include "include/math/main_cache.h" 
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
MatExpVector::MatExpVector(const SumKronProd& skp, const cx_vec& v, const vec& time_list, MatExpVectorMethod method)
{
    _method = method;
    _skp = skp;
    _vector = v;
    //_prefactor = prefactor;
    _prefactor = cx_double(0.0, -1.0);
    _time_list = time_list;

    _nTime = time_list.n_elem;
    _dim = skp.getDim();

    _klim = 10;//  Lanczos factorization length;
    _krylov_m = 30;// _krylov_m = 30, optimized in Expokit;
    _krylov_tol = 1e-12;
    _itrace = 0;
    _is_print_skp = false;
}
MatExpVector::MatExpVector(const cx_mat& m, const cx_vec& v, const cx_double prefactor, const vec& time_list)
{
    _method = Explicit;
    _prefactor = prefactor;
    _matrix = _prefactor*m;
    _vector = v;
    _time_list = time_list;

    _nTime = time_list.n_elem;
    _dim = m.n_cols;

    _klim = 10;//  Lanczos factorization length;
    _krylov_m = 30;// _krylov_m = 30, optimized in Expokit;
    _krylov_tol = 1e-12;
    _itrace = 0;
    _is_print_skp = false;
}
MatExpVector::MatExpVector(const sp_cx_mat& m, const cx_vec& v, const cx_double prefactor, const vec& time_list)
{
    _method = ExplicitSparse;
    _prefactor = prefactor;
    _sp_matrix = _prefactor*m;
    _vector = v;
    _time_list = time_list;

    _nTime = time_list.n_elem;
    _dim = m.n_cols;

    _klim = 10;//  Lanczos factorization length;
    _krylov_m = 30;// _krylov_m = 30, optimized in Expokit;
    _krylov_tol = 1e-12;
    _itrace = 0;
    _is_print_skp = false;
}

cx_mat MatExpVector::run()
{
    cx_mat res;
    switch (_method) {
        case Explicit:
            res = runExplicit();
            break;
        case ExplicitSparse:
            res = runExplicitSparse();
            break;
        case Inexplicit:
            res = runInexplicit();
            break;
        case InexplicitGPU:
            res = runInexplicitGPU();
            break;
        default:
            cout << "Exp method not sopport." << endl;
            assert(0);
    }
    return res;
}
    
cx_mat MatExpVector::runExplicit()
{/*{{{*/
    std::complex<double> * mat = _matrix.memptr();
    std::complex<double> * vecC = _vector.memptr();

    // krylov_zgexpv and krylov_zcooexpv;
    std::complex<double> *w_seq = new std::complex<double> [_nTime * _dim];
    
    int err = krylov_zgexpv(_dim, (double _Complex *)mat, (double _Complex *)vecC, _nTime, &_time_list[0], (double _Complex *)w_seq, _klim, _krylov_m, _krylov_tol,  _itrace);

    cx_mat res(&w_seq[0], _dim, _nTime);
    _resVectorList = res;

    delete[] w_seq;
    return _resVectorList;
}/*}}}*/

cx_mat MatExpVector::runExplicitSparse()
{/*{{{*/
    sp_cx_mat::const_iterator start = _sp_matrix.begin();
    sp_cx_mat::const_iterator end   = _sp_matrix.end();

    std::complex<double> * vecC = _vector.memptr();

    int nz = distance(start , end);

    int * ia = new int [nz];
    int * ja = new int [nz];
    std::complex<double> * a = new std::complex<double> [nz];
    //int i=0;
    for(sp_cx_mat::const_iterator it = start; it != end; ++it)
    {
        int i = distance(start, it);
        ia[i] = it.row()+1;
        ja[i] = it.col()+1;
        a[i] = (*it);
        //i++;
    }

    // krylov_zgexpv and krylov_zcooexpv;
    std::complex<double> *w_seq = new std::complex<double> [_nTime * _dim];
    
    int err = krylov_zcooexpv(_dim, nz, ia, ja, (double _Complex *)a, (double _Complex *)vecC, _nTime, &_time_list[0], (double _Complex *)w_seq, _klim, _krylov_m, _krylov_tol, _itrace);

    cx_mat res(&w_seq[0], _dim, _nTime);
    _resVectorList = res;

    delete[] ia;
    delete[] ja;
    delete[] a;
    delete[] w_seq;
    return _resVectorList;
}/*}}}*/

cx_mat MatExpVector::runInexplicit()
{/*{{{*/
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
    if(_is_print_skp)
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
    
    main_mkl_(  &nSpin,
                &nTerm, 
                coeff_list, 
                nBody_list, 
                pos_offset, 
                pos_list, 
                dim_list, 
                mat_offset, 
                matC, 
                &nDim, 
                spin_dim, 
                vecC, 
                &total_nbody, 
                &total_dim, 
                &_klim, 
                &nt, 
                tlist,
                &_krylov_m, 
                &_krylov_tol, 
                &_itrace,
                w_seq,
                &w_seq_len );
    cx_mat res(w_seq, nDim, nt);
    _resVectorList = res;

    delete[] pos_offset;
    delete[] dim2;
    delete[] mat_offset;
    delete[] matC;
    delete[] w_seq;

    return _resVectorList;
}/*}}}*/

cx_mat MatExpVector::runInexplicitGPU()
{/*{{{*/
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
    if(_is_print_skp)
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
    
    size_t  _maxThreadsPerBlock  = 256;
    size_t  _maxGridSize[3]      = {2147483647,65535,65535};
    
    main_cache_(  &nSpin,
                &nTerm, 
                coeff_list, 
                nBody_list, 
                pos_offset, 
                pos_list, 
                dim_list, 
                mat_offset, 
                matC, 
                &nDim, 
                spin_dim, 
                vecC, 
                &total_nbody, 
                &total_dim, 
                &_klim, 
                &_maxThreadsPerBlock,
                _maxGridSize,
                &nt, 
                tlist,
                &_krylov_m, 
                &_krylov_tol, 
                &_itrace,
                w_seq,
                &w_seq_len );
    cx_mat res(w_seq, nDim, nt);
    _resVectorList = res;

    delete[] pos_offset;
    delete[] dim2;
    delete[] mat_offset;
    delete[] matC;
    delete[] w_seq;

    return _resVectorList;
}/*}}}*/
//}}}
////////////////////////////////////////////////////////////////////////////////

