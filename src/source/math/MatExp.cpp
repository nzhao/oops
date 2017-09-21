#include "include/math/MatExp.h" 
#include "expv/include/expv.h" 
#include <complex.h>
#undef complex

////////////////////////////////////////////////////////////////////////////////
//{{{  MatExp
MatExp::MatExp(const cx_mat& m, cx_double prefactor, MatExpMethod method)
{
    _matrix = m;
    _prefactor = prefactor;
    _method = method;
}

MatExp::MatExp(const int nspin, const vec coeff_list, const cx_vec local_phi0,\
		       const vector<uvec> op_list, const cx_double PREFACTOR,\ 
		       const vec TIME_LIST, MatExpMethod method)
{
    _nspin = nspin;
    _coeff_list = coeff_list;
    _v = local_phi0;
    _op_list = op_list;
    _time_list = TIME_LIST;
    _prefactor = PREFACTOR;
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
	    case LargeMatExpv:
	        _resMatrix = large_mat_expv(); 
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

cx_mat MatExp::large_mat_expv()
{
	long            n;
    long            m;
    double _Complex *v;
    double _Complex *w;
    double          tol;
    double          anorm;
    double _Complex *wsp;
    long            lwsp;
    long            *iwsp;
    long            liwsp;
    long            itrace;
    long            iflag;
    long            lnspin;
    double          *coeff_lst;
    double          *op_coeff_lst;
    long            nterm;
    double          *tlst;
    long            tn;
    double _Complex *oplst;
    long            n_op;
    MPI_Request     *reqs;
    MPI_Status      *status;
    long            nchunk;
    long			chunkmax;
    double _Complex *cache;
    ham             h, **hlst, *alst;
    double _Complex prefactor;
    int             rank,nprocess;
    
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&nprocess);
    
    lnspin = _nspin;
    n = _v.size();
    m = 30;// _krylov_m;
    tol = 1.0e-14;// _krylov_tol;
    itrace = 0;// for debug;
    iflag = 0;
    prefactor = 1.0 * I;

    tn = _time_list.size();
    vec temp = _time_list;
    tlst = temp.memptr();
    
    lwsp = n * (m + 2) + 5 * (m + 2) * (m + 2) + 7;
    liwsp = m + 2;
    nterm = 3 * lnspin + 9 * lnspin * (lnspin - 1) / 2;
    chunkmax = 1048576;
    nchunk = (n - 1) / chunkmax + 1;
    n_op = _op_list.size();
    
    v = (double _Complex*)malloc(sizeof(double _Complex) * n);
    w = (double _Complex*)malloc(sizeof(double _Complex) * n);
    cache = (double _Complex*)malloc(sizeof(double _Complex) * n * 2);
    wsp = (double _Complex*)malloc(sizeof(double _Complex) * lwsp);
    iwsp = (long*)malloc(sizeof(long) * liwsp);
    coeff_lst = (double*)malloc(sizeof(double) * nterm);
    reqs = (MPI_Request*)malloc(sizeof(MPI_Request) * nchunk * 2);
    status = (MPI_Status*)malloc(sizeof(MPI_Status) * nchunk * 2);
    op_coeff_lst = (double*)malloc(sizeof(double) * nterm * n_op);
    oplst = (double _Complex*)malloc(sizeof(double _Complex) * tn * n_op);
    
    //coeff transform(attention on difference between sigma)
    int idx = 0;
    for (int i = 0; i < lnspin - 1; i++){
      for (int j = i + 1; j < lnspin; j++){
        for (int a = 0; a < 3; a++){
          for (int b = 0; b < 3; b++){
            coeff_lst[getPairIdx(lnspin, i, a, j, b)] = _coeff_list[idx] / 4;
            op_coeff_lst[getPairIdx(lnspin, i, a, j, b)] = 0.0;
            idx++;
          }
        }
      }
    }
    for (int i = 0; i < lnspin; i++){
      for (int a = 0; a < 3; a++){
       coeff_lst[getSingleIdx(lnspin, i, a)] = _coeff_list[idx] / 2;
       op_coeff_lst[getSingleIdx(lnspin, i, a)] = 0.0;
       idx++;
      }
    }
    
    anorm = 0.0;
    for (int i = 0; i < nterm; i++)
      anorm += fabs(coeff_lst[i]);
 
    printf("[Data] anorm = %e, typical_time = %e\n",anorm,1.0 / anorm);
    
    #pragma omp parallel for
    for (long i = 0; i < n; i++){
      v[i] = 0.0;
      w[i] = 0.0;
      cache[i] = 0.0;
      cache[i + n] = 0.0;
      for (long j = 0; j < m + 2; j++)
        wsp[n * j + i] = 0.0;
    }
    
    #pragma omp parallel for
    for (int i = 0; i < n; i++)
      v[i] = real(_v[i]) + I*imag(_v[i]);
    
    //operator:E [*] A
    for (int i = 0; i < n_op; i++){
        if( _op_list[i][0] == 1)
            op_coeff_lst[ i * nterm + getSingleIdx(lnspin,_op_list[i][1],_op_list[i][2])] = 1.0 / 2;
        else if( _op_list[i][0] == 2)
            op_coeff_lst[ i * nterm + getPairIdx(lnspin,_op_list[i][1],_op_list[i][2],_op_list[i][3],_op_list[i][4])] = 1.0 / 4;
        else{
            printf("wrong input!!!");
            exit(0);}
    } 
    
    
    h = hamGen(lnspin, coeff_lst);
    
    double  *op_temp;
    op_temp = (double*)malloc(sizeof(double) * nterm);
    alst = (ham*)malloc(sizeof(ham) * n_op);
    for (int i = 0; i < n_op; i++){
        for (int j = 0; j < nterm; j++){
            op_temp[j] = op_coeff_lst[ j + (i * nterm) ];}
        alst[i] = hamGen(lnspin, op_temp);
    }
    
    double  tic, tac;
    tic = omp_get_wtime();   
	
	bool space = false;//later
    printf("[expv] n = %d\n",n);
    printf("[expv] nspin = %d\n",lnspin);
    printf("[expv] lwsp = %d\n",lwsp);
    printf("[expv] nterm = %d\n",nterm);
    printf("[expv] n = %p\n",&n);
    printf("[expv] nspin = %p\n",&lnspin);
    printf("[expv] lwsp = %p\n",&lwsp);
    printf("[expv] nterm = %p\n",&nterm);
    mpi_hzcsrkexpv_(&n, &m, v, w, &tol, &anorm, wsp, &lwsp, iwsp, &liwsp, &itrace,\
                  &iflag, &lnspin, &h, &prefactor, &nterm,tlst, &tn, alst, &n_op, \
                  oplst, reqs, status, cache, &nchunk);
                
    tac = omp_get_wtime();
    
    printf("[debug] rank %4d: mpi_Hzcsrkexpv works %f second!\n",rank, tac-tic);// for debug;
    
    cx_mat _resMatrix = zeros<cx_mat>(tn,n_op);
    for (int j =0; j < n_op; j++){
        cx_vec res = zeros<cx_vec>(tn,1);
        for (int i = 0; i < tn; i++){
            res[i] = cx_double(creal(oplst[i+j*tn]),cimag(oplst[i+j*tn]));
        }
        _resMatrix.col(j) = res;
    }


    for (int i = 0; i < n_op; i++) hamFree(alst[i]);
    hamFree(h);
    
    free(op_temp);
    free(op_coeff_lst);
    free(oplst);
    free(status);
    free(reqs);
    free(coeff_lst);
    free(iwsp);
    free(wsp);
    free(cache);
    free(w);
    free(v);
    
    return _resMatrix;
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

MatExpVector::MatExpVector(const int nspin, const vector<vec> coeff_list, const cx_vec local_phi0,const vector<uvec> op_list,\
                           const cx_vec final_vec, const cx_double PREFACTOR, const vec TIME_LIST, MatExpVectorMethod method)
{
    _nspin = nspin;//Liouville space
    _coeff_list = coeff_list;
    _v = local_phi0;
    _op_list = op_list;
    _f = final_vec;
    _time_list = TIME_LIST;
    _prefactor = PREFACTOR;
    _method = method;
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
        case LargeVecExpv:
	        res = runLargeVecExpv();
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

cx_mat MatExpVector::runLargeVecExpv()
{/*{{{*/
	long            n;
    long            m;
    double _Complex *v;
    double _Complex *f;
    double _Complex *w;
    double          tol;
    double          anorm;
    double _Complex *wsp;
    long            lwsp;
    long            *iwsp;
    long            liwsp;
    long            itrace;
    long            iflag;
    long            lnspin;
    double          *coeff_lst1, *coeff_lst2;
    double          *op_coeff_lst;
    long            nterm;
    double          *tlst;
    long            tn;
    double _Complex *oplst;
    long            n_op;
    MPI_Request     *reqs;
    MPI_Status      *status;
    long            nchunk;
    long			chunkmax;
    double _Complex *cache;
    ham             h1, h2, *alst;
    double          h3;
    double _Complex prefactor;
    int             rank,nprocess;
    
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&nprocess);

    lnspin = _nspin;
    n = _v.size();
    m = 30;// _krylov_m;
    tol = 1.0e-14;// _krylov_tol;
    itrace = 0;// for debug;
    iflag = 0;
    prefactor = real(_prefactor) + imag(_prefactor) * I;
    
    tn = _time_list.size();
    vec temp = _time_list;
    tlst = temp.memptr();
    
    lwsp = n * (m + 2) + 5 * (m + 2) * (m + 2) + 7;
    liwsp = m + 2;
    nterm = 3 * lnspin + 9 * lnspin * (lnspin - 1) / 2;
    chunkmax = 1048576;
    nchunk = (n - 1) / chunkmax + 1;
    n_op = _op_list.size();
    
    v = (double _Complex*)malloc(sizeof(double _Complex) * n);
    f = (double _Complex*)malloc(sizeof(double _Complex) * n);
    w = (double _Complex*)malloc(sizeof(double _Complex) * n);
    cache = (double _Complex*)malloc(sizeof(double _Complex) * n * 2);
    wsp = (double _Complex*)malloc(sizeof(double _Complex) * lwsp);
    iwsp = (long*)malloc(sizeof(long) * liwsp);
    coeff_lst1 = (double*)malloc(sizeof(double) * nterm);
    coeff_lst2 = (double*)malloc(sizeof(double) * nterm);
    reqs = (MPI_Request*)malloc(sizeof(MPI_Request) * nchunk * 2);
    status = (MPI_Status*)malloc(sizeof(MPI_Status) * nchunk * 2);
    op_coeff_lst = (double*)malloc(sizeof(double) * nterm * n_op);
    oplst = (double _Complex*)malloc(sizeof(double _Complex) * tn * n_op);
    
    int idx = 0;
    for (int i = 0; i < lnspin - 1; i++){
      for (int j = i + 1; j < lnspin; j++){
        for (int a = 0; a < 3; a++){
          for (int b = 0; b < 3; b++){
            coeff_lst1[getPairIdx(lnspin, i, a, j, b)] =  _coeff_list[0][idx] / 4;
            coeff_lst2[getPairIdx(lnspin, i, a, j, b)] =  _coeff_list[1][idx] / 4;
            idx++;
          }
        }
      }
    }
    for (int i = 0; i < lnspin; i++){
      for (int a = 0; a < 3; a++){
       coeff_lst1[getSingleIdx(lnspin, i, a)] =  _coeff_list[0][idx] / 2;
       coeff_lst2[getSingleIdx(lnspin, i, a)] =  _coeff_list[1][idx] / 2;
       idx++;
      }
    }
    

    for(int i = 0; i < nterm * n_op; i++)
        op_coeff_lst[i] = 0.0;

    anorm = 0.0;
    for (int i = 0; i < nterm; i++)
      anorm += fabs(coeff_lst1[i]);
    
    if( rank == 0)
        printf("[Data] MatExp,anorm = %e, typical_time = %e\n",anorm,1.0 / anorm);
    
    #pragma omp parallel for
    for (long i = 0; i < n; i++){
      v[i] = 0.0;
      f[i] = 0.0;
      w[i] = 0.0;
      cache[i] = 0.0;
      cache[i + n] = 0.0;
      for (long j = 0; j < m + 2; j++)
        wsp[n * j + i] = 0.0;
    }
    

    #pragma omp parallel for
    for (int i = 0; i < n; i++){
      v[i] = real(_v[i]) + I*imag(_v[i]);
      f[i] = real(_f[i]) + I*imag(_f[i]);}
    
    //operator:E [*] A
    for (int i = 0; i < n_op; i++){
        if( _op_list[i][0] == 1){
            op_coeff_lst[ i * nterm + getSingleIdx(lnspin, _op_list[i][1],_op_list[i][2])] = 1.0 / 2;}
        else if( _op_list[i][0] == 2){
            op_coeff_lst[ i * nterm + getPairIdx(lnspin,_op_list[i][1],_op_list[i][2], _op_list[i][3], _op_list[i][4])] = 1.0 / 4;}
        else{
            printf("[Debug] MatExp,wrong input!!!");
            exit(0);}
    } 
    
    h1 = hamGen(lnspin, coeff_lst1);
    h2 = hamGen(lnspin, coeff_lst2);
    h3 = _coeff_list[2][0];
    
    double  *op_temp;
    op_temp = (double*)malloc(sizeof(double) * nterm);
    alst = (ham*)malloc(sizeof(ham) * n_op);
    for (int i = 0; i < n_op; i++){
        for (int j = 0; j < nterm; j++){
            op_temp[j] = op_coeff_lst[ j + (i * nterm) ];}
        alst[i] = hamGen(lnspin, op_temp);
    }
    
    double  tic, tac;
    tic = omp_get_wtime();
    
    //{{{ 
    //printf("[Debug] MatExp,n = %d, m = %d, tol = %e, anorm = %e\n",n,m,tol,anorm);   
    //printf("[Debug] MatExp,lwsp = %d, liwsp = %d, itrace = %d, iflag = %d\n",lwsp,liwsp,itrace,iflag);
    //printf("[Debug] MatExp,lnspin = %d, nterm = %d, tn = %d, n_op = %d\n",lnspin,nterm,tn,n_op);
    //}}}
    printf("[Hint] MatExp,Working ...\n");
    mpi_zcsrkexpv_(&n, &m, v, f, w, &tol, &anorm, wsp, &lwsp, iwsp, &liwsp, &itrace,\
                  &iflag, &lnspin, &h1, &h2, &h3, &prefactor, &nterm, tlst, &tn, \
                  alst, &n_op, oplst, reqs, status, cache, &nchunk);
                
    tac = omp_get_wtime();
    
    printf("[Debug] MatExp,rank %4d: mpi_zcsrkexpv works %f second!\n",rank, tac-tic);// for debug;
    
    _resVectorList = zeros<cx_mat>(tn,n_op);
    for (int j =0; j < n_op; j++){
        cx_vec res = zeros<cx_vec>(tn,1);
        for (int i = 0; i < tn; i++){
            res[i] = cx_double(creal(oplst[i+j*tn]),cimag(oplst[i+j*tn]));
        }
        _resVectorList.col(j) = res;
    }

    for (int i = 0; i < n_op; i++) hamFree(alst[i]);
    hamFree(h1);
    hamFree(h2);
    
    free(op_temp);
    free(op_coeff_lst);
    free(oplst);
    free(status);
    free(reqs);
    free(coeff_lst1);
    free(coeff_lst2);
    free(iwsp);
    free(wsp);
    free(cache);
    free(w);
    free(v);
    free(f);
    
    return _resVectorList;
}/*}}}*/

//}}}
////////////////////////////////////////////////////////////////////////////////

