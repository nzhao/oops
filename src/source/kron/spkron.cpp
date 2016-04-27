#include "include/kron/spkron.h"

////////////////////////////////////////////////////////////////////////////////
//{{{ CooSpMat
template<typename T>
cooSpMat<T>::cooSpMat(const int dim)
{/*{{{*/
    _nnz = dim;
    _n_cols = dim;
    _n_rows = dim;
    _ia = new unsigned int [_nnz];
    _ja = new unsigned int [_nnz];
    _a = new T [_nnz];
    for(int i=0; i<dim; ++i)
    {
        _ia[i] = i;
        _ja[i] = i;
        _a[i] = 1;
    }
}/*}}}*/

template<typename T>
cooSpMat<T>::cooSpMat(const Mat<T>& full_mat)
{/*{{{*/
    SpMat<T> m(full_mat);
    typename SpMat<T>::const_iterator start = m.begin();
    typename SpMat<T>::const_iterator end   = m.end();
    _nnz = m.n_nonzero;
    _n_cols = m.n_cols;
    _n_rows = m.n_rows;
    _ia = new unsigned int [_nnz];
    _ja = new unsigned int [_nnz];
    _a = new T [_nnz];
    for(typename SpMat<T>::const_iterator it = start; it != end; ++it)
    {
        int i = distance(start, it);
        _ia[i] = it.row();
        _ja[i] = it.col();
        _a[i] = (*it);
    }
}/*}}}*/

template<typename T>
cooSpMat<T>::cooSpMat(const SpMat<T>& m)
{/*{{{*/
    typename SpMat<T>::const_iterator start = m.begin();
    typename SpMat<T>::const_iterator end   = m.end();
    _nnz = m.n_nonzero;
    _n_cols = m.n_cols;
    _n_rows = m.n_rows;
    _ia = new unsigned int [_nnz];
    _ja = new unsigned int [_nnz];
    _a = new T [_nnz];
    for(typename SpMat<T>::const_iterator it = start; it != end; ++it)
    {
        int i = distance(start, it);
        _ia[i] = it.row();
        _ja[i] = it.col();
        _a[i] = (*it);
    }
}/*}}}*/

template<typename T>
cooSpMat<T>::cooSpMat(const int nnz, const int nrow, const int ncol, const unsigned int * ia, const unsigned int * ja, const T * a)
{/*{{{*/
    _nnz = nnz;
    _n_rows = nrow;
    _n_cols = ncol;
    _ia = new unsigned int [_nnz];
    _ja = new unsigned int [_nnz];
    _a = new T [_nnz];
    for(int i=0; i<nnz; ++i)
    {
        _ia[i] = ia[i];
        _ja[i] = ja[i];
        _a[i] = a[i];
    }
}/*}}}*/

template<typename T>
cooSpMat<T>::~cooSpMat()
{/*{{{*/
    if(_ia != NULL) delete[] _ia;
    if(_ja != NULL) delete[] _ja;
    if(_a != NULL) delete[] _a;
}/*}}}*/

template<typename T>
SpMat<T> cooSpMat<T>::mat() const
{/*{{{*/
    urowvec rIdx(_ia, _nnz);
    urowvec cIdx(_ja, _nnz);
    umat locations = join_vert(rIdx, cIdx);
    Col<T>  values(_a, _nnz);
    SpMat<T> res(locations, values, _n_rows, _n_cols);
    return res;
}/*}}}*/

template<typename T>
Mat<T> cooSpMat<T>::full_mat() const
{/*{{{*/
    Mat<T> res= conv_to< Mat<T> >::from( mat() );
    return res;
}/*}}}*/

template<typename T>
void cooSpMat<T>::copy(const cooSpMat<T>& m)
{/*{{{*/
    if(_ia != NULL) delete[] _ia;
    if(_ja != NULL) delete[] _ja;
    if(_a != NULL) delete[] _a;

    _n_cols = m.ncol();
    _n_rows = m.nrow();
    _nnz = m.nnz();
    _ia = new unsigned int [_nnz];
    _ja = new unsigned int [_nnz];
    _a = new T [_nnz];

    memcpy(_ia, m.rIdx(), _nnz*sizeof(unsigned int));
    memcpy(_ja, m.cIdx(), _nnz*sizeof(unsigned int));
    memcpy(_a, m.pVal(), _nnz*sizeof(T));
}/*}}}*/

//}}}
////////////////////////////////////////////////////////////////////////////////

template class cooSpMat<cx_double>;

////{{{ application
//sp_mat randomMat(int dim, double threshold)
//{[>{{{<]
    //mat A = randu<mat>(dim,dim);
    //uvec posA = find(A<threshold);
    //A.elem(posA) = zeros(1, posA.size() );
    //return sp_mat(A);
//}[>}}}<]

//sp_cx_mat randomMat(int dim, cx_double threshold)
//{[>{{{<]
    //sp_mat A = randomMat(dim, threshold.real() );
    //sp_mat B = randomMat(dim, threshold.real() );
    //return sp_cx_mat(A, B);
//}[>}}}<]

//typedef cx_double Tp;
//int main()
//{
    //int dim = 20; Tp threshold = 0.9;

    //////////////////////////////////////////////////////////////////////////////////
    ////{{{ test kron 1
    //SpMat<Tp> spA = randomMat(dim, threshold);
    //SpMat<Tp> spB = randomMat(dim, threshold);
    //Mat<Tp> A= conv_to< Mat<Tp> >::from(spA);
    //Mat<Tp> B= conv_to< Mat<Tp> >::from(spB);
    //Mat<Tp> Id; Id.eye(dim, dim);
    
    //cout << "dir kron" << endl;
    //SpMat<Tp> spC(kron(A, B));
    //SpMat<Tp> spCI(kron(A, Id));
    //SpMat<Tp> spIC(kron(Id, B));


    //cooSpMat<Tp> spDataA(spA);
    //cooSpMat<Tp> spDataB(spB);
    //cout << "sp kron" << endl;
    //cooSpMat<Tp> spDataC = spkron(spDataA, spDataB);
    //cooSpMat<Tp> spDataCI = spkron(spDataA, dim);
    //cooSpMat<Tp> spDataIC = spkron(dim, spDataB);
    ////spDataC.print();

    //cout << "diff = " << endl;
    //cout << spC-spDataC.mat() << endl;
    //cout << spCI-spDataCI.mat() << endl;
    //cout << spIC-spDataIC.mat() << endl;
    ////}}}
    //////////////////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////////////////
    ////{{{ test kron 2
    //int nspin = 6; vector<int> dim_list;
    //for(int i=0; i<nspin; ++i)
        //dim_list.push_back( 2 );

    //for(int i=0; i<nspin; ++i)
        //cout << i << ", ";
    //cout << endl;
    //for(int i=0; i<nspin;++i)
        //cout << dim_list[i] << ", ";
    //cout << endl;

    //Mat<Tp> m; m << 0.0 << 2.0 << endr << 2.0 << 0.0;
    //vector<int> idx_list;
    //vector< Mat<Tp> > mat_list;
    //idx_list.push_back( 0 ); idx_list.push_back( 5 );
    //mat_list.push_back( m ); mat_list.push_back( m );

    //vector<int> dimAcc = dim_acc( dim_list, idx_list);

    //for(int i=0; i<dimAcc.size(); ++i)
        //cout << dimAcc[i] << ", ";
    //cout << endl;

    //cooSpMat<Tp> res = spkron(dim_list, idx_list, mat_list);
    //cout << res.mat() << endl;

    ////}}}
    //////////////////////////////////////////////////////////////////////////////////
    
    
    //return 0;
//}
////}}}

