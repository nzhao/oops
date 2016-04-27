#ifndef SPKRON_H
#define SPKRON_H
#include <iostream>
#include <armadillo>
//#include "KronProd.h"
#include <functional>   // std::multiplies
#include <numeric>      // std::accumulate

using namespace std;
using namespace arma;

////////////////////////////////////////////////////////////////////////////////
//{{{ cooSpMat
template <typename T>
class cooSpMat
{
public:
    cooSpMat() {};
    cooSpMat(const int dim);
    cooSpMat(const Mat<T>& m);
    cooSpMat(const SpMat<T>& m);
    cooSpMat(const int nnz, const int nrow, const int ncol, const unsigned int * ia, const unsigned int * ja, const T * a);
    ~cooSpMat();

    int nnz() const {return _nnz;}
    int nrow() const {return _n_rows;}
    int ncol() const {return _n_cols;}
    unsigned int * rIdx() const {return _ia;}
    unsigned int * cIdx() const {return _ja;}
    unsigned int row(int i) const {return _ia[i];}
    unsigned int col(int i) const {return _ja[i];}
    T val(int i) const {return _a[i];}
    T * pVal() const {return _a;}

    void copy(const cooSpMat<T>& m);

    SpMat<T> mat() const;
    Mat<T>   full_mat() const;

    void print() const { for(int i=0; i<_nnz; ++i) cout << i << " " << _ia[i] << ", " << _ja[i] << " : " << _a[i] << endl;}
private:
    int _n_cols;
    int _n_rows;
    int _nnz;
    unsigned int * _ia;
    unsigned int * _ja;
    T * _a;
};
//}}}
////////////////////////////////////////////////////////////////////////////////

//template<typename T> 
//cooSpMat<T> spkron(const cooSpMat<T>& A, const cooSpMat<T>& B);

//template<typename T>
//cooSpMat<T> spkron(const int dim, const cooSpMat<T>& B);

//template<typename T> 
//cooSpMat<T> spkron(const cooSpMat<T>& A, const int dim);

//vector<size_t> dim_acc(const vector<size_t>& dim_list, const vector<size_t>& idx_list);

//template<typename T>
//cooSpMat<T> spkron(const vector<size_t>& dim_list, const vector<size_t>& idx_list, const vector< Mat<T> >& mat_list);


template<typename T>  inline
cooSpMat<T> spkron(const cooSpMat<T>& A, const cooSpMat<T>& B)
{/*{{{*/
    int nrow = A.nrow() * B.nrow();
    int ncol = A.ncol() * B.ncol();
    int nnz = A.nnz() * B.nnz(); 

    unsigned int * ia = new unsigned int [nnz];
    unsigned int * ja = new unsigned int [nnz];
    T * a = new T [nnz];

    int idx = 0;
    for(int i=0; i<A.nnz(); ++i)
    {
        for(int j=0; j<B.nnz(); ++j)
        {
            ia[idx] = A.row(i)*B.nrow() + B.row(j);
            ja[idx] = A.col(i)*B.ncol() + B.col(j);
            a[idx] = A.val(i)*B.val(j);
            idx ++;
        }
    }
    cooSpMat<T> res(nnz, nrow, ncol, ia, ja, a);

    delete[] ia;
    delete[] ja;
    delete[] a;

    return res;
}/*}}}*/

template<typename T> inline
cooSpMat<T> spkron(const int dim, const cooSpMat<T>& B)
{/*{{{*/
    int nrow = dim * B.nrow();
    int ncol = dim * B.ncol();
    int nnz = dim * B.nnz(); 

    unsigned int * ia = new unsigned int [nnz];
    unsigned int * ja = new unsigned int [nnz];
    T * a = new T [nnz];

    int idx = 0;
    for(int i=0; i<dim; ++i)
    {
        for(int j=0; j<B.nnz(); ++j)
        {
            ia[idx] = i*B.nrow() + B.row(j);
            ja[idx] = i*B.ncol() + B.col(j);
            a[idx] = B.val(j);
            idx ++;
        }
    }
    cooSpMat<T> res(nnz, nrow, ncol, ia, ja, a);

    delete[] ia;
    delete[] ja;
    delete[] a;

    return res;
}/*}}}*/

template<typename T>  inline
cooSpMat<T> spkron(const cooSpMat<T>& A, const int dim)
{/*{{{*/
    int nrow = A.nrow() * dim;
    int ncol = A.ncol() * dim;
    int nnz = A.nnz() * dim; 

    unsigned int * ia = new unsigned int [nnz];
    unsigned int * ja = new unsigned int [nnz];
    T * a = new T [nnz];

    int idx = 0;
    for(int i=0; i<A.nnz(); ++i)
    {
        for(int j=0; j<dim; ++j)
        {
            ia[idx] = A.row(i)*dim + j;
            ja[idx] = A.col(i)*dim + j;
            a[idx] = A.val(i);
            idx ++;
        }
    }
    cooSpMat<T> res(nnz, nrow, ncol, ia, ja, a);

    delete[] ia;
    delete[] ja;
    delete[] a;

    return res;
}/*}}}*/

inline vector<size_t> dim_acc(const vector<size_t>& dim_list, const vector<size_t>& idx_list)
{/*{{{*/
    vector<size_t> res;
    
    vector<size_t> idx_list1(idx_list);
    idx_list1.push_back(dim_list.size());

    int p0 = 0;
    for(int i=0; i<idx_list1.size(); ++i)
    {
        res.push_back( accumulate(dim_list.begin() + p0, dim_list.begin() + idx_list1[i], 1, multiplies<size_t>()) );
        p0 = idx_list1[i]+1;
    }
    return res;
}/*}}}*/

template<typename T> inline
cooSpMat<T> spkron(const vector<size_t>& dim_list, const vector<size_t>& idx_list, const vector< Mat<T> >& mat_list)
{/*{{{*/
    vector<size_t> id_dim_list = dim_acc(dim_list, idx_list);
    
    cooSpMat<T> res(id_dim_list[0]);
    for(int i=0; i<mat_list.size(); ++i)
    {
        cooSpMat<T> kronM1 = spkron( spkron( res, cooSpMat<T>(mat_list[i])), id_dim_list[i+1] );
        res.copy( kronM1 );
    }
    return res;
}/*}}}*/
#endif
