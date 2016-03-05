#ifndef MATEXP_H
#define MATEXP_H

#include <vector>
#include <cassert>
#include <armadillo>
#include "include/misc/misc.h"
#include "include/kron/KronProd.h"
#include <numeric>      // std::partial_sum

using namespace arma;

////////////////////////////////////////////////////////////////////////////////
//{{{  MatExp
class MatExp
{
public:
    enum MatExpMethod {ArmadilloExpMat, PadeApproximation};

    MatExp(){};
    MatExp(const cx_mat& m, cx_double prefactor, MatExpMethod method); 
    ~MatExp(){};

    void   run();
    cx_mat getResultMatrix() const {return _resMatrix;};
protected:
private:
    cx_mat       _matrix;
    cx_double    _prefactor;
    MatExpMethod _method;

    cx_mat       _resMatrix;

    cx_mat armadillo_exp_mat() {return expmat(_prefactor*_matrix);};
    cx_mat pade_exp_mat();
};
//}}}
////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
//{{{  MatExpVector
class MatExpVector
{
public:
    MatExpVector() {};
    MatExpVector(const SumKronProd& skp, cx_double prefactor, const vec& time_list);
    ~MatExpVector() {};

    void run();
    vector<cx_vec> getResult() const {return _resVectorList;}; 
protected:
private:
    SumKronProd _skp;
    cx_double   _prefactor;
    vec         _time_list;
    size_t         _nTime;
    size_t         _dim; 
    
    vector<cx_vec>  _resVectorList;

    size_t _klim;
    size_t _krylov_m;
    double _krylov_tol;
    size_t _itrace;
};
//}}}
////////////////////////////////////////////////////////////////////////////////
#endif
