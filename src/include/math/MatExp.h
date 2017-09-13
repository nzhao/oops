#ifndef MATEXP_H
#define MATEXP_H

#include <iostream>
#include <complex>
#include <math.h>

#include <mpi.h>
#include <omp.h>
#include <complex>

#include <armadillo>

#include <vector>
#include <cassert>
#include <armadillo>
#include "include/misc/misc.h"
#include "include/kron/KronProd.h"
#include <numeric>      // std::partial_sum
#include "expv/include/expv.h"

using namespace arma;

////////////////////////////////////////////////////////////////////////////////
//{{{  MatExp
class MatExp
{
public:
    enum MatExpMethod {ArmadilloExpMat, PadeApproximation, LargeMatExpv};

    MatExp(){};
    MatExp(const cx_mat& m, cx_double prefactor, MatExpMethod method);
    MatExp(const int nspin, const vec coeff_list, const cx_vec local_phi0,\
		   const vector<uvec> op_list, const cx_double PREFACTOR,\ 
		   const vec TIME_LIST, MatExpMethod method);
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
    
    //external variable in LargeMatExpv
    cx_mat large_mat_expv();
    vec _coeff_list;
    vector<uvec> _op_list;
    int _nspin;
    cx_vec _v;
    vec _time_list;
};
//}}}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
//{{{  MatExpVector
class MatExpVector
{
public:
    enum MatExpVectorMethod {Explicit, ExplicitSparse, Inexplicit, InexplicitGPU, LargeVecExpv};

    MatExpVector() {};
    MatExpVector(const SumKronProd& skp, const cx_vec& v, const vec& time_list, MatExpVectorMethod method);
    MatExpVector(const cx_mat& m, const cx_vec& v, const cx_double prefactor, const vec& time_list);
    MatExpVector(const sp_cx_mat& m, const cx_vec& v, const cx_double prefactor, const vec& time_list);
    MatExpVector(const int nspin, const vec coeff_list, const cx_vec local_phi0,const vector<uvec> op_list,\
           const cx_vec final_vec, const cx_double PREFACTOR,const vec TIME_LIST, MatExpVectorMethod method);
    ~MatExpVector() {};

    cx_mat run();
    cx_mat runExplicit();
    cx_mat runExplicitSparse();
    cx_mat runInexplicit();
    cx_mat runInexplicitGPU();
    cx_mat runLargeVecExpv();
    cx_mat getResult() const {return _resVectorList;}; 

    void enable_step_print() {_itrace = 1;}
    void enable_skp_print() {_is_print_skp = true;}
    void disable_step_print() {_itrace = 0;}
    void disable_skp_print() {_is_print_skp = false;}
protected:
private:
    MatExpVectorMethod _method;
    SumKronProd _skp;
    cx_mat      _matrix;
    sp_cx_mat   _sp_matrix;
    cx_vec      _vector;
    cx_double   _prefactor;
    vec         _time_list;
    size_t         _nTime;
    size_t         _dim; 
    
    cx_mat  _resVectorList;

    size_t _klim;
    size_t _krylov_m;
    double _krylov_tol;
    size_t _itrace;

    bool _is_print_skp;
    
    //external variable in LargeMatExpv
    vec _coeff_list;
    vector<uvec> _op_list;
    int _nspin;
    cx_vec _v;
    cx_vec _f;
};
//}}}
////////////////////////////////////////////////////////////////////////////////
#endif
