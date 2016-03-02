#ifndef MATEXP_H
#define MATEXP_H

#include <vector>
#include <cassert>
#include <armadillo>
#include "include/kron/KronProd.h"

using namespace arma;

class MatExp
{
public:
    enum MatExpMethod {ArmadilloExpMat, PadeApproximation};

    MatExp(){};
    MatExp(const cx_mat& m, cx_double prefactor, MatExpMethod method); 
    ~MatExp(){};

    void   run();
    cx_mat getResultMatrix() {return _resMatrix;};
protected:
private:
    cx_mat       _matrix;
    cx_double    _prefactor;
    MatExpMethod _method;

    cx_mat       _resMatrix;

    cx_mat armadillo_exp_mat() {return expmat(_prefactor*_matrix);};
    cx_mat pade_exp_mat();
};
#endif
