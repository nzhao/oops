#ifndef QUANTUMOPERATOR_H
#define QUANTUMOPERATOR_H
#include <armadillo>
#include "mat.h"
#include "include/easylogging++.h"
#include "include/kron/KronProd.h"

using namespace std;

/// \defgroup Quantum Quantum
/// @{

/// \defgroup QuantumOperator QuantumOperator
/// @{

////////////////////////////////////////////////////////////////////////////////
//{{{ QuantumOperator
class QuantumOperator
{
public:
    QuantumOperator();
    ~QuantumOperator();

    cx_mat       getMatrix() {return _kron_form.full();};
    SumKronProd& getKronProdForm(){return _kron_form;};
//    void         saveMatrix();
protected:
    int         _dimension;
    DIM_LIST    _dim_list;
    SumKronProd _kron_form;
private:
};
//}}}
////////////////////////////////////////////////////////////////////////////////

/// @}
/// @}
#endif
