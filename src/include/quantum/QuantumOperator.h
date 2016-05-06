#ifndef QUANTUMOPERATOR_H
#define QUANTUMOPERATOR_H
#ifdef HAS_MATLAB
#include <mat.h>
#endif

#include <armadillo>
#include "include/kron/KronProd.h"

using namespace std;
extern string DEBUG_PATH;

/// \defgroup Quantum Quantum
/// @{

/// \defgroup QuantumOperator QuantumOperator
/// @{

////////////////////////////////////////////////////////////////////////////////
//{{{ QuantumOperator
class QuantumOperator
{
public:
    QuantumOperator() {};
    ~QuantumOperator() {};

    cx_mat       getMatrix() {return _kron_form.full();};
    SumKronProd  getKronProdForm() const  {return _kron_form;};
    DIM_LIST     getDimList() const {return _dim_list;};
    int          getDimension() const {return _dimension;};
    void         saveMatrix(string filename);
    QuantumOperator& scale(double factor) {_kron_form.scale(factor); return *this;};

    friend QuantumOperator operator + (const QuantumOperator& op1, const QuantumOperator& op2);
    friend QuantumOperator operator - (const QuantumOperator& op1, const QuantumOperator& op2);
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
