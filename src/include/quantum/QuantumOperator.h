#ifndef QUANTUMOPERATOR_H
#define QUANTUMOPERATOR_H
#include <armadillo>
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

    virtual void makeKronForm()=0;
    
    cx_mat       fullMatrix() {return _kron_form.full();};
    SumKronProd& getKronProdForm(){return _kron_form;};
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
