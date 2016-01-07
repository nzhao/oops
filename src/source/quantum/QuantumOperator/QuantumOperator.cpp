#include "include/quantum/QuantumOperator.h"

////////////////////////////////////////////////////////////////////////////////
//{{{ QuantumOperator
QuantumOperator::QuantumOperator()
{ LOG(INFO) << "Default constructor: QuantumOperator."; }

QuantumOperator::~QuantumOperator()
{ LOG(INFO) << "Default destructor: QuantumOperator."; }

cx_mat QuantumOperator::fullMatrix()
{
    return _kron_form.full();
}
//}}}
////////////////////////////////////////////////////////////////////////////////
