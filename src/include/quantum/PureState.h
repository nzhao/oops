#ifndef PURESTATE_H
#define PURESTATE_H
#include <armadillo>
#include "include/quantum/QuantumState.h"

using namespace arma;

/// \addtogroup QuantumState PureState
/// @{

/// \defgroup PureState PureState
/// @{
////////////////////////////////////////////////////////////////////////////////
//{{{ PureState
class PureState:public QuantumState
{
public:
    PureState();
    PureState(const int& dim);
    PureState(const cx_vec& v);
    ~PureState();

    void setVector(const cx_vec& v) {_vector = v;};
    void setComponent(int i, cx_double val) {_vector(i) = val;};
protected:
private:
};
//}}}
////////////////////////////////////////////////////////////////////////////////

/// @}
/// @}
#endif
