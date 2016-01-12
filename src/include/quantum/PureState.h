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
    PureState(const cx_vec& v) {_vector = v; _is_pure=true;};
    ~PureState();

    cx_vec getVector(){return _vector;};

    void setVector(const cx_vec& v) {_vector = v;};
    void setComponent(int i, cx_double val) {_vector(i) = val;};
protected:
private:
    cx_vec _vector;
};
//}}}
////////////////////////////////////////////////////////////////////////////////

/// @}
/// @}
#endif
