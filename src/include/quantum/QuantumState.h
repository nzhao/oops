#ifndef QUANTUMSTATE_H
#define QUANTUMSTATE_H

#include <armadillo>
#include "include/easylogging++.h"

using namespace arma;
/// \addtogroup Quantum 
/// @{

/// \defgroup QuantumState QuantumState
/// @{


////////////////////////////////////////////////////////////////////////////////
//{{{ QuantumState
class QuantumState
{
public:
    QuantumState();
    ~QuantumState();

    cx_vec getVector(){return _vector;};
protected:
    bool _is_pure;
    int _dimension;
    cx_vec _vector;
private:
};
//}}}
////////////////////////////////////////////////////////////////////////////////

/// @}
#endif
