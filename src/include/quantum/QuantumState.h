#ifndef QUANTUMSTATE_H
#define QUANTUMSTATE_H

#include <armadillo>

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

    cx_vec getVector() const {return _vector;};
    size_t    getDimension() const {return _dimension;};
protected:
    bool _is_pure;
    size_t _dimension;
    cx_vec _vector;
private:
};
//}}}
////////////////////////////////////////////////////////////////////////////////

/// @}
#endif
