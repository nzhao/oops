#ifndef QUANTUMSTATE_H
#define QUANTUMSTATE_H

#include "include/easylogging++.h"
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
protected:
    bool _is_pure;
    int _dimension;
private:
};
//}}}
////////////////////////////////////////////////////////////////////////////////

/// @}
#endif
