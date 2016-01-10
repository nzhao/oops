#ifndef LIOUVILLESPACEOPERATOR_H
#define LIOUVILLESPACEOPERATOR_H

#include "include/quantum/QuantumOperator.h"
#include "Include/quantum/HilbertSpaceOperator.h"

/// \addtogroup QuantumOperator 
/// @{

/// \defgroup LiouvilleSpaceOperator LiouvilleSpaceOperator
/// @{

struct EXPAN 
{
    HilbertSpaceOperator op;
    MatExpanFunc* func;
};


////////////////////////////////////////////////////////////////////////////////
//{{{ LiouvilleSpaceOperator
class LiouvilleSpaceOperator:public QuantumOperator
{
public:
    LiouvilleSpaceOperator();
    ~LiouvilleSpaceOperator();

    LiouvilleSpaceOperator(const HilbertSpaceOperator& op, MatExpanFunc* func);
protected:
private:
    bool _is_expanded;
    EXPAN _expan;
};
//}}}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
//{{{ Liouvillian
class Liouvillian:public LiouvilleSpaceOperator
{
public:
    Liouvillian();
    Liouvillian(const Hamiltonian& hm):LiouvilleSpaceOperator(hm, CIRCLEC){};
    ~Liouvillian();
protected:
private:
};
//}}}
////////////////////////////////////////////////////////////////////////////////
/// @}
/// @}
#endif
