#ifndef LIOUVILLESPACEOPERATOR_H
#define LIOUVILLESPACEOPERATOR_H

#include "include/quantum/QuantumOperator.h"
#include "Include/quantum/HilbertSpaceOperator.h"

/// \addtogroup QuantumOperator 
/// @{

/// \defgroup LiouvilleSpaceOperator LiouvilleSpaceOperator
/// @{

enum OperatorExapnsionMethod { SHARP = 0, FLATTEN = 1, CircleC = 2};

////////////////////////////////////////////////////////////////////////////////
//{{{ LiouvilleSpaceOperator
class LiouvilleSpaceOperator:public QuantumOperator
{
public:
    LiouvilleSpaceOperator();
    ~LiouvilleSpaceOperator();

    LiouvilleSpaceOperator(const HilbertSpaceOperator& op, FUNC* func);
    void makeKronForm(){};
protected:
private:
    bool _is_expanded;
    HilbertSpaceOperator _hilbert_op;
    FUNC * _func;
};
//}}}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
//{{{ Liouvillian
class Liouvillian
{
public:
    Liouvillian();
    ~Liouvillian();
protected:
private:
};
//}}}
////////////////////////////////////////////////////////////////////////////////
/// @}
/// @}
#endif
