#ifndef LIOUVILLESPACEOPERATOR_H
#define LIOUVILLESPACEOPERATOR_H

#include "include/quantum/QuantumOperator.h"
#include "Include/quantum/HilbertSpaceOperator.h"

/// \addtogroup QuantumOperator 
/// @{

/// \defgroup LiouvilleSpaceOperator LiouvilleSpaceOperator
/// @{

enum OperatorExapnsionMethod { SHARP = 0, FLATTEN = 1};

////////////////////////////////////////////////////////////////////////////////
//{{{ LiouvilleSpaceOperator
class LiouvilleSpaceOperator:public QuantumOperator
{
public:
    LiouvilleSpaceOperator();
    ~LiouvilleSpaceOperator();

    LiouvilleSpaceOperator(const HilbertSpaceOperator& lOp, const OperatorExapnsionMethod method);
    LiouvilleSpaceOperator(const HilbertSpaceOperator& lOp, const HilbertSpaceOperator& rOp);
protected:
private:
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
