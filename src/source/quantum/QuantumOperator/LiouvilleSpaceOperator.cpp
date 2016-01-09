#include "include/quantum/LiouvilleSpaceOperator.h"


////////////////////////////////////////////////////////////////////////////////
//{{{ LiouvilleSpaceOperator
LiouvilleSpaceOperator::LiouvilleSpaceOperator()
{ LOG(INFO) << "Default constructor: LiouvilleSpaceOperator";}

LiouvilleSpaceOperator::~LiouvilleSpaceOperator()
{ LOG(INFO) << "Default destructor: LiouvilleSpaceOperator";}

LiouvilleSpaceOperator::LiouvilleSpaceOperator(const HilbertSpaceOperator& op, const OperatorExapnsionMethod method)
{
    _is_expanded = true;
    //_hilbert_op = op;
    _expansion_method = method;
}
//}}}
////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
//{{{ Liouvillian
Liouvillian::Liouvillian()
{ LOG(INFO) << "Default constructor: Liouvillian";}

Liouvillian::~Liouvillian()
{ LOG(INFO) << "Default destructor: Liouvillian";}

//}}}
////////////////////////////////////////////////////////////////////////////////
