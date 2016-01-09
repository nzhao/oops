#include "include/quantum/LiouvilleSpaceOperator.h"


////////////////////////////////////////////////////////////////////////////////
//{{{ LiouvilleSpaceOperator
LiouvilleSpaceOperator::LiouvilleSpaceOperator()
{ LOG(INFO) << "Default constructor: LiouvilleSpaceOperator";}

LiouvilleSpaceOperator::~LiouvilleSpaceOperator()
{ LOG(INFO) << "Default destructor: LiouvilleSpaceOperator";}

LiouvilleSpaceOperator::LiouvilleSpaceOperator(const HilbertSpaceOperator& op, FUNC* func)
{
    _is_expanded = true;
    _hilbert_op = op;
    _func = func;
    _kron_form = _func(_hilbert_op.getKronProdForm());
    _dim_list = _kron_form.getDimList();
    _dimension = 1;
    for( auto d : _dim_list)
        _dimension *= d;
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
