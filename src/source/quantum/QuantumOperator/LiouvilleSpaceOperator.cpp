#include "include/quantum/LiouvilleSpaceOperator.h"


////////////////////////////////////////////////////////////////////////////////
//{{{ LiouvilleSpaceOperator
LiouvilleSpaceOperator::LiouvilleSpaceOperator()
{ LOG(INFO) << "Default constructor: LiouvilleSpaceOperator";}

LiouvilleSpaceOperator::~LiouvilleSpaceOperator()
{ LOG(INFO) << "Default destructor: LiouvilleSpaceOperator";}

LiouvilleSpaceOperator::LiouvilleSpaceOperator(const HilbertSpaceOperator& op, MatExpanFunc* func)
{
    _is_expanded = true;
    _expan.op = op;
    _expan.func = func;

    _kron_form = Expand(_expan.op.getKronProdForm(), func);
    _dim_list = _kron_form.getDimList();
    _dimension = 1;
    //for( auto d : _dim_list)
    //    _dimension *= d;
    for(int i=0; i<_dim_list.size(); ++i)
        _dimension *= _dim_list[i];
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
