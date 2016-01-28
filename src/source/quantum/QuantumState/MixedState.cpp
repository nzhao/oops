#include "include/quantum/MixedState.h"

////////////////////////////////////////////////////////////////////////////////
//{{{ MixedState 

MixedState::MixedState()
{ LOG(INFO) << "Default constructor: MixedState";}

MixedState::~MixedState()
{ LOG(INFO) << "Default destructor: MixedState";}
//}}}
////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
//{{{ DensityOperator
DensityOperator::DensityOperator()
{ LOG(INFO) << "Default constructor: DensityOperator";}

DensityOperator::~DensityOperator()
{ LOG(INFO) << "Default destructor: DensityOperator";}

DensityOperator::DensityOperator(const vector<cSPIN>& spin_list)
{
    _op = HilbertSpaceOperator(spin_list);
    _dimension = _op.getDimension();
}

//}}}
////////////////////////////////////////////////////////////////////////////////
