#include "include/quantum/PureState.h"

////////////////////////////////////////////////////////////////////////////////
//{{{ PureState
PureState::PureState()
{ LOG(INFO) << "Default constructor: PureState";
    _is_pure = true;
}

PureState::~PureState()
{ LOG(INFO) << "Default destructor: PureState";}

PureState::PureState(const int& dim)
{
    _is_pure=true;
    _dimension = dim; 
    _vector = zeros<cx_vec> (dim);
}

PureState::PureState(const cx_vec& v)
{
    _is_pure=true;
    _vector = v; 
    _dimension=v.size();
};
//}}}
////////////////////////////////////////////////////////////////////////////////
