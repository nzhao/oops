#include "include/quantum/PureState.h"

////////////////////////////////////////////////////////////////////////////////
//{{{ PureState
PureState::PureState()
{ //LOG(INFO) << "Default constructor: PureState";
    _is_pure = true;
}

PureState::~PureState()
{ //LOG(INFO) << "Default destructor: PureState";
}

PureState::PureState(const int dim)
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

PureState::PureState(const cSPIN& spin)
{
    _is_pure=true;
    _dimension = spin.get_multiplicity();
    _vector= zeros<cx_vec> (_dimension);
}
PureState::PureState(const vector<cSPIN>& spin_list)
{
    _is_pure=true;
    _dimension = 1;
    for(int i=0; i<spin_list.size(); ++i)
        _dimension *= spin_list[i].get_multiplicity();
    _vector= zeros<cx_vec> (_dimension);
}
//}}}
////////////////////////////////////////////////////////////////////////////////
