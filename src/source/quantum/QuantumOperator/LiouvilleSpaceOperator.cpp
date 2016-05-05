#include "include/quantum/LiouvilleSpaceOperator.h"
#include "include/misc/misc.h"


////////////////////////////////////////////////////////////////////////////////
//{{{ LiouvilleSpaceOperator
LiouvilleSpaceOperator::LiouvilleSpaceOperator()
{ //LOG(INFO) << "Default constructor: LiouvilleSpaceOperator";
}

LiouvilleSpaceOperator::~LiouvilleSpaceOperator()
{ //LOG(INFO) << "Default destructor: LiouvilleSpaceOperator";
}

LiouvilleSpaceOperator::LiouvilleSpaceOperator(const vector<cSPIN>& spin_list)
{
    _spin_list = spin_list;
    _dimension = 1;
    for(int i=0; i<_spin_list.size(); ++i)
    {
        int dim = _spin_list[i].get_dimension2(); 
        _dim_list.push_back( dim );
        _dimension *= dim;
    }
}

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

void LiouvilleSpaceOperator::addInteraction(cSpinInteraction& spin_interaction)
{
    //cout << "dim_list = ";
    //print_vector(_dim_list);
    //cout << endl;
    //cout << "interaction dim_list = ";
    //print_vector(spin_interaction.getDimList());
    //cout << endl;
    assert( this->_dim_list == spin_interaction.getDimList() );
    _interaction_list.push_back(spin_interaction);
}

void LiouvilleSpaceOperator::make()
{
    _kron_form = _interaction_list[0].getSumKronProd();
    for(int i=1; i<_interaction_list.size(); ++i)
        _kron_form = _kron_form + _interaction_list[i].getSumKronProd();
}
//}}}
////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
//{{{ Liouvillian
Liouvillian::Liouvillian()
{ //LOG(INFO) << "Default constructor: Liouvillian";
}

Liouvillian::Liouvillian(const QuantumOperator& op)
{
    _is_expanded = false;
    _dimension = op.getDimension();
    _dim_list = op.getDimList();
    _kron_form = op.getKronProdForm();
}

Liouvillian::~Liouvillian()
{ //LOG(INFO) << "Default destructor: Liouvillian";
}

//}}}
////////////////////////////////////////////////////////////////////////////////
