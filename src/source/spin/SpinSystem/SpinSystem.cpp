#include <vector>
#include <assert.h>
#include "include/spin/SpinSystem.h"


////////////////////////////////////////////////////////////////////////////////
//{{{ cSpinSystem
cSpinSystem::cSpinSystem()
{
    cout << "default constructor of cSpinSystem." << endl;
}
cSpinSystem::cSpinSystem(const vector<cSPIN>& spin_list)
{
    cout << "constructor of cSpinSystem with spin_list." << endl;
    _spin_list = spin_list;

    for(auto spin: _spin_list)
        _dim_list.push_back( spin.get_dimension() );
}
cSpinSystem::~cSpinSystem()
{
    cout << "destructor of cSpinSystem." << endl;
}

void cSpinSystem::addSpinInteraction(cSpinInteraction& spin_interaction)
{
    assert( this->_dim_list == spin_interaction.getDimList() );
    _spin_interaction_list.push_back(spin_interaction);
}

SumKronProd cSpinSystem::getSumKronOperator()
{
    SumKronProd res=_spin_interaction_list[0].getSumKronProd();
    for(int i=1; i<_spin_interaction_list.size(); ++i)
        res=res+_spin_interaction_list[i].getSumKronProd();
    return res;
}

//}}}
////////////////////////////////////////////////////////////////////////////////
