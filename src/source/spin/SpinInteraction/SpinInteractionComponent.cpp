#include "include/spin/SpinInteractionComponent.h"

using namespace std;

////////////////////////////////////////////////////////////////////////////////
//{{{ cSpinInteractionDomain
cSpinInteractionDomain::cSpinInteractionDomain()
{
    cout << "cSpinInteractionDomain" << endl;
}
cSpinInteractionDomain::~cSpinInteractionDomain()
{
    cout << "cSpinInteractionDomain, destructed." << endl;
}

ostream&  operator << (ostream& outs, const cSpinInteractionDomain& dm)
{
    int i=0;
    for(auto idx: dm._index_list)
    {
        outs << "interaction domain[" << i << "]: ";
        for(int j=0; j<dm._nbody; ++j)
        {
            outs << idx[j] << ", ";
        }
        outs << endl;
        i++;
    }
    return outs;
}
//}}}
//----------------------------------------------------------------------------//
//{{{ SpinPair
SpinPair::SpinPair(const vector<cSPIN>& spin_list)
{
    cout << "SpinPair constructor" << endl;
    int nspin=spin_list.size();

    _nbody = 2;

    for(int i=0; i<nspin; ++i)
        for(int j=i+1; j<nspin; ++j)
            _index_list.push_back(vector<int> {i, j});

    for(auto idx:_index_list)
        _spin_aggregate.push_back( vector<cSPIN> { spin_list[ idx[0] ] , spin_list[ idx[1] ] });
}
SpinPair::~SpinPair()
{
    cout << "SpinPair, destructed." << endl;
}
//}}}
////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
//{{{ cSpinInteractionForm
cSpinInteractionForm::cSpinInteractionForm()
{
    cout << "cSpinInteractionForm" << endl;
}
cSpinInteractionForm::~cSpinInteractionForm()
{
    cout << "cSpinInteractionForm, destructed." << endl;
}
//}}}
//----------------------------------------------------------------------------//
//{{{ TwoSpinInteractionFrom
TwoSpinInteractionFrom::TwoSpinInteractionFrom(cSpinInteractionDomain& domain)
{
    auto sag = domain.getSpinAggregate();
    for(auto it=sag.begin(); it!=sag.end(); ++it)
    {
        cSPIN spin0=(*it)[0];    cSPIN spin1=(*it)[1];
        cout << spin0.sx() * spin1.sy() << endl;
    }
}
TwoSpinInteractionFrom::~TwoSpinInteractionFrom()
{
    cout << "desctructor TwoSpinInteractionFrom" << endl;
}
//}}}
////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
//{{{ cSpinInteractionCoeff
cSpinInteractionCoeff::cSpinInteractionCoeff()
{
    cout << "cSpinInteractionCoeff" << endl;
}
cSpinInteractionCoeff::~cSpinInteractionCoeff()
{
    cout << "cSpinInteractionCoeff, destructed." << endl;
}
//}}}
////////////////////////////////////////////////////////////////////////////////
