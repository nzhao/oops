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
SpinPair::SpinPair(int nspin)
{
    cout << "SpinPair constructor" << endl;

    _nbody = 2;

    for(int i=0; i<nspin; ++i)
        for(int j=i+1; j<nspin; ++j)
            _index_list.push_back(vector<int> {i, j});
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
