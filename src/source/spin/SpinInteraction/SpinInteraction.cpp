#include "include/spin/SpinInteraction.h"

using namespace std;

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// cSpinInteraction
cSpinInteraction::cSpinInteraction()
{
    cout << "default construction function of SpinInteraction." << endl;
}

cSpinInteraction::cSpinInteraction(const vector<cSPIN>& spin_list)
{
    cout << "construction function of SpinInteraction with spin_list." << endl;
    _spin_list=spin_list;
}

cSpinInteraction::~cSpinInteraction()
{
    cout << "default destruction function of SpinInteraction." << endl;
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// SpinDipolarInteraction 
SpinDipolarInteraction::SpinDipolarInteraction()
{
    cout << "default SpinDipolarInteraction constructor. " << endl;
}

SpinDipolarInteraction::SpinDipolarInteraction(const vector<cSPIN>& spin_list)
{
    cout << "constructor of SpinDipolarInteraction with spin_list." << endl;
    _spin_list=spin_list;
    _domain=SpinPair(_spin_list.size());

    cout << _domain << endl;
}

SpinDipolarInteraction::~SpinDipolarInteraction()
{
    cout << "default destruction function of SpinDipolarInteraction." << endl;
}

void SpinDipolarInteraction::make()
{
    cout << "make is called." << endl;
}
