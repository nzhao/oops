#include "include/spin/SpinInteraction.h"

using namespace std;

////////////////////////////////////////////////////////////////////////////////
//{{{ cSpinInteraction
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
//}}}
////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
//{{{ SpinDipolarInteraction 
SpinDipolarInteraction::SpinDipolarInteraction()
{
    cout << "default SpinDipolarInteraction constructor. " << endl;
}

SpinDipolarInteraction::SpinDipolarInteraction(const vector<cSPIN>& spin_list)
{
    cout << "constructor of SpinDipolarInteraction with spin_list." << endl;

    int nspin;
    _spin_list=spin_list;

    _domain=SpinPair(spin_list);
    _form=TwoSpinInteractionForm(_domain);
    _coeff=DipolarInteractionCoeff(_domain);

//    cout << _domain << endl;
//    cout << _form << endl;
//    cout << _coeff << endl;
}

SpinDipolarInteraction::~SpinDipolarInteraction()
{
    cout << "default destruction function of SpinDipolarInteraction." << endl;
}

void SpinDipolarInteraction::make()
{
    cout << "make is called." << endl;
}
//}}}
////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
//{{{ SpinZeemanInteraction 
SpinZeemanInteraction::SpinZeemanInteraction()
{
    cout << "default SpinZeemanInteraction constructor. " << endl;
}

SpinZeemanInteraction::SpinZeemanInteraction(const vector<cSPIN>& spin_list, const vec& magB)
{
    cout << "constructor of SpinZeemanInteraction with spin_list." << endl;

    int nspin;
    _spin_list=spin_list;

    _domain=SingleSpin(spin_list);
    _form=SingleSpinInteractionForm(_domain);
    _coeff=ZeemanInteractionCoeff(_domain, magB);

    cout << _domain << endl;
    cout << _form << endl;
    cout << _coeff << endl;
}

SpinZeemanInteraction::~SpinZeemanInteraction()
{
    cout << "default destruction function of SpinZeemanInteraction." << endl;
}

void SpinZeemanInteraction::make()
{
    cout << "make is called." << endl;
}
//}}}
////////////////////////////////////////////////////////////////////////////////
