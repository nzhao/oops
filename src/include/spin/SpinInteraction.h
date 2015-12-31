#ifndef SPININTERACTION_H
#define SPININTERACTION_H

#include <vector>
#include "include/spin/SpinInteractionComponent.h"
#include "include/spin/Spin.h"

////////////////////////////////////////////////////////////////////////////////
//{{{ cSpinInteraction
class cSpinInteraction
{
public:
     cSpinInteraction();
     cSpinInteraction(const vector<cSPIN>& spin_list);
    ~cSpinInteraction();

    virtual void make() = 0;

protected:
    vector<cSPIN> _spin_list;
    cSpinInteractionDomain _domain;
    cSpinInteractionForm   _form;
    cSpinInteractionCoeff  _coeff;

private:
};
//}}}
////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
//{{{ SpinDipolarInteraction
class SpinDipolarInteraction:public cSpinInteraction
{
public:
    SpinDipolarInteraction();
    SpinDipolarInteraction(const vector<cSPIN>& spin_list);
    ~SpinDipolarInteraction();
    void make();
protected:
private:
};
//}}}
////////////////////////////////////////////////////////////////////////////////
#endif
