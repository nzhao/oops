#ifndef SPININTERACTION_H
#define SPININTERACTION_H

#include <vector>
#include "include/spin/SpinInteractionComponent.h"
#include "include/spin/Spin.h"
////////////////////////////////////////////////////////////////////////////////
//{{{ KronProdForm
struct KronProdForm
{
    INDICES spin_index;
    double coeff;
    TERM mat;
};
//}}}
////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
//{{{ cSpinInteraction
class cSpinInteraction
{
public:
     cSpinInteraction();
     cSpinInteraction(const vector<cSPIN>& spin_list);
    ~cSpinInteraction();

    void make();
    vector<KronProdForm> getKronProdList(){return _kronProd_list;};

protected:
    vector<cSPIN> _spin_list;
    cSpinInteractionDomain _domain;
    cSpinInteractionForm   _form;
    cSpinInteractionCoeff  _coeff;

    vector<KronProdForm> _kronProd_list;
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
protected:
private:
};
//}}}
////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
//{{{ SpinZeemanInteraction
class SpinZeemanInteraction:public cSpinInteraction
{
public:
    SpinZeemanInteraction();
    SpinZeemanInteraction(const vector<cSPIN>& spin_list, const vec& magB);
    ~SpinZeemanInteraction();
protected:
private:
};
//}}}
////////////////////////////////////////////////////////////////////////////////
#endif
