#ifndef SPININTERACTION_H
#define SPININTERACTION_H

#include <vector>
#include "include/el/easylogging++.h"
#include "include/spin/Spin.h"
#include "include/spin/SpinInteractionComponent.h"
#include "include/spin/SpinInteractionDefine.h"
#include "include/kron/KronProd.h"
#include "include/quantum/PureState.h"

/// \addtogroup Spin
/// @{

/// \defgroup SpinInteraction SpinInteraction
/// @{

/// \defgroup SpinInteractionClasses SpinInteractionClasses
/// @{
/// This class defines an abstract class. 
/// cSpinInteraction.
///
////////////////////////////////////////////////////////////////////////////////
//{{{ cSpinInteraction
class cSpinInteraction
{
public:
    enum Space {Hilbert, Liouville};

     cSpinInteraction();
     cSpinInteraction(const vector<cSPIN>& spin_list);
     cSpinInteraction(const vector<cSPIN>& spin_list, Space sp);
    ~cSpinInteraction();

    void make();
    SumKronProd& getSumKronProd(){return _sum_kron_prod;};
    DIM_LIST getDimList() {return _dim_list;};

    friend ostream&  operator << (ostream& outs, cSpinInteraction& spin_interaction);
protected:
    vector<cSPIN> _spin_list;
    cSpinInteractionDomain _domain;
    cSpinInteractionForm   _form;
    cSpinInteractionCoeff  _coeff;
    SumKronProd _sum_kron_prod;
    DIM_LIST _dim_list;    
    Space _space;
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



////////////////////////////////////////////////////////////////////////////////
//{{{ DipolarField 
class DipolarField:public cSpinInteraction
{
public:
    DipolarField();
    DipolarField(const vector<cSPIN>& spin_list, const cSPIN& center_spin, const PureState& state);
    DipolarField(const vector<cSPIN>& spin_list, const vector<cSPIN>& source_list, const vector<PureState>& state_list);
    DipolarField(const vector<cSPIN>& spin_list, const vector<cSPIN>& source_list, const vector<PureState>& state_list, const uvec& exclude_idx);
    ~DipolarField();
protected:
private:
};
//}}}
////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
//{{{ SpinDephasing
class SpinDephasing:public cSpinInteraction
{
public:
    SpinDephasing() {};
    SpinDephasing(const vector<cSPIN>& spin_list, const double dephasing_rate, const vec& axis);
    ~SpinDephasing() {};
protected:
private:
};
//}}}
////////////////////////////////////////////////////////////////////////////////

/// @}

/// @}
/// @}
#endif
