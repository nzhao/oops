#ifndef SPININTERACTION_H
#define SPININTERACTION_H

#include <vector>
#include "include/easylogging++.h"
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
     cSpinInteraction();
     cSpinInteraction(const vector<cSPIN>& spin_list);
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
/// @}

/// @}
/// @}
#endif
