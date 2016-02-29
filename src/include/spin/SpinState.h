#ifndef SPINSTATE_H
#define SPINSTATE_H

#include "include/spin/Spin.h"
#include "include/spin/SpinInteraction.h"

/// \addtogroup Spin SpinState
/// @{
/// \defgroup SpinState SpinState
/// @{

////////////////////////////////////////////////////////////////////////////////
//{{{ cSpinState
class cSpinState:public cSpinInteraction
{
public:
    cSpinState();
    cSpinState(const vector<cSPIN>& spin_list):cSpinInteraction(spin_list){};
    ~cSpinState();

    void make();
protected:
private:
};
//}}}
////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
//{{{ SpinPolarization
class SpinPolarization:public cSpinState
{
public:
    SpinPolarization();
    SpinPolarization(const vector<cSPIN>& spin_list, const vector<int>& index_list, const vector<vec>& pol_list);
    SpinPolarization(const vector<cSPIN>& spin_list, const vec& pol);
    ~SpinPolarization();
protected:
private:
};
//}}}
////////////////////////////////////////////////////////////////////////////////



/// @}
/// @}
#endif
