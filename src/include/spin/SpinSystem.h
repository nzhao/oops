#ifndef SPINSYSTEM_H
#define SPINSYSTEM_H
#include "include/spin/spin.h"
#include "include/spin/SpinInteraction.h"
#include "include/spin/SpinInteractionDefine.h"
#include "include/kron/KronProd.h"

////////////////////////////////////////////////////////////////////////////////
//{{{ cSpinSystem
class cSpinSystem
{
public:
    cSpinSystem();
    cSpinSystem(const vector<cSPIN>& spin_list);
    ~cSpinSystem();

    void addSpinInteraction(cSpinInteraction& spin_interaction);
    SumKronProd getSumKronOperator();

protected:
private:
    DIM_LIST _dim_list;
    vector<cSPIN> _spin_list;
    vector<cSpinInteraction> _spin_interaction_list;
};
//}}}
////////////////////////////////////////////////////////////////////////////////
#endif
