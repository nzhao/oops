#ifndef HILBERTSPACEOPERATOR_H
#define HILBERTSPACEOPERATOR_H
#include <armadillo>
#include "include/spin/spin.h"
#include "include/easylogging++.h"
#include "include/quantum/QuantumOperator.h"
#include "include/spin/SpinInteraction.h"
/// \addtogroup QuantumOperator
/// @{
//
/// \defgroup HilbertSpaceOperator HilbertSpaceOperator
/// @{


////////////////////////////////////////////////////////////////////////////////
//{{{ HilbertSpaceOperator
class HilbertSpaceOperator:public QuantumOperator
{
public:
    HilbertSpaceOperator();
    ~HilbertSpaceOperator();

protected:
private:
};
#endif
//}}}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
//{{{ Hamiltonian
class Hamiltonian:public HilbertSpaceOperator
{
public:
    Hamiltonian();
    Hamiltonian(const vector<cSPIN>& spin_list);
    ~Hamiltonian();

    void addInteraction(cSpinInteraction& spin_interaction);
    void makeKronForm();

private:
    vector<cSPIN> _spin_list;
    vector<cSpinInteraction> _interaction_list;
};
//}}}
////////////////////////////////////////////////////////////////////////////////

/// @}
/// @}
