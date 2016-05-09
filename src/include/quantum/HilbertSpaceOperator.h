#ifndef HILBERTSPACEOPERATOR_H
#define HILBERTSPACEOPERATOR_H
#include <armadillo>
#include "include/spin/Spin.h"
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
    HilbertSpaceOperator(const vector<cSPIN>& spin_list);
    ~HilbertSpaceOperator();

    void addInteraction(cSpinInteraction& spin_interaction);
    void make();
protected:
    vector<cSPIN> _spin_list;
    vector<cSpinInteraction> _interaction_list;
private:
};
//}}}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
//{{{ Hamiltonian
class Hamiltonian:public HilbertSpaceOperator
{
public:
    Hamiltonian();
    Hamiltonian(const vector<cSPIN>& spin_list):HilbertSpaceOperator(spin_list){};
    ~Hamiltonian();


private:
};
//}}}
////////////////////////////////////////////////////////////////////////////////

/// @}
/// @}
#endif
