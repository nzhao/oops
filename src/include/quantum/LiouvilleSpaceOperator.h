#ifndef LIOUVILLESPACEOPERATOR_H
#define LIOUVILLESPACEOPERATOR_H

#include "include/quantum/QuantumOperator.h"
#include "include/quantum/HilbertSpaceOperator.h"

/// \addtogroup QuantumOperator 
/// @{

/// \defgroup LiouvilleSpaceOperator LiouvilleSpaceOperator
/// @{

struct EXPAN 
{
    HilbertSpaceOperator op;
    MatExpanFunc* func;
};


////////////////////////////////////////////////////////////////////////////////
//{{{ LiouvilleSpaceOperator
class LiouvilleSpaceOperator:public QuantumOperator
{
public:
    LiouvilleSpaceOperator();
    ~LiouvilleSpaceOperator();

    LiouvilleSpaceOperator(const vector<cSPIN>& spin_list);
    LiouvilleSpaceOperator(const HilbertSpaceOperator& op, MatExpanFunc* func);

    void addInteraction(cSpinInteraction& spin_interaction);
    void make();
protected:
    bool _is_expanded;
    EXPAN _expan;

    vector<cSPIN>            _spin_list;
    vector<cSpinInteraction> _interaction_list;
private:
};
//}}}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
//{{{ Liouvillian
class Liouvillian:public LiouvilleSpaceOperator
{
public:
    Liouvillian();
    Liouvillian(const QuantumOperator& op);
    Liouvillian(const Hamiltonian& hm, MatExpanFunc* func):LiouvilleSpaceOperator(hm, func){};
    Liouvillian(const Hamiltonian& hm):LiouvilleSpaceOperator(hm, CIRCLEC){};
    ~Liouvillian();
protected:
private:
};
//}}}
////////////////////////////////////////////////////////////////////////////////
/// @}
/// @}
#endif
