#ifndef MIXEDSTATE_H
#define MIXEDSTATE_H
#include <armadillo>
#include "include/spin/Spin.h"
#include "include/spin/SpinInteraction.h"
#include "include/quantum/HilbertSpaceOperator.h"
#include "include/quantum/QuantumState.h"

/// \addtogroup QuantumState QuantumState
/// @{

/// \defgroup MixedState MixedState
/// @{

////////////////////////////////////////////////////////////////////////////////
//{{{ MixedState 
class MixedState:public QuantumState
{
public:
    MixedState();
    ~MixedState();

    virtual void make()=0;
    virtual void addStateComponent(cSpinInteraction& spin_interaction)=0;
protected:
private:
};
//}}}
////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
//{{{ DensityOperator 
class DensityOperator:public MixedState
{
public:
    DensityOperator();
    ~DensityOperator();

    DensityOperator(const vector<cSPIN>& spin_list);

    void make() {_op.make();};
    void makeVector(){ _vector = this->getKronProdForm().vecterize(); };
    //void makeVector(){ _vector = vectorise( this->getKronProdForm().full() ); };
    void addStateComponent(cSpinInteraction& spin_interaction){_op.addInteraction(spin_interaction);};

    SumKronProd getKronProdForm(){return _op.getKronProdForm();};
    cx_mat getMatrix(){return _op.getMatrix();};
protected:
private:
    HilbertSpaceOperator _op;
};
//}}}
////////////////////////////////////////////////////////////////////////////////

/// @}
/// @}
#endif
