#ifndef QUANTUMEVOLUTIONALGORITHM_H
#define QUANTUMEVOLUTIONALGORITHM_H
#include "include/misc/misc.h"
#include "include/quantum/QuantumOperator.h"
#include "include/quantum/QuantumState.h"

////////////////////////////////////////////////////////////////////////////////
//{{{ QuantumEvolutionAlgorithm
class QuantumEvolutionAlgorithm
{
public:
    QuantumEvolutionAlgorithm();
    QuantumEvolutionAlgorithm(QuantumOperator& op, QuantumState& st);
    ~QuantumEvolutionAlgorithm();

    void setTimeSequence(const vec& tlist){_time_list = tlist;};
    vec  getTimeSequence() const {return _time_list;};
    size_t  getStateDim() const {return _state_dimension;};
    vector<cx_vec> getResult() const {return _vector_list;};
    cx_vec getInitalState() const {return _init_state_ptr->getVector();}; 

    virtual void perform()=0;
protected:
    vec              _time_list;
    QuantumOperator* _operator_ptr;
    QuantumState*    _init_state_ptr;
    size_t           _state_dimension;

    cx_mat         _matrix;
    vector<cx_vec> _vector_list;
private:
};
//}}}
////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
//{{{ SimpleFullMatrixVectorEvolution
class SimpleFullMatrixVectorEvolution:public QuantumEvolutionAlgorithm
{
public:
    SimpleFullMatrixVectorEvolution();
    SimpleFullMatrixVectorEvolution(QuantumOperator& op, QuantumState& st):QuantumEvolutionAlgorithm(op, st){};
    ~SimpleFullMatrixVectorEvolution();


    void perform();
protected:
private:
};
//}}}
////////////////////////////////////////////////////////////////////////////////
#endif
