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

    virtual void perform()=0;
protected:
    vec              _time_list;
    QuantumOperator* _operator_ptr;
    QuantumState*    _init_state_ptr;
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
    cx_mat         _matrix;
    vector<cx_vec> _vector_list;
};
//}}}
////////////////////////////////////////////////////////////////////////////////
#endif
