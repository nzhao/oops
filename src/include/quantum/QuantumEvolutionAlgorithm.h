#ifndef QUANTUMEVOLUTIONALGORITHM_H
#define QUANTUMEVOLUTIONALGORITHM_H
#include "include/misc/misc.h"
#include "include/quantum/QuantumOperator/QuantumOperator.h"
#include "include/quantum/QuantumState/QuantumState.h"

////////////////////////////////////////////////////////////////////////////////
//{{{ QuantumEvolutionAlgorithm
class QuantumEvolutionAlgorithm
{
public:
    QuantumEvolutionAlgorithm();
    ~QuantumEvolutionAlgorithm();
protected:
    bool _is_equidistant;
    vector<double> _time_list;
private:
};
//}}}
////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
//{{{ FullMatrixVectorEvolution
class FullMatrixVectorEvolution:public QuantumEvolutionAlgorithm
{
public:
    FullMatrixVectorEvolution();
    FullMatrixVectorEvolution(const QuantumOperator& op, const QuantumState& st);
    ~FullMatrixVectorEvolution();

    void perform();
protected:
private:
    cx_mat         _U0;
    cx_mat         _matrix;
    vector<cx_vec> _vector_list;

    cx_vec step(const cx_vec& vec_old, double dt);
    cx_vec step(const cx_vec& vec_old);
};
//}}}
////////////////////////////////////////////////////////////////////////////////
#endif
