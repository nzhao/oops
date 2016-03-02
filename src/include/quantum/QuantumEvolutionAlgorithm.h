#ifndef QUANTUMEVOLUTIONALGORITHM_H
#define QUANTUMEVOLUTIONALGORITHM_H
#include "include/misc/misc.h"
#include "include/quantum/QuantumOperator.h"
#include "include/quantum/QuantumState.h"
#include "include/quantum/MixedState.h"
#include "include/quantum/PureState.h"
#include "include/math/MatExp.h"

////////////////////////////////////////////////////////////////////////////////
//{{{ QuantumEvolutionAlgorithm
class QuantumEvolutionAlgorithm
{
public:
    QuantumEvolutionAlgorithm() {};
    QuantumEvolutionAlgorithm(QuantumOperator& op, QuantumState& st);
    ~QuantumEvolutionAlgorithm() {};

    void setTimeSequence(double t0, double t1, int nt) { _time_list = linspace<vec>(t0, t1, nt); };
    //void setTimeSequence(const vec& tlist){_time_list = tlist;};
    vec  getTimeSequence() const {return _time_list;};
    size_t  getStateDim() const {return _state_dimension;};
    vector<cx_vec> getResult() const {return _vector_list;};
    vector<cx_mat> getResultMat() const {return _state_mat_list;};
    cx_vec getInitalState() const {return _init_state.getVector();}; 
    size_t getMatrixDim() const {return _init_state.getDimension();}

    virtual void perform()=0;
protected:
    vec              _time_list;
    QuantumOperator  _operator;
    QuantumState     _init_state;
    size_t           _state_dimension;

    cx_mat         _matrix;
    vector<cx_vec> _vector_list;
    vector<cx_mat> _state_mat_list;
private:
};
//}}}
////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
//{{{ SimpleFullMatrixVectorEvolution
class SimpleFullMatrixVectorEvolution:public QuantumEvolutionAlgorithm
{
public:
    SimpleFullMatrixVectorEvolution() {};
    SimpleFullMatrixVectorEvolution(QuantumOperator& op, QuantumState& st):QuantumEvolutionAlgorithm(op, st){};
    ~SimpleFullMatrixVectorEvolution() {};

    void perform();
protected:
private:
};
//}}}
////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
//{{{ PiecewiseFullMatrixVectorEvolution
class PiecewiseFullMatrixVectorEvolution:public QuantumEvolutionAlgorithm
{
public:
    PiecewiseFullMatrixVectorEvolution() {};
    PiecewiseFullMatrixVectorEvolution(const vector<QuantumOperator>& op_list, const vector<double>& time_segment, const QuantumState& st);
    ~PiecewiseFullMatrixVectorEvolution() {};

    void perform();
protected:
private:
    vector<QuantumOperator> _op_list;
    vector<double> _time_segment;
};
//}}}
////////////////////////////////////////////////////////////////////////////////




////////////////////////////////////////////////////////////////////////////////
//{{{  PiecewiseFullMatrixMatrixEvolution
class PiecewiseFullMatrixMatrixEvolution:public QuantumEvolutionAlgorithm
{
public:
    PiecewiseFullMatrixMatrixEvolution() {};
    PiecewiseFullMatrixMatrixEvolution(
            const vector<QuantumOperator>& left_op_list, 
            const vector<QuantumOperator>& right_op_list, const vector<double>& time_segment, const DensityOperator& ds);
    ~PiecewiseFullMatrixMatrixEvolution() {};

    void perform();
protected:
private:
    DensityOperator _density_matrix;
    vector<QuantumOperator> _left_op_list;
    vector<QuantumOperator> _right_op_list;
    vector<double> _time_segment;
};

//}}}
////////////////////////////////////////////////////////////////////////////////
#endif
