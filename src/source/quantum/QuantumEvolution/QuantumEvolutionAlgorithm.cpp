#include "include/quantum/QuantumEvolutionAlgorithm.h"


////////////////////////////////////////////////////////////////////////////////
//{{{ QuantumEvolutionAlgorithm
QuantumEvolutionAlgorithm::QuantumEvolutionAlgorithm()
{ LOG(INFO) << "Default constructor: QuantumEvolutionAlgorithm";}

QuantumEvolutionAlgorithm::QuantumEvolutionAlgorithm(QuantumOperator& op, QuantumState& st)
{
    _operator_ptr = &op;
    _init_state_ptr = &st;
}

QuantumEvolutionAlgorithm::~QuantumEvolutionAlgorithm()
{ LOG(INFO) << "Default destructor: QuantumEvolutionAlgorithm";}
//}}}
////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
//{{{ SimpleFullMatrixVectorEvolution
SimpleFullMatrixVectorEvolution::SimpleFullMatrixVectorEvolution()
{ LOG(INFO) << "Default constructor: SimpleFullMatrixVectorEvolution";}

void SimpleFullMatrixVectorEvolution::perform()
{
    _matrix = _operator_ptr->getMatrix();
    _vector_list.push_back(_init_state_ptr->getVector());
    
    ////////////////////////////////////////////////////////////////////////////////
    //begin evolution
    double dt;
    double t_now     = _time_list[0];
    cx_vec state_now = _vector_list[0];
    cx_vec state_next;

    for(double t_next: _time_list)
    {
        cout << "Evolving form t=" << t_now << " to t=" << t_next << " ... " << endl;
        dt = t_next - t_now;
        state_next = expmat( -1.0*dt*II*_matrix ) * state_now;

        _vector_list.push_back( state_next );

        t_now = t_next;
        state_now = state_next;
    }
}

SimpleFullMatrixVectorEvolution::~SimpleFullMatrixVectorEvolution()
{ LOG(INFO) << "Default destructor: SimpleFullMatrixVectorEvolution";}
//}}}
////////////////////////////////////////////////////////////////////////////////
