#include "include/quantum/QuantumEvolutionAlgorithm.h"

////////////////////////////////////////////////////////////////////////////////
//{{{ QuantumEvolutionAlgorithm
QuantumEvolutionAlgorithm::QuantumEvolutionAlgorithm()
{ LOG(INFO) << "Default constructor: QuantumEvolutionAlgorithm";}

QuantumEvolutionAlgorithm::~QuantumEvolutionAlgorithm()
{ LOG(INFO) << "Default destructor: QuantumEvolutionAlgorithm";}
//}}}
////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
//{{{ MatrixVectorEvolution
FullMatrixVectorEvolution::MatrixVectorEvolution()
{ LOG(INFO) << "Default constructor: FullMatrixVectorEvolution";}

FullMatrixVectorEvolution::FullMatrixVectorEvolution(const QuantumOperator& op, const QuantumState& st)
{
    _matrix = op.getMatrix();

    _vector_list.reserve( _time_list.size() );
    _vector_list.push_back(st.getVector());

    if( _is_equidistant )
    {
        dt = _time_list[1] - _time_list[0];
        _U0 = expmat( -1.0 * dt * II * _matrix);
    }
}

cx_vec FullMatrixVectorEvolution::step(const cx_vec& vec_old, double dt)
{ return expmat(-1.0 *  dt * II * _matrix) * vec_old; }`

cx_vec FullMatrixVectorEvolution::step(const cx_vec& vec_old)
{ return _U0 * vec_old; }`

void FullMatrixVectorEvolution::perform()
{
    double t_now = _time_list[0];
    cx_vec state_now = _vector_list[0];
    cx_vec state_new;

    for(int i=1; i<_time_list.size(); ++i)
    {
        if( _is_equidistant)
            state_new = step(state_now);
        else
        {
            dt = _time_list[i] - t_now;
            state_new = step(state_now, dt);
        }
        _vector_list.push_back( state_new );

        t_now = _time_list[i];
        state_now = state_new;
    }
}

FullMatrixVectorEvolution::~FullMatrixVectorEvolution()
{ LOG(INFO) << "Default destructor: FullMatrixVectorEvolution";}
//}}}
////////////////////////////////////////////////////////////////////////////////
