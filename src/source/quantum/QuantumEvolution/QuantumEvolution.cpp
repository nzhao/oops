#include "include/quantum/QuantumEvolution.h"

////////////////////////////////////////////////////////////////////////////////
//{{{ QuantumEvolution 
QuantumEvolution::QuantumEvolution()
{ //LOG(INFO) << "Default constructor: QuantumEvolution";
}

QuantumEvolution::~QuantumEvolution()
{ //LOG(INFO) << "Default destructor: QuantumEvolution";
}

//}}}
////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
//{{{ ClusterCoherenceEvolution
ClusterCoherenceEvolution::ClusterCoherenceEvolution()
{ //LOG(INFO) << "Default constructor: ClusterCoherenceEvolution";
}

ClusterCoherenceEvolution::~ClusterCoherenceEvolution()
{ //LOG(INFO) << "Default destructor: ClusterCoherenceEvolution";
}

vec ClusterCoherenceEvolution::calc_obs()
{
    _time_list = _kernel->getTimeSequence();
    size_t n_time = _time_list.size();
    
    vector<cx_vec>  state = _kernel->getResult();
    cx_mat stateMat( n_time, _kernel->getStateDim(), fill::zeros);
    for(int i=0; i<n_time; ++i)
        stateMat.row(i) = trans( state[i] );

    cx_vec init_state_vect =  _kernel->getInitalState();
    vec res = _kernel->getMatrixDim() * real(stateMat*init_state_vect);
    return res;
}

//}}}
////////////////////////////////////////////////////////////////////////////////
