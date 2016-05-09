#ifndef QUANTUMEVOLUTION_H
#define QUANTUMEVOLUTION_H
#include "include/quantum/QuantumEvolutionAlgorithm.h"
#include <armadillo>

/// \addtogroup Quantum
/// @{
/// \defgroup QuantumEvolution QuantumEvolution
/// @{


////////////////////////////////////////////////////////////////////////////////
//{{{ QuantumEvllution 
class QuantumEvolution
{
public:
    QuantumEvolution() {};
    QuantumEvolution(QuantumEvolutionAlgorithm* kernel){ _kernel = kernel; };
    ~QuantumEvolution() {};

    void run() {_kernel->perform();};

protected:
    QuantumEvolutionAlgorithm* _kernel;
    vec _time_list;
private:
};
//}}}
////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
//{{{ ClusterCoherenceEvolution
class ClusterCoherenceEvolution:public QuantumEvolution
{
public:
    ClusterCoherenceEvolution() {};
    ClusterCoherenceEvolution(QuantumEvolutionAlgorithm* kernel):QuantumEvolution(kernel) {};
    ~ClusterCoherenceEvolution() {};

    vec calc_obs();
protected:
private:
};
//}}}
////////////////////////////////////////////////////////////////////////////////
/// @}
/// @}
#endif
