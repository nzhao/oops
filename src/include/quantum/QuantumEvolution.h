#ifndef QUANTUMEVOLUTION_H
#define QUANTUMEVOLUTION_H
#include "include/easylogging++.h"
#include "include/quantum/QuantumEvolutionAlgorithm.h"


/// \addtogroup Quantum
/// @{
/// \defgroup QuantumEvolution QuantumEvolution
/// @{


////////////////////////////////////////////////////////////////////////////////
//{{{ QuantumEvllution 
class QuantumEvolution
{
public:
    QuantumEvolution();
    QuantumEvolution(QuantumEvolutionAlgorithm* kernel){ _kernel = kernel; };
    ~QuantumEvolution();

    void run() {_kernel->perform();};
protected:
private:
    QuantumEvolutionAlgorithm* _kernel;
};
//}}}
////////////////////////////////////////////////////////////////////////////////
/// @}
/// @}
#endif
