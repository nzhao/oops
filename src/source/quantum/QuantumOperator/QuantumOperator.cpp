#include "include/quantum/QuantumOperator.h"

////////////////////////////////////////////////////////////////////////////////
//{{{ QuantumOperator
QuantumOperator::QuantumOperator()
{ LOG(INFO) << "Default constructor: QuantumOperator."; }

QuantumOperator::~QuantumOperator()
{ LOG(INFO) << "Default destructor: QuantumOperator."; }

void QuantumOperator::saveMatrix()
{
    cx_mat m= this->fullMatrix();
    mat m_r = real(m).t();
    mat m_i = -imag(m).t();

    const char *file = "../src/debug/Operator.mat";
    mxArray *pArray = mxCreateDoubleMatrix(_dimension,_dimension,mxCOMPLEX);

    int dim2=_dimension*_dimension;
    memcpy((void *)(mxGetPr(pArray)), (void *) m_r.memptr(), dim2*sizeof(double));
    memcpy((void *)(mxGetPi(pArray)), (void *) m_i.memptr(), dim2*sizeof(double));
    
    MATFile *mFile = matOpen(file, "w");
    matPutVariableAsGlobal(mFile, "OperatorMat", pArray);
    matClose(mFile);

    mxDestroyArray(pArray);
}
//}}}
////////////////////////////////////////////////////////////////////////////////
