#include "include/quantum/QuantumOperator.h"

////////////////////////////////////////////////////////////////////////////////
//{{{ QuantumOperator
QuantumOperator::QuantumOperator()
{ //LOG(INFO) << "Default constructor: QuantumOperator.";
}

QuantumOperator::~QuantumOperator()
{ //LOG(INFO) << "Default destructor: QuantumOperator.";
}

#ifdef HAS_MATLAB
void QuantumOperator::saveMatrix(string filename)
{
    cx_mat m= this->getMatrix();
    mat m_r = real(m).t();
    mat m_i = -imag(m).t();

    mxArray *pArray = mxCreateDoubleMatrix(_dimension,_dimension,mxCOMPLEX);

    int dim2=_dimension*_dimension;
    memcpy((void *)(mxGetPr(pArray)), (void *) m_r.memptr(), dim2*sizeof(double));
    memcpy((void *)(mxGetPi(pArray)), (void *) m_i.memptr(), dim2*sizeof(double));
    
    MATFile *mFile = matOpen(filename.c_str(), "w");
    matPutVariableAsGlobal(mFile, "OperatorMat", pArray);
    matClose(mFile);

    mxDestroyArray(pArray);
}
#else
void QuantumOperator::saveMatrix(string filename)
{
    cout << "MATLAB not installed." << endl;
}
#endif

QuantumOperator operator + (const QuantumOperator& op1, const QuantumOperator& op2)
{
    QuantumOperator res;
    res._dim_list = op1._dim_list;
    res._dimension = op1._dimension;
    SumKronProd skp1 =  op1.getKronProdForm();
    SumKronProd skp2 =  op2.getKronProdForm();
    res._kron_form = skp1 + skp2;
    return res;
}

QuantumOperator operator - (const QuantumOperator& op1, const QuantumOperator& op2)
{
    QuantumOperator temp = op2;
    QuantumOperator res = op1 + temp.scale(-1.0);
    return res;
}
///}}}
////////////////////////////////////////////////////////////////////////////////
