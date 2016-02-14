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
void QuantumOperator::saveMatrix(string name)
{
    cx_mat m= this->getMatrix();
    mat m_r = real(m).t();
    mat m_i = -imag(m).t();

    mxArray *pArray = mxCreateDoubleMatrix(_dimension,_dimension,mxCOMPLEX);

    int dim2=_dimension*_dimension;
    memcpy((void *)(mxGetPr(pArray)), (void *) m_r.memptr(), dim2*sizeof(double));
    memcpy((void *)(mxGetPi(pArray)), (void *) m_i.memptr(), dim2*sizeof(double));
    
    char dbg_filename[500];
    strcpy(dbg_filename, PROJECT_PATH);
    strcat(dbg_filename, "/dat/debug/");
    strcat(dbg_filename, name.c_str());
    strcat(dbg_filename, ".mat");
    cout << dbg_filename << endl;
    MATFile *mFile = matOpen(dbg_filename, "w");
    matPutVariableAsGlobal(mFile, name.c_str(), pArray);
    matClose(mFile);

    mxDestroyArray(pArray);
}
#else
void QuantumOperator::saveMatrix(string name)
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
