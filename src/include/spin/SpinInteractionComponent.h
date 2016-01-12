#ifndef SPININTERACTIONCOMPONENT_H
#define SPININTERACTIONCOMPONENT_H 

#include <vector>
#include <armadillo>
#include "include/easylogging++.h"
#include "include/spin/Spin.h"
#include "include/Spin/SpinInteractionDefine.h"

using namespace std;
using namespace arma;

/// \addtogroup SpinInteraction 
/// @{

/// \defgroup SpinInteractionComponent SpinInteractionComponent
/// @{

/// \defgroup SpinInteractoinDomain 
/// @{
////////////////////////////////////////////////////////////////////////////////
//{{{ cSpinInteractionDomain
class cSpinInteractionDomain
{
public:
     cSpinInteractionDomain();
    ~cSpinInteractionDomain();

    INDEX_LIST getIndexList(){return _index_list;};
    vector< vector<cSPIN> > getSpinAggregate(){return _spin_aggregate;};
    int getLength(){return _index_list.size();};
    int get_nBody(){return _nbody;};

    friend ostream&  operator << (ostream& outs, const cSpinInteractionDomain& dm);
protected:
    int _nbody;
    vector< vector<cSPIN> > _spin_aggregate;
    INDEX_LIST _index_list;
};
//}}}
//----------------------------------------------------------------------------//
//{{{ SpinPair
class SpinPair:public cSpinInteractionDomain
{
public:
    SpinPair(const vector<cSPIN>& spin_list);
    ~SpinPair();
};
//}}}
//----------------------------------------------------------------------------//
//{{{ SingleSpin
class SingleSpin:public cSpinInteractionDomain
{
public:
    SingleSpin(const vector<cSPIN>& spin_list);
    ~SingleSpin();
};
//}}}
//----------------------------------------------------------------------------//
//{{{ IndividualSpin
class IndividualSpin:public cSpinInteractionDomain
{
public:
    IndividualSpin(const vector<cSPIN>& spin_list, int i);
    ~IndividualSpin();
protected:
private:
};
//}}}
////////////////////////////////////////////////////////////////////////////////
/// @}

/// \defgroup SpinInteractionForm
/// @{
////////////////////////////////////////////////////////////////////////////////
//{{{ cSpinInteractionForm
class cSpinInteractionForm
{
public:
     cSpinInteractionForm();
    ~cSpinInteractionForm();

    MAT_LIST getMatList(){return _mat_list;};
    int getLength(){return _mat_list.size();};
    int get_nTerm(){return _nterm;};

    friend ostream&  operator << (ostream& outs, cSpinInteractionForm& form);
protected:
    int _nterm;
    MAT_LIST _mat_list;
};
//}}}
//----------------------------------------------------------------------------//
//{{{ TwoSpinInteractionForm
class TwoSpinInteractionForm:public cSpinInteractionForm
{
public:
    TwoSpinInteractionForm(cSpinInteractionDomain& domain);
    ~TwoSpinInteractionForm();
};
//}}}
//----------------------------------------------------------------------------//
//{{{ SingleSpinInteractionForm
class SingleSpinInteractionForm:public cSpinInteractionForm
{
public:
    SingleSpinInteractionForm(cSpinInteractionDomain& domain);
    ~SingleSpinInteractionForm();
};
//}}}
//----------------------------------------------------------------------------//
//{{{ IndividualSpinForm
class IndividualSpinForm:public cSpinInteractionForm
{
public:
    IndividualSpinForm(cSpinInteractionDomain& domain);
    ~IndividualSpinForm();
};
//}}}
////////////////////////////////////////////////////////////////////////////////
/// @}


/// \defgroup SpinInteractionCoeff 
/// @{
////////////////////////////////////////////////////////////////////////////////
//{{{ cSpinInteractionCoeff
class cSpinInteractionCoeff
{
public:
     cSpinInteractionCoeff();
    ~cSpinInteractionCoeff();

    COEFF_LIST getCoeffList(){return _coeff_list;};
    int getLength(){return _coeff_list.size();};
    int get_nCoeff(){return _nCoeff;};

    friend ostream&  operator << (ostream& outs, cSpinInteractionCoeff& coef);
protected:
    int _nCoeff;
    COEFF_LIST _coeff_list;
};
//}}}
//----------------------------------------------------------------------------//
//{{{ DipolarInteractionCoeff
class DipolarInteractionCoeff:public cSpinInteractionCoeff
{
public:
    DipolarInteractionCoeff(cSpinInteractionDomain& domain);
    ~DipolarInteractionCoeff();
};
//}}}
//----------------------------------------------------------------------------//
//{{{ ZeemanInteractionCoeff
class ZeemanInteractionCoeff:public cSpinInteractionCoeff
{
public:
    ZeemanInteractionCoeff(cSpinInteractionDomain& domain, const vec& magB);
    ~ZeemanInteractionCoeff();
};
//}}}
//----------------------------------------------------------------------------//
//{{{ PolarizationCoeff
class PolarizationCoeff:public cSpinInteractionCoeff
{
public:
    PolarizationCoeff(cSpinInteractionDomain& domain, const vec& pol);
    ~PolarizationCoeff();
};
//}}}
////////////////////////////////////////////////////////////////////////////////
///@}

///@}
///@}

#endif
