#ifndef SPININTERACTIONCOMPONENT_H
#define SPININTERACTIONCOMPONENT_H 

#include <vector>
#include <armadillo>
#include "include/spin/Spin.h"
#include "include/spin/SpinInteractionDefine.h"
#include "include/quantum/PureState.h"

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

    INDEX_LIST getIndexList() const {return _index_list;};
    vector< vector<cSPIN> > getSpinAggregate() const {return _spin_aggregate;};
    size_t getLength() const {return _index_list.size();};
    int get_nBody() const {return _nbody;};

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
    SingleSpin(const vector<cSPIN>& spin_list, const vector<int>& pick_up_spins);
    ~SingleSpin();
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
    size_t getLength(){return _mat_list.size();};
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
    TwoSpinInteractionForm(const cSpinInteractionDomain& domain);
    ~TwoSpinInteractionForm();
};
//}}}
//----------------------------------------------------------------------------//
//{{{ SingleSpinInteractionForm
class SingleSpinInteractionForm:public cSpinInteractionForm
{
public:
    SingleSpinInteractionForm(const cSpinInteractionDomain& domain);
    ~SingleSpinInteractionForm();
};
//}}}
//----------------------------------------------------------------------------//
//{{{ SingleSpinDephasing
class SingleSpinDephasing:public cSpinInteractionForm
{
public:
    SingleSpinDephasing(const cSpinInteractionDomain& domain, const vec& axis);
    ~SingleSpinDephasing(){};
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
    size_t getLength(){return _coeff_list.size();};
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
    DipolarInteractionCoeff(const cSpinInteractionDomain& domain);
    ~DipolarInteractionCoeff();
};
//}}}
//----------------------------------------------------------------------------//
//{{{ ZeemanInteractionCoeff
class ZeemanInteractionCoeff:public cSpinInteractionCoeff
{
public:
    ZeemanInteractionCoeff(const cSpinInteractionDomain& domain, const vec& magB);
    ~ZeemanInteractionCoeff();
};
//}}}
//----------------------------------------------------------------------------//
////////////////////////////////////////////////////////////////////////////////
//{{{ DipolarFieldInteractionCoeff
class DipolarFieldInteractionCoeff:public cSpinInteractionCoeff
{
public:
    DipolarFieldInteractionCoeff(const cSpinInteractionDomain& domain, const cSPIN& center_spin, const PureState& state);
    DipolarFieldInteractionCoeff(const cSpinInteractionDomain& domain, const vector<cSPIN>& spin_list, const vector<PureState>& state_list);
    DipolarFieldInteractionCoeff(const cSpinInteractionDomain& domain, const vector<cSPIN>& spin_list, const vector<PureState>& state_list, const vec& pre_factor_list);
    ~DipolarFieldInteractionCoeff();
protected:
private:
};
//}}}
//----------------------------------------------------------------------------//
//{{{ PolarizationCoeff
class PolarizationCoeff:public cSpinInteractionCoeff
{
public:
    PolarizationCoeff(const cSpinInteractionDomain& domain, const vector<vec>& pol);
    ~PolarizationCoeff();
};
//}}}
//----------------------------------------------------------------------------//
//{{{ SpinDephasingRate
class SpinDephasingRate:public cSpinInteractionCoeff
{
public:
    SpinDephasingRate(const cSpinInteractionDomain& domain, const double dephasing_rate);
    ~SpinDephasingRate() {};
};
//}}}
////////////////////////////////////////////////////////////////////////////////
///@}

///@}
///@}

#endif
