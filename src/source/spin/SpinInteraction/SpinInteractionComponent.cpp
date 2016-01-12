#include "include/misc/misc.h"
#include "include/spin/SpinInteractionComponent.h"

using namespace std;

////////////////////////////////////////////////////////////////////////////////
//{{{ cSpinInteractionDomain
cSpinInteractionDomain::cSpinInteractionDomain()
{ LOG(INFO) << "Default constructor: cSpinInteractionDomain.";
}
cSpinInteractionDomain::~cSpinInteractionDomain()
{ LOG(INFO) << "Default destructor: cSpinInteractionDomain.";
}

ostream&  operator << (ostream& outs, const cSpinInteractionDomain& dm)
{
    int i=0;
    for(auto idx: dm._index_list)
    {
        outs << "interaction domain[" << i << "]: ";
        for(int j=0; j<dm._nbody; ++j)
        {
            outs << idx[j] << ", ";
        }
        outs << endl;
        i++;
    }
    return outs;
}
//}}}
//----------------------------------------------------------------------------//
//{{{ SpinPair
SpinPair::SpinPair(const vector<cSPIN>& spin_list)
{ LOG(INFO) << "Constructor: SpinPair with spin_list";
    int nspin=spin_list.size();

    _nbody = 2;

    for(int i=0; i<nspin; ++i)
        for(int j=i+1; j<nspin; ++j)
            _index_list.push_back(vector<int> {i, j});

    for(auto idx:_index_list)
        _spin_aggregate.push_back( vector<cSPIN> { spin_list[ idx[0] ] , spin_list[ idx[1] ] });
}
SpinPair::~SpinPair()
{ LOG(INFO) << "Default destructor: SpinPair.";
}
//}}}
//----------------------------------------------------------------------------//
//{{{ SingleSpin
SingleSpin::SingleSpin(const vector<cSPIN>& spin_list)
{ LOG(INFO) << "Constructor: SingleSpin with spin_list";
    int nspin=spin_list.size();

    _nbody = 1;

    for(int i=0; i<nspin; ++i)
        _index_list.push_back( vector<int> {i} );

    for(auto idx:_index_list)
        _spin_aggregate.push_back( vector<cSPIN> { spin_list[ idx[0] ] });
}
SingleSpin::~SingleSpin()
{ LOG(INFO) << "Default destructor: SingleSpin.";
}
//}}}
//----------------------------------------------------------------------------//
//{{{ IndividualSpin 
IndividualSpin::IndividualSpin(const vector<cSPIN>& spin_list, int i)
{ LOG(INFO) << "Default constructor: IndividualSpin";

    int nspin=spin_list.size();
    _nbody = 1;
    _index_list.push_back( vector<int> {i} );
    _spin_aggregate.push_back( vector<cSPIN> { spin_list[i] });
}

IndividualSpin::~IndividualSpin()
{ LOG(INFO) << "Default destructor: IndividualSpin";}
//}}}
////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
//{{{ cSpinInteractionForm
cSpinInteractionForm::cSpinInteractionForm()
{  LOG(INFO) << "Default constructor: cSpinInteractionForm.";
}
cSpinInteractionForm::~cSpinInteractionForm()
{ LOG(INFO) << "Default destructor: cSpinInteractionForm.";
}
ostream&  operator << (ostream& outs, cSpinInteractionForm& form)
{
    int i = 0;
    MAT_LIST matlist=form.getMatList();
    for(auto term_list : matlist)
    {
        int j = 0;
        for(auto term : term_list)
        {
            outs << "form[" << i << "], " << "term[" << j << "]= " << endl << endl;
            int q = 0;
            for(auto mat : term)
            {
                outs << mat << endl;
                if(q<term.size()-1)
                    outs << " \t \t \t * " << endl << endl;
                q++;
            }
            outs << "--------------------------------------------------------"<< endl;
            outs << endl;
            j++;
        }
        i++;
    }
    return outs;
}
//}}}
//----------------------------------------------------------------------------//
//{{{ TwoSpinInteractionFrom
TwoSpinInteractionForm::TwoSpinInteractionForm(cSpinInteractionDomain& domain)
{
    _nterm = 9;

    auto sag = domain.getSpinAggregate();
    for(auto it=sag.begin(); it!=sag.end(); ++it)
    {
        cSPIN spin0=(*it)[0];    cSPIN spin1=(*it)[1];

        vector<TERM> term_list;

        term_list.push_back( TERM { spin0.sx(), spin1.sx() } );
        term_list.push_back( TERM { spin0.sx(), spin1.sy() } );
        term_list.push_back( TERM { spin0.sx(), spin1.sz() } );

        term_list.push_back( TERM { spin0.sy(), spin1.sx() } );
        term_list.push_back( TERM { spin0.sy(), spin1.sy() } );
        term_list.push_back( TERM { spin0.sy(), spin1.sz() } );

        term_list.push_back( TERM { spin0.sz(), spin1.sx() } );
        term_list.push_back( TERM { spin0.sz(), spin1.sy() } );
        term_list.push_back( TERM { spin0.sz(), spin1.sz() } );

        _mat_list.push_back( term_list );
    }
}
TwoSpinInteractionForm::~TwoSpinInteractionForm()
{ LOG(INFO) << "Default destructor: TwoSpinInteractionForm.";
}
//}}}
//----------------------------------------------------------------------------//
//{{{ SingleSpinInteractionForm
SingleSpinInteractionForm::SingleSpinInteractionForm(cSpinInteractionDomain& domain)
{
    _nterm = 6;

    auto sag = domain.getSpinAggregate();
    for(auto it=sag.begin(); it!=sag.end(); ++it)
    {
        cSPIN spin0=(*it)[0];

        vector<TERM> term_list;

        term_list.push_back( TERM { spin0.sx() } );
        term_list.push_back( TERM { spin0.sy() } );
        term_list.push_back( TERM { spin0.sz() } );

        term_list.push_back( TERM { spin0.sx()*spin0.sx() } );
        term_list.push_back( TERM { spin0.sy()*spin0.sy() } );
        term_list.push_back( TERM { spin0.sz()*spin0.sz() } );

        _mat_list.push_back( term_list );
    }
}
SingleSpinInteractionForm::~SingleSpinInteractionForm()
{ LOG(INFO) << "Default destructor: SingleSpinInteractionForm.";
}
//}}}
//----------------------------------------------------------------------------//
//{{{ IndividualSpinForm
IndividualSpinForm::IndividualSpinForm(cSpinInteractionDomain& domain)
{
    _nterm = 3;

    auto sag = domain.getSpinAggregate();
    for(auto it=sag.begin(); it!=sag.end(); ++it)
    {
        cSPIN spin0=(*it)[0];

        vector<TERM> term_list;

        term_list.push_back( TERM { spin0.sx() } );
        term_list.push_back( TERM { spin0.sy() } );
        term_list.push_back( TERM { spin0.sz() } );

        _mat_list.push_back( term_list );
    }
}
IndividualSpinForm::~IndividualSpinForm()
{ LOG(INFO) << "Default destructor: IndividualSpinForm.";
}
//}}}
////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
//{{{ cSpinInteractionCoeff
cSpinInteractionCoeff::cSpinInteractionCoeff()
{ LOG(INFO) << "Default constructor: cSpinInteractionCoeff.";
}
cSpinInteractionCoeff::~cSpinInteractionCoeff()
{ LOG(INFO) << "Default destructor: cSpinInteractionCoeff.";
}
ostream&  operator << (ostream& outs, cSpinInteractionCoeff& coef)
{
    int i = 0;
    COEFF_LIST coefflist = coef.getCoeffList(); 
    for(auto co : coefflist)
    {
        int j = 0;
        for(auto val : co)
        {
            cout << val << ", ";
        }
        cout << endl;
    }
   return outs;
}
//}}}
//----------------------------------------------------------------------------//
//{{{ DipolarInteractionCoeff
DipolarInteractionCoeff::DipolarInteractionCoeff(cSpinInteractionDomain& domain)
{
    _nCoeff = 9;

    auto sag = domain.getSpinAggregate();
    for(auto it=sag.begin(); it!=sag.end(); ++it)
    {
        cSPIN spin0=(*it)[0];    cSPIN spin1=(*it)[1];
        vec coeffs = dipole(spin0, spin1);
        _coeff_list.push_back(coeffs);
    }
}
DipolarInteractionCoeff::~DipolarInteractionCoeff()
{ LOG(INFO) << "Default destructor: DipolarInteractionCoeff.";
}
//}}}
//----------------------------------------------------------------------------//
//{{{ ZeemanInteractionCoeff
ZeemanInteractionCoeff::ZeemanInteractionCoeff(cSpinInteractionDomain& domain, const vec& magB)
{
    _nCoeff = 6;

    auto sag = domain.getSpinAggregate();
    for(auto it=sag.begin(); it!=sag.end(); ++it)
    {
        cSPIN spin0=(*it)[0];
        vec coeffs = zeeman(spin0, magB);
        _coeff_list.push_back(coeffs);
    }
}
ZeemanInteractionCoeff::~ZeemanInteractionCoeff()
{ LOG(INFO) << "Default destructor: ZeemanInteractionCoeff.";
}
//}}}
//----------------------------------------------------------------------------//
//{{{ PolarizationCoeff
PolarizationCoeff::PolarizationCoeff(cSpinInteractionDomain& domain, const vec& pol)
{
    _nCoeff = 3;
    _coeff_list.push_back(pol);
}
PolarizationCoeff::~PolarizationCoeff()
{ LOG(INFO) << "Default destructor: PolarizationCoeff.";
}
//}}}
////////////////////////////////////////////////////////////////////////////////

