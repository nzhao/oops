#include "include/misc/misc.h"
#include "include/spin/SpinInteractionComponent.h"

using namespace std;

////////////////////////////////////////////////////////////////////////////////
//{{{ cSpinInteractionDomain
cSpinInteractionDomain::cSpinInteractionDomain()
{ //LOG(INFO) << "Default constructor: cSpinInteractionDomain.";
}
cSpinInteractionDomain::~cSpinInteractionDomain()
{ //LOG(INFO) << "Default destructor: cSpinInteractionDomain.";
}

ostream&  operator << (ostream& outs, const cSpinInteractionDomain& dm)
{
    int i=0;
    //for(auto idx: dm._index_list)
    //{
    //    outs << "interaction domain[" << i << "]: ";
    //    for(int j=0; j<dm._nbody; ++j)
    //    {
    //        outs << idx[j] << ", ";
    //    }
    //    outs << endl;
    //    i++;
    //}
    for(int q=0; q<dm._index_list.size(); ++q)
    {
        outs << "interaction domain[" << i << "]: ";
        for(int j=0; j<dm._nbody; ++j)
        {
            outs << dm._index_list[q][j] << ", ";
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
{ //LOG(INFO) << "Constructor: SpinPair with spin_list";
    size_t nspin=spin_list.size();

    _nbody = 2;

    for(int i=0; i<nspin; ++i)
        for(int j=i+1; j<nspin; ++j)
        {
            vector<int> x; x.push_back(i); x.push_back(j);
            _index_list.push_back(x);
        }

    //for(auto idx:_index_list)
    //    _spin_aggregate.push_back( vector<cSPIN> { spin_list[ idx[0] ] , spin_list[ idx[1] ] });
    for(int i=0; i<_index_list.size(); ++i)
    {
        vector<cSPIN> x;
        x.push_back( spin_list[ _index_list[i][0] ] );
        x.push_back( spin_list[ _index_list[i][1] ] );
        _spin_aggregate.push_back(x);
    }
}
SpinPair::~SpinPair()
{ //LOG(INFO) << "Default destructor: SpinPair.";
}
//}}}
//----------------------------------------------------------------------------//
//{{{ SingleSpin
SingleSpin::SingleSpin(const vector<cSPIN>& spin_list)
{ //LOG(INFO) << "Constructor: SingleSpin with spin_list";
    size_t nspin=spin_list.size();

    _nbody = 1;

    //for(int i=0; i<nspin; ++i)
    //    _index_list.push_back( vector<int> {i} );
    for(int i=0; i<nspin; ++i)
    {
        vector<int> x; x.push_back(i);
        _index_list.push_back(x);
    }

    //for(auto idx:_index_list)
    //    _spin_aggregate.push_back( vector<cSPIN> { spin_list[ idx[0] ] });
    for(int i=0; i<_index_list.size(); ++i)
    {
        vector<cSPIN> x;
        x.push_back( spin_list[ _index_list[i][0] ] );
        _spin_aggregate.push_back(x);
    }
}

SingleSpin::SingleSpin(const vector<cSPIN>& spin_list, const vector<int>& pick_up_spins)
{ //LOG(INFO) << "Constructor: SingleSpin with spin_list and pick_up list.";
    //size_t nspin=spin_list.size();
    _nbody = 1;

    //for(int i : pick_up_spins)
    //    _index_list.push_back( vector<int> {i} );
    for(int i=0; i<pick_up_spins.size(); ++i)
    {
        vector<int> x; x.push_back(i);
        _index_list.push_back(x);
    }

    //for(auto idx:_index_list)
    //    _spin_aggregate.push_back( vector<cSPIN> { spin_list[ idx[0] ] });
    for(int i=0; i<_index_list.size(); ++i)
    {
        vector<cSPIN> x;
        x.push_back( spin_list[ _index_list[i][0] ] );
        _spin_aggregate.push_back(x);
    }
}

SingleSpin::~SingleSpin()
{ //LOG(INFO) << "Default destructor: SingleSpin.";
}
//}}}
////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
//{{{ cSpinInteractionForm
cSpinInteractionForm::cSpinInteractionForm()
{  //LOG(INFO) << "Default constructor: cSpinInteractionForm.";
}
cSpinInteractionForm::~cSpinInteractionForm()
{ //LOG(INFO) << "Default destructor: cSpinInteractionForm.";
}
ostream&  operator << (ostream& outs, cSpinInteractionForm& form)
{
    int i = 0;
    MAT_LIST matlist=form.getMatList();
    //for(auto term_list : matlist)
    //{
    //    int j = 0;
    //    for(auto term : term_list)
    //    {
    //        outs << "form[" << i << "], " << "term[" << j << "]= " << endl << endl;
    //        int q = 0;
    //        for(auto mat : term)
    //        {
    //            outs << mat << endl;
    //            if(q<term.size()-1)
    //                outs << " \t \t \t * " << endl << endl;
    //            q++;
    //        }
    //        outs << "--------------------------------------------------------"<< endl;
    //        outs << endl;
    //        j++;
    //    }
    //    i++;
    //}
    for(int kk=0; kk<matlist.size(); ++kk)
    {
        int j = 0;
        for(int qq=0; qq<matlist[kk].size(); ++qq)
        {
            outs << "form[" << i << "], " << "term[" << j << "]= " << endl << endl;
            int q = 0;
            //for(auto mat : term)
            for(int zz=0; zz<matlist[kk][qq].size(); ++zz)
            {
                outs << matlist[kk][qq][zz] << endl;
                if(q<matlist[kk][qq].size()-1)
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
TwoSpinInteractionForm::TwoSpinInteractionForm(const cSpinInteractionDomain& domain)
{
    _nterm = 9;

    //auto sag = domain.getSpinAggregate();
    //for(auto it=sag.begin(); it!=sag.end(); ++it)
    //{
    //    cSPIN spin0=(*it)[0];    cSPIN spin1=(*it)[1];

    //    vector<TERM> term_list;

    //    term_list.push_back( TERM { spin0.sx(), spin1.sx() } );
    //    term_list.push_back( TERM { spin0.sx(), spin1.sy() } );
    //    term_list.push_back( TERM { spin0.sx(), spin1.sz() } );

    //    term_list.push_back( TERM { spin0.sy(), spin1.sx() } );
    //    term_list.push_back( TERM { spin0.sy(), spin1.sy() } );
    //    term_list.push_back( TERM { spin0.sy(), spin1.sz() } );

    //    term_list.push_back( TERM { spin0.sz(), spin1.sx() } );
    //    term_list.push_back( TERM { spin0.sz(), spin1.sy() } );
    //    term_list.push_back( TERM { spin0.sz(), spin1.sz() } );

    //    _mat_list.push_back( term_list );
    //}
    vector< vector<cSPIN> > sag;
    sag = domain.getSpinAggregate();
    vector< vector<cSPIN> >::iterator it;
    for(it=sag.begin(); it!=sag.end(); ++it)
    {
        cSPIN spin0=(*it)[0];    cSPIN spin1=(*it)[1];

        vector<TERM> term_list;

        TERM t;  t.reserve(2);
        t.push_back( spin0.sx() ); t.push_back( spin1.sx() ); term_list.push_back( t ); t.clear();
        t.push_back( spin0.sx() ); t.push_back( spin1.sy() ); term_list.push_back( t ); t.clear();
        t.push_back( spin0.sx() ); t.push_back( spin1.sz() ); term_list.push_back( t ); t.clear();
        t.push_back( spin0.sy() ); t.push_back( spin1.sx() ); term_list.push_back( t ); t.clear();
        t.push_back( spin0.sy() ); t.push_back( spin1.sy() ); term_list.push_back( t ); t.clear();
        t.push_back( spin0.sy() ); t.push_back( spin1.sz() ); term_list.push_back( t ); t.clear();
        t.push_back( spin0.sz() ); t.push_back( spin1.sx() ); term_list.push_back( t ); t.clear();
        t.push_back( spin0.sz() ); t.push_back( spin1.sy() ); term_list.push_back( t ); t.clear();
        t.push_back( spin0.sz() ); t.push_back( spin1.sz() ); term_list.push_back( t ); t.clear();

        _mat_list.push_back( term_list );
    }
}
TwoSpinInteractionForm::~TwoSpinInteractionForm()
{ //LOG(INFO) << "Default destructor: TwoSpinInteractionForm.";
}
//}}}
//----------------------------------------------------------------------------//
//{{{ SingleSpinInteractionForm
SingleSpinInteractionForm::SingleSpinInteractionForm(const cSpinInteractionDomain& domain)
{
    _nterm = 6;

    //auto sag = domain.getSpinAggregate();
    //for(auto it=sag.begin(); it!=sag.end(); ++it)
    //{
    //    cSPIN spin0=(*it)[0];

    //    vector<TERM> term_list;

    //    term_list.push_back( TERM { spin0.sx() } );
    //    term_list.push_back( TERM { spin0.sy() } );
    //    term_list.push_back( TERM { spin0.sz() } );

    //    term_list.push_back( TERM { spin0.sx()*spin0.sx() } );
    //    term_list.push_back( TERM { spin0.sy()*spin0.sy() } );
    //    term_list.push_back( TERM { spin0.sz()*spin0.sz() } );

    //    _mat_list.push_back( term_list );
    //}
    vector< vector<cSPIN> > sag;
    sag = domain.getSpinAggregate();
    vector< vector<cSPIN> >::iterator it;
    for(it=sag.begin(); it!=sag.end(); ++it)
    {
        cSPIN spin0=(*it)[0];

        vector<TERM> term_list;

        TERM t;  t.reserve(1);
        t.push_back( spin0.sx() ); term_list.push_back( t ); t.clear();
        t.push_back( spin0.sy() ); term_list.push_back( t ); t.clear();
        t.push_back( spin0.sz() ); term_list.push_back( t ); t.clear();
        t.push_back( spin0.sx()*spin0.sx() ); term_list.push_back( t ); t.clear();
        t.push_back( spin0.sy()*spin0.sy() ); term_list.push_back( t ); t.clear();
        t.push_back( spin0.sz()*spin0.sz() ); term_list.push_back( t ); t.clear();

        _mat_list.push_back( term_list );
    }
}
SingleSpinInteractionForm::~SingleSpinInteractionForm()
{ //LOG(INFO) << "Default destructor: SingleSpinInteractionForm.";
}
//}}}
////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
//{{{ cSpinInteractionCoeff
cSpinInteractionCoeff::cSpinInteractionCoeff()
{ //LOG(INFO) << "Default constructor: cSpinInteractionCoeff.";
}
cSpinInteractionCoeff::~cSpinInteractionCoeff()
{ //LOG(INFO) << "Default destructor: cSpinInteractionCoeff.";
}
ostream&  operator << (ostream& outs, cSpinInteractionCoeff& coef)
{
    //int i = 0;
    COEFF_LIST coefflist = coef.getCoeffList(); 
    //for(auto co : coefflist)
    //{
    //    int j = 0;
    //    for(auto val : co)
    //    {
    //        cout << val << ", ";
    //    }
    //    cout << endl;
    //}
    for(int kk=0; kk<coefflist.size(); ++kk)
    {
        //int j = 0;
        for(int qq=0; qq<coefflist[kk].size(); ++qq)
        {
            cout << coefflist[kk][qq] << ", ";
        }
        cout << endl;
    }
   return outs;
}
//}}}
//----------------------------------------------------------------------------//
//{{{ DipolarInteractionCoeff
DipolarInteractionCoeff::DipolarInteractionCoeff(const cSpinInteractionDomain& domain)
{
    _nCoeff = 9;

    //auto sag = domain.getSpinAggregate();
    //for(auto it=sag.begin(); it!=sag.end(); ++it)
    //{
    //    cSPIN spin0=(*it)[0];    cSPIN spin1=(*it)[1];
    //    vec coeffs = dipole(spin0, spin1);
    //    _coeff_list.push_back(coeffs);
    //}
    vector< vector<cSPIN> > sag;
    sag = domain.getSpinAggregate();
    vector< vector<cSPIN> >::iterator it;
    for(it=sag.begin(); it!=sag.end(); ++it)
    {
        cSPIN spin0=(*it)[0];    cSPIN spin1=(*it)[1];
        vec coeffs = dipole(spin0, spin1);
        _coeff_list.push_back(coeffs);
    }
}
DipolarInteractionCoeff::~DipolarInteractionCoeff()
{ //LOG(INFO) << "Default destructor: DipolarInteractionCoeff.";
}
//}}}
//----------------------------------------------------------------------------//
//{{{ ZeemanInteractionCoeff
ZeemanInteractionCoeff::ZeemanInteractionCoeff(const cSpinInteractionDomain& domain, const vec& magB)
{
    _nCoeff = 6;

    //auto sag = domain.getSpinAggregate();
    //for(auto it=sag.begin(); it!=sag.end(); ++it)
    //{
    //    cSPIN spin0=(*it)[0];
    //    vec coeffs = zeeman(spin0, magB);
    //    _coeff_list.push_back(coeffs);
    //}
    vector< vector<cSPIN> > sag;
    sag = domain.getSpinAggregate();
    vector< vector<cSPIN> >::iterator it;
    for(it=sag.begin(); it!=sag.end(); ++it)
    {
        cSPIN spin0=(*it)[0];
        vec coeffs = zeeman(spin0, magB);
        _coeff_list.push_back(coeffs);
    }
}
ZeemanInteractionCoeff::~ZeemanInteractionCoeff()
{ //LOG(INFO) << "Default destructor: ZeemanInteractionCoeff.";
}
//}}}
//----------------------------------------------------------------------------//
//{{{ DipolarFieldInteractionCoeff
DipolarFieldInteractionCoeff::DipolarFieldInteractionCoeff(const cSpinInteractionDomain& domain, const cSPIN& center_spin, const PureState& state)
{ 
    _nCoeff = 6;

    vector< vector<cSPIN> > sag;
    sag = domain.getSpinAggregate();
    vector< vector<cSPIN> >::iterator it;
    for(it=sag.begin(); it!=sag.end(); ++it)
    {
        cSPIN spin0=(*it)[0];
        vec dip_field = dipole_field(spin0, center_spin, state.getVector() );
        vec coeffs; coeffs << dip_field[0] << dip_field[1] << dip_field[2] << 0.0 << 0.0 << 0.0;
        _coeff_list.push_back(coeffs);
    }
}
DipolarFieldInteractionCoeff::DipolarFieldInteractionCoeff(const cSpinInteractionDomain& domain, const vector<cSPIN>& spin_list, const vector<PureState>& state_list)
{
    _nCoeff = 6;

    vector< vector<cSPIN> > sag;
    sag = domain.getSpinAggregate();
    vector< vector<cSPIN> >::iterator it;
    for(it=sag.begin(); it!=sag.end(); ++it)
    {
        cSPIN spin0=(*it)[0];

        vec dip_field = zeros<vec>(3);
        for(int i=0; i<spin_list.size(); ++i)
            dip_field += dipole_field(spin0, spin_list[i], state_list[i].getVector() );
        
        vec coeffs; coeffs << dip_field[0] << dip_field[1] << dip_field[2] << 0.0 << 0.0 << 0.0;
        _coeff_list.push_back(coeffs);
    }
}
DipolarFieldInteractionCoeff::DipolarFieldInteractionCoeff(const cSpinInteractionDomain& domain, const vector<cSPIN>& spin_list, const vector<PureState>& state_list, const vec& pre_factor_list)
{
    _nCoeff = 6;

    vector< vector<cSPIN> > sag;
    sag = domain.getSpinAggregate();
    vector< vector<cSPIN> >::iterator it;
    for(it=sag.begin(); it!=sag.end(); ++it)
    {
        cSPIN spin0=(*it)[0];

        vec dip_field = zeros<vec>(3);
        for(int i=0; i<spin_list.size(); ++i)
            dip_field += pre_factor_list[i] * dipole_field(spin0, spin_list[i], state_list[i].getVector() );
        
        vec coeffs; coeffs << dip_field[0] << dip_field[1] << dip_field[2] << 0.0 << 0.0 << 0.0;
        _coeff_list.push_back(coeffs);
    }
}

DipolarFieldInteractionCoeff::~DipolarFieldInteractionCoeff()
{ //LOG(INFO) << "Default destructor: DipolarFieldInteractionCoeff";
}
//}}}
//----------------------------------------------------------------------------//
//{{{ PolarizationCoeff
PolarizationCoeff::PolarizationCoeff(const cSpinInteractionDomain& domain, const vector<vec>& pol)
{
    _nCoeff = 6;
    //for(auto pol_i:pol)
    //{
    //    vec coeffs = {pol_i[0], pol_i[1], pol_i[2], 0.0, 0.0, 0.0};
    //    _coeff_list.push_back(coeffs);
    //}
    for(int ii=0; ii< pol.size(); ++ii)
    {
        vec pol_i = pol[ii];
        vec coeffs;
        coeffs<< pol_i[0] << pol_i[1] << pol_i[2] << 0.0 << 0.0 << 0.0;
        _coeff_list.push_back(coeffs);
    }
}
PolarizationCoeff::~PolarizationCoeff()
{ //LOG(INFO) << "Default destructor: PolarizationCoeff.";
}
//}}}
////////////////////////////////////////////////////////////////////////////////

