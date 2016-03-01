#include <assert.h>
#include "include/spin/SpinInteraction.h"
#include "include/kron/KronProd.h"

using namespace std;


////////////////////////////////////////////////////////////////////////////////
//{{{ cSpinInteraction
cSpinInteraction::cSpinInteraction()
{ //LOG(INFO) << "Defaut constructor: cSpinInteraction.";
}

cSpinInteraction::cSpinInteraction(const vector<cSPIN>& spin_list)
{ //LOG(INFO) << "Constructor: cSpinInteraction with spin_list.";
    _spin_list=spin_list;
}

cSpinInteraction::~cSpinInteraction()
{ //LOG(INFO) << "Default destructor: cSpinInteraction.";
}

void cSpinInteraction::make()
{
    assert(_domain.getLength() == _form.getLength());
    assert(_domain.getLength() == _coeff.getLength());
    assert(_form.get_nTerm() == _coeff.get_nCoeff() );

    if( !_spin_list.empty() )
        for(int i=0; i<_spin_list.size(); ++i)
        {
            cSPIN spin=_spin_list[i];
            _dim_list.push_back( spin.get_dimension() );
        }
    //if( !_spin_list.empty() )
    //    for(auto spin: _spin_list)
    //        _dim_list.push_back( spin.get_dimension() );

    INDEX_LIST  idxList=_domain.getIndexList();
    MAT_LIST    matList=_form.getMatList();
    COEFF_LIST coefList=_coeff.getCoeffList();

    vector<KronProd> kronProd_list;
    size_t domainSize = _domain.getLength();
    int nTerm = _form.get_nTerm();
    for(int i=0; i<domainSize; ++i)
        for(int j=0; j<nTerm; ++j)
            if(coefList[i][j] != 0.0)
            {
                KronProd kp=KronProd(_dim_list);
                kp.fill( idxList[i], coefList[i][j], matList[i][j] );
                kronProd_list.push_back( kp );
            }
    _sum_kron_prod=SumKronProd(kronProd_list);
}

ostream&  operator << (ostream& outs, cSpinInteraction& interaction)
{
    cout << interaction._sum_kron_prod << endl;
    return outs;
}
//}}}
////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
//{{{ SpinDipolarInteraction 
SpinDipolarInteraction::SpinDipolarInteraction()
{ //LOG(INFO) << "Default constructor: SpinDipolarInteraction.";
}

SpinDipolarInteraction::SpinDipolarInteraction(const vector<cSPIN>& spin_list)
{ //LOG(INFO) << "Constructor: SpinDipolarInteraction with spin_list";
    _spin_list=spin_list;

    _domain=SpinPair(spin_list);
    _form=TwoSpinInteractionForm(_domain);
    _coeff=DipolarInteractionCoeff(_domain);
    
    make();
}

SpinDipolarInteraction::~SpinDipolarInteraction()
{ //LOG(INFO) << "Default destructor: SpinDipolarInteraction.";
}
//}}}
////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
//{{{ SpinZeemanInteraction 
SpinZeemanInteraction::SpinZeemanInteraction()
{ //LOG(INFO) << "Default constructor: SpinZeemanInteraction.";
}

SpinZeemanInteraction::SpinZeemanInteraction(const vector<cSPIN>& spin_list, const vec& magB)
{ //LOG(INFO) << "Constructor: SpinZeemanInteraction with spin_list and magB";

    _spin_list=spin_list;

    _domain=SingleSpin(spin_list);
    _form=SingleSpinInteractionForm(_domain);
    _coeff=ZeemanInteractionCoeff(_domain, magB);
    
    make();
}

SpinZeemanInteraction::~SpinZeemanInteraction()
{ //LOG(INFO) << "Default destructor: SpinZeemanInteraction.";
}
//}}}
////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
//{{{ DipolarField 
DipolarField::DipolarField()
{ //LOG(INFO) << "Default constructor: DipolarField";
}

DipolarField::DipolarField(const vector<cSPIN>& spin_list, const cSPIN& center_spin, const PureState& state)
{ //LOG(INFO) << "Constructor: DipolarField with center spin and spin state";

    _spin_list=spin_list;

    _domain=SingleSpin(spin_list);
    _form=SingleSpinInteractionForm(_domain);
    _coeff=DipolarFieldInteractionCoeff(_domain, center_spin, state);

    make();
}
DipolarField::DipolarField(const vector<cSPIN>& spin_list, const vector<cSPIN>& source_list, const vector<PureState>& state_list)
{
    _spin_list=spin_list;

    _domain=SingleSpin(spin_list);
    _form=SingleSpinInteractionForm(_domain);
    _coeff=DipolarFieldInteractionCoeff(_domain, source_list, state_list);

    make();
}
DipolarField::DipolarField(const vector<cSPIN>& spin_list, const vector<cSPIN>& source_list, const vector<PureState>& state_list, const uvec& exclude_idx)
{
    _spin_list=spin_list;

    _domain=SingleSpin(spin_list);
    _form=SingleSpinInteractionForm(_domain);

    vec mask = ones<vec>(source_list.size() );
    for(int i=0; i< exclude_idx.n_elem; ++i)
        mask( exclude_idx[i]) = 0.0;
    _coeff=DipolarFieldInteractionCoeff(_domain, source_list, state_list, mask);

    make();
}

DipolarField::~DipolarField()
{ //LOG(INFO) << "Default destructor: DipolarField";
}
//}}}
////////////////////////////////////////////////////////////////////////////////
