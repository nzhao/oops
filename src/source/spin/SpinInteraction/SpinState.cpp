#include "include/spin/SpinState.h"

////////////////////////////////////////////////////////////////////////////////
//{{{ cSpinState
cSpinState::cSpinState()
{ LOG(INFO) << "Default constructor: cSpinState";}

void cSpinState::make()
{
    this->cSpinInteraction::make();

    DIM_LIST dim_list = _sum_kron_prod.getKronProdList()[0].getDimList();
    MULTIPLIER c = 1.0;
    for(auto d:dim_list) c*=d;
    KronProd identity = KronProd(dim_list);
    identity.fill( INDICES {}, 1.0/c, TERM {});

    _sum_kron_prod.append(identity);
}

cSpinState::~cSpinState()
{ LOG(INFO) << "Default destructor: cSpinState";}
//}}}
////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
//{{{ SpinPolarization 
SpinPolarization::SpinPolarization()
{ LOG(INFO) << "Default constructor: SpinPolarization";}

SpinPolarization::SpinPolarization(const vector<cSPIN>& spin_list, const vector<int>& index_list, const vector<vec>& pol_list)
{
    int nspin;
    _spin_list=spin_list;

    _domain=SingleSpin(spin_list, index_list);
    _form=SingleSpinInteractionForm(_domain);
    _coeff=PolarizationCoeff(_domain, pol_list);

    make();
}

SpinPolarization::~SpinPolarization()
{ LOG(INFO) << "Default destructor: SpinPolarization";}
//}}}
////////////////////////////////////////////////////////////////////////////////
