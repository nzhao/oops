#include "include/spin/SpinState.h"

////////////////////////////////////////////////////////////////////////////////
//{{{ cSpinState
cSpinState::cSpinState()
{ //LOG(INFO) << "Default constructor: cSpinState";
}

void cSpinState::make()
{
    this->cSpinInteraction::make();

    _space = Hilbert;
    //DIM_LIST dim_list = _sum_kron_prod.getKronProdList()[0].getDimList();
    MULTIPLIER c = 1.0;
    //for(auto d:dim_list) c*=d;
    for(int i=0; i<_dim_list.size(); ++i) c*=_dim_list[i];
    KronProd identity = KronProd(_dim_list);
    INDICES emptyIdx; TERM emptyTERM;
    identity.fill(emptyIdx, 1.0/c, emptyTERM);

    _sum_kron_prod.append(identity);
}

cSpinState::~cSpinState()
{ //LOG(INFO) << "Default destructor: cSpinState";
}
//}}}
////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
//{{{ SpinPolarization 
SpinPolarization::SpinPolarization()
{ //LOG(INFO) << "Default constructor: SpinPolarization";
}

SpinPolarization::SpinPolarization(const vector<cSPIN>& spin_list, const vector<int>& index_list, const vector<vec>& pol_list)
{
    _space = Hilbert;
    _spin_list=spin_list;

    _domain=SingleSpin(spin_list, index_list);
    _form=SingleSpinInteractionForm(_domain);
    _coeff=PolarizationCoeff(_domain, pol_list);

    make();
}

SpinPolarization::SpinPolarization(const vector<cSPIN>& spin_list, const vec& pol)
{
    _space = Hilbert;
    _spin_list=spin_list;
    vector<vec> pol_list;
    for(int i=0; i<spin_list.size(); ++i)
        pol_list.push_back(pol);

    _domain=SingleSpin(spin_list);
    _form=SingleSpinInteractionForm(_domain);
    _coeff=PolarizationCoeff(_domain, pol_list);

    make();
}

SpinPolarization::~SpinPolarization()
{ //LOG(INFO) << "Default destructor: SpinPolarization";
}
//}}}
////////////////////////////////////////////////////////////////////////////////
