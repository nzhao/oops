#include "include/app/cce.h"


////////////////////////////////////////////////////////////////////////////////
//{{{  SingleSampleCCE
class SingleSampleCCE:public CCE
{
public:
    SingleSampleCCE(int my_rank, int worker_num, const ConfigXML& cfg):CCE(my_rank, worker_num, cfg) { set_parameters(); }

    void set_parameters();
    void prepare_bath_state();
    vec cluster_evolution(int cce_order, int index);
protected:
private:
    Hamiltonian create_spin_hamiltonian(const cSPIN& espin, const PureState& center_spin_stat, const vector<cSPIN>& spin_list, const cClusterIndex& clstIndex);
    PureState create_cluster_state(const cClusterIndex& clstIndex);
    void cache_dipole_field();
    vec calc_observables(QuantumEvolutionAlgorithm* ker1, QuantumEvolutionAlgorithm* ker2);

    int _bath_state_seed;
    vector< vector<vec> > _dipole_field_data;
    vector<PureState>     _bath_state_list;

};
//}}}
////////////////////////////////////////////////////////////////////////////////
