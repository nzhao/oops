#include "include/app/single_sample_cce.h"



////////////////////////////////////////////////////////////////////////////////
//{{{  SingleSampleCCE
void SingleSampleCCE::set_parameters()
{
    string input_filename  = _cfg.getStringParameter("Data",       "input_file");
    string output_filename = _cfg.getStringParameter("Data",       "output_file");

    _state_idx0            = _cfg.getIntParameter   ("CenterSpin", "state_index0");
    _state_idx1            = _cfg.getIntParameter   ("CenterSpin", "state_index1");
    _center_spin_name      = _cfg.getStringParameter("CenterSpin", "name");

    _cut_off_dist          = _cfg.getDoubleParameter("SpinBath",   "cut_off_dist");
    _bath_state_seed       = _cfg.getIntParameter   ("SpinBath",   "bath_state_seed");
    _max_order             = _cfg.getIntParameter   ("SpinBath",   "max_order");

    _nTime                 = _cfg.getIntParameter   ("Dynamics",   "nTime");
    _t0                    = _cfg.getDoubleParameter("Dynamics",   "t0"); 
    _t1                    = _cfg.getDoubleParameter("Dynamics",   "t1"); 

    _pulse_name            = _cfg.getStringParameter("Condition",  "pulse_name");
    _pulse_num             = _cfg.getIntParameter   ("Condition",  "pulse_number");

    _magB                  = _cfg.getVectorParameter("Condition",  "magnetic_field");

    _bath_spin_filename = INPUT_PATH + input_filename;
    _result_filename    = OUTPUT_PATH + output_filename;

    _time_list = linspace<vec>(_t0, _t1, _nTime);
}

void SingleSampleCCE::prepare_bath_state()
{
    vector<cSPIN> sl = _bath_spins.getSpinList();
    srand(_bath_state_seed);
    for(int i=0; i<sl.size(); ++i)
    {
        PureState psi_i(sl[i]);
        psi_i.setComponent( rand()%2, 1.0);
        _bath_state_list.push_back(psi_i);
    }
    //cache_dipole_field();
}

vec SingleSampleCCE::cluster_evolution(int cce_order, int index)
{
    vector<cSPIN> spin_list = _my_clusters.getCluster(cce_order, index);
    cClusterIndex clstIndex = _my_clusters.getClusterIndex(cce_order, index);

    Hamiltonian hami0 = create_spin_hamiltonian(_center_spin, _state_pair.first, spin_list, clstIndex);
    Hamiltonian hami1 = create_spin_hamiltonian(_center_spin, _state_pair.second, spin_list, clstIndex);

    vector<QuantumOperator> hm_list1 = riffle((QuantumOperator) hami0, (QuantumOperator) hami1, _pulse_num);
    vector<QuantumOperator> hm_list2 = riffle((QuantumOperator) hami1, (QuantumOperator) hami0, _pulse_num);
    vector<double> time_segment = Pulse_Interval(_pulse_name, _pulse_num);

    PureState psi = create_cluster_state(clstIndex);

    PiecewiseFullMatrixVectorEvolution kernel1(hm_list1, time_segment, psi);
    PiecewiseFullMatrixVectorEvolution kernel2(hm_list2, time_segment, psi);
    kernel1.setTimeSequence( _t0, _t1, _nTime);
    kernel2.setTimeSequence( _t0, _t1, _nTime);

    ClusterCoherenceEvolution dynamics1(&kernel1);
    ClusterCoherenceEvolution dynamics2(&kernel2);
    dynamics1.run();
    dynamics2.run();

    return calc_observables(&kernel1, &kernel2);
}

Hamiltonian SingleSampleCCE::create_spin_hamiltonian(const cSPIN& espin, const PureState& center_spin_state, const vector<cSPIN>& spin_list, const cClusterIndex& clstIndex )
{/*{{{*/
    SpinDipolarInteraction dip(spin_list);

    SpinZeemanInteraction zee(spin_list, _magB);

    DipolarField hf_field(spin_list, espin, center_spin_state);

    DipolarField bath_field(spin_list, _bath_spins.getSpinList(), _bath_state_list, clstIndex.getIndex() );

    Hamiltonian hami(spin_list);
    hami.addInteraction(dip);
    hami.addInteraction(zee);
    hami.addInteraction(hf_field);
    hami.addInteraction(bath_field);
    hami.make();
    return hami;
}/*}}}*/

PureState SingleSampleCCE::create_cluster_state(const cClusterIndex& clstIndex)
{/*{{{*/
    uvec idx = clstIndex.getIndex();
    cx_vec state_vec = _bath_state_list[ idx[0] ].getVector();
    for(int i=1; i<idx.n_elem; ++i)
        state_vec = kron( state_vec, _bath_state_list[ idx[i] ].getVector() );

    PureState res( state_vec );
    return res;
}/*}}}*/

void SingleSampleCCE::cache_dipole_field()
{/*{{{*/
    vector<cSPIN> sl = _bath_spins.getSpinList();
    for(int i=0; i<sl.size(); ++i)
    {
        vector<vec> dip_i;
        for(int j=0; j<sl.size(); ++j)
        {
            vec v_i = dipole_field(sl[i], sl[j], _bath_state_list[j].getVector() ); 
            dip_i.push_back( v_i );
        }
        _dipole_field_data.push_back(dip_i);
    }
}/*}}}*/

vec SingleSampleCCE::calc_observables(QuantumEvolutionAlgorithm* kernel1, QuantumEvolutionAlgorithm* kernel2)
{/*{{{*/
    vector<cx_vec>  state1 = kernel1->getResult();
    vector<cx_vec>  state2 = kernel2->getResult();
    cx_vec res = ones<cx_vec>(_nTime);
    for(int i=0; i<_nTime; ++i)
        res(i) =  cdot(state1[i], state2[i]);
    return real(res);
}/*}}}*/
//}}}
////////////////////////////////////////////////////////////////////////////////
