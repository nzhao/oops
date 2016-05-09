#include "include/app/ensemble_cce.h"


////////////////////////////////////////////////////////////////////////////////
//{{{  EnsembleCCE
void EnsembleCCE::set_parameters()
{
    string input_filename  = _para["input"].as<string>();
    string output_filename = _para["output"].as<string>();

    _state_idx0            = _para["state0"].as<int>();
    _state_idx1            = _para["state1"].as<int>();

    _max_order             = _para["cce"].as<int>();
    _cut_off_dist          = _para["cutoff"].as<double>();
    _bath_polarization     = vec( _para["polarization"].as<string>() );
    _bath_dephasing_rate   = _para["dephasing_rate"].as<double>();
    _bath_dephasing_axis   = vec( _para["dephasing_axis"].as<string>() );

    _nTime                 = _para["nTime"].as<int>();
    _t0                    = _para["start"].as<double>(); 
    _t1                    = _para["finish"].as<double>(); 

    _pulse_name            = _para["pulse"].as<string>();
    _pulse_num             = _para["pulse_num"].as<int>();
    _magB                  = vec( _para["magnetic_field"].as<string>() );
    
    _bath_spin_filename = INPUT_PATH + input_filename;
    _result_filename    = OUTPUT_PATH + output_filename + ".mat";
    _time_list = linspace<vec>(_t0, _t1, _nTime);
}

vec EnsembleCCE::cluster_evolution(int cce_order, int index)
{
    vector<cSPIN> spin_list = _my_clusters.getCluster(cce_order, index);
    
    Hamiltonian hami0 = create_spin_hamiltonian(_center_spin, _state_pair.first, spin_list);
    Hamiltonian hami1 = create_spin_hamiltonian(_center_spin, _state_pair.second, spin_list);
    LiouvilleSpaceOperator dephase = create_incoherent_operator(spin_list);

    QuantumOperator lvA = create_spin_liouvillian(hami0, hami1) + dephase;
    QuantumOperator lvB = create_spin_liouvillian(hami1, hami0) + dephase;

    vector<QuantumOperator> lv_list = riffle( lvA,  lvB, _pulse_num);
    DensityOperator ds = create_spin_density_state(spin_list);
    vector<double> time_segment = Pulse_Interval(_pulse_name, _pulse_num);

    PiecewiseFullMatrixVectorEvolution kernel(lv_list, time_segment, ds);
    kernel.setTimeSequence( _t0, _t1, _nTime);
    
    ClusterCoherenceEvolution dynamics(&kernel);
    dynamics.run();
    
    return calc_observables(&kernel);
}

Hamiltonian EnsembleCCE::create_spin_hamiltonian(const cSPIN& espin, const PureState& center_spin_state, const vector<cSPIN>& spin_list)
{/*{{{*/
    SpinDipolarInteraction dip(spin_list);

    SpinZeemanInteraction zee(spin_list, _magB);

    DipolarField hf_field(spin_list, espin, center_spin_state);

    Hamiltonian hami(spin_list);
    hami.addInteraction(dip);
    hami.addInteraction(zee);
    hami.addInteraction(hf_field);
    hami.make();
    return hami;
}/*}}}*/

LiouvilleSpaceOperator EnsembleCCE::create_incoherent_operator(const vector<cSPIN>& spin_list)
{/*{{{*/
    double rate = 2.0*datum::pi*_bath_dephasing_rate;
    vec axis = normalise(_bath_dephasing_axis);

    SpinDephasing dephasing(spin_list, rate, axis);
    LiouvilleSpaceOperator dephaseOperator(spin_list);
    dephaseOperator.addInteraction(dephasing);
    dephaseOperator.make();
    return dephaseOperator;

}/*}}}*/
Liouvillian EnsembleCCE::create_spin_liouvillian(const Hamiltonian& hami0, const Hamiltonian hami1)
{/*{{{*/
    Liouvillian lv0(hami0, SHARP);
    Liouvillian lv1(hami1, FLAT);
    Liouvillian lv = lv0 - lv1;
    return lv;
}/*}}}*/

DensityOperator EnsembleCCE::create_spin_density_state(const vector<cSPIN>& spin_list)
{/*{{{*/
    SpinPolarization p(spin_list, _bath_polarization);

    DensityOperator ds(spin_list);
    ds.addStateComponent(p);
    ds.make();
    ds.makeVector();
    return ds;
}/*}}}*/

vec EnsembleCCE::calc_observables(QuantumEvolutionAlgorithm* kernel)
{/*{{{*/
    int dim = kernel->getMatrixDim();
    cx_vec init_st = kernel->getInitalState();
    vector<cx_vec>  state = kernel->getResult();
    vec res = ones<vec>(_nTime);
    for(int i=0; i<_nTime; ++i)
        res(i) = real( cdot( dim*init_st, state[i]) );
    return res;
}/*}}}*/
//}}}
////////////////////////////////////////////////////////////////////////////////
