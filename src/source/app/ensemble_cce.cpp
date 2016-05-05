#include "include/app/ensemble_cce.h"


////////////////////////////////////////////////////////////////////////////////
//{{{  EnsembleCCE
void EnsembleCCE::set_parameters()
{
    string input_filename  = _cfg.getStringParameter("Data",       "input_file");
    string output_filename = _cfg.getStringParameter("Data",       "output_file");

    _state_idx0            = _cfg.getIntParameter   ("CenterSpin", "state_index0");
    _state_idx1            = _cfg.getIntParameter   ("CenterSpin", "state_index1");
    _center_spin_name      = _cfg.getStringParameter("CenterSpin", "name");

    _cut_off_dist          = _cfg.getDoubleParameter("SpinBath",   "cut_off_dist");
    _max_order             = _cfg.getIntParameter   ("SpinBath",   "max_order");
    _bath_polarization     = _cfg.getVectorParameter("SpinBath",   "bath_polarization");
    _bath_dephasing_rate   = _cfg.getDoubleParameter("SpinBath",   "dephasing_rate");
    _bath_dephasing_axis   = _cfg.getVectorParameter("SpinBath",   "dephasing_axis");

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
