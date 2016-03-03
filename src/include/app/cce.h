#include "include/oops.h"
#include "include/app/DefectCenter.h"

extern string INPUT_PATH;
extern string OUTPUT_PATH;

////////////////////////////////////////////////////////////////////////////////
//{{{  CCE
class CCE
{
public:
    CCE(): _my_rank(0), _worker_num(1) {};
    CCE(int my_rank, int worker_num, DefectCenter* defect, const ConfigXML& cfg);
    ~CCE() {};
	void run();

    cSPIN           getCenterSpin() const {return _center_spin;}
    cSpinCollection getSpinCollecion() const {return _bath_spins;}
    cSpinCluster    getSpinClusters() const {return _my_clusters;}
    vector<mat>     getResultMatrix() const {return _cce_evovle_result;}
protected:
    ConfigXML        _cfg;
    DefectCenter*    _defect_center;
    string           _center_spin_name;
    int              _state_idx0;
    int              _state_idx1;
    pair<PureState, 
         PureState>  _state_pair;
    string             _bath_spin_filename;
    string             _result_filename;
    double           _t0;
    double           _t1;
    int              _nTime;
    vec              _time_list;
    double           _cut_off_dist;
    int              _max_order;
    vec              _magB;
    string           _pulse_name;
    int              _pulse_num;

    int              _my_rank;
    int              _worker_num;

    cSPIN            _center_spin;
    cSpinCollection  _bath_spins;
    cSpinCluster     _my_clusters;

    cSpinCluster     _spin_clusters;
    vector<mat>      _cce_evovle_result;
    vector<mat>      _cce_evovle_result_tilder;
    mat              _final_result;
    mat              _final_result_each_order;

private:
    virtual void     set_parameters()=0;
    void             prepare_center_spin();
    void             create_bath_spins();
    virtual void     prepare_bath_state()=0;
    void             create_spin_clusters();
    void             job_distribution();
    void             run_each_clusters();
    void             DataGathering(mat& resMat, int cce_order, int clst_num);

    virtual vec      cluster_evolution(int cce_order, int index)=0;
    //virtual vec      calc_observables(QuantumEvolutionAlgorithm* ker)=0;
    void             post_treatment();
    void             cce_coherence_reduction();
    void             compuate_final_coherence();
    void             export_mat_file();
};
//}}}
////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
//{{{  EnsembleCCE
class EnsembleCCE:public CCE
{
public:
    EnsembleCCE(){};
    EnsembleCCE(int my_rank, int worker_num, DefectCenter* defect, const ConfigXML& cfg):CCE(my_rank, worker_num, defect, cfg) {};
    ~EnsembleCCE(){};
protected:

private:
    void set_parameters();
    void prepare_bath_state();
    vec cluster_evolution(int cce_order, int index);
    Hamiltonian create_spin_hamiltonian(const cSPIN& espin, const PureState& center_spin_stat, const vector<cSPIN>& spin_list);
    Liouvillian create_spin_liouvillian(const Hamiltonian& hami0, const Hamiltonian hami1);
    DensityOperator create_spin_density_state(const vector<cSPIN>& spin_list);

    vec _bath_polarization;
    vec calc_observables(QuantumEvolutionAlgorithm* ker);
};
//}}}
////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
//{{{  SingleSampleCCE
class SingleSampleCCE:public CCE
{
public:
    SingleSampleCCE() {};
    SingleSampleCCE(int my_rank, int worker_num, DefectCenter* defect, const ConfigXML& cfg):CCE(my_rank, worker_num, defect, cfg) {};
    ~SingleSampleCCE() {};
protected:
private:
    void set_parameters();
    void prepare_bath_state();
    vec cluster_evolution(int cce_order, int index);
    Hamiltonian create_spin_hamiltonian(const cSPIN& espin, const PureState& center_spin_stat, const vector<cSPIN>& spin_list, const cClusterIndex& clstIndex);
    Liouvillian create_spin_liouvillian(const Hamiltonian& hami0, const Hamiltonian hami1);
    PureState create_cluster_state(const cClusterIndex& clstIndex);

    void cache_dipole_field();
    vec calc_observables(QuantumEvolutionAlgorithm* ker1, QuantumEvolutionAlgorithm* ker2);

    int _bath_state_seed;
    vector< vector<vec> > _dipole_field_data;
    vector<PureState>     _bath_state_list;

};
//}}}
////////////////////////////////////////////////////////////////////////////////

