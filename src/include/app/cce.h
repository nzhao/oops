#include "include/oops.h"

extern char PROJECT_PATH[];

class CCE
{
public:
    CCE(): _my_rank(0), _worker_num(1) {};
    CCE(int my_rank, int worker_num, const string& config_file);
    ~CCE() {};
	void run();

    cSPIN           getCenterSpin() const {return _center_spin;}
    cSpinCollection getSpinCollecion() const {return _bath_spins;}
    cSpinCluster    getSpinClusters() const {return _my_clusters;}
    vector<mat>     getResultMatrix() const {return _cce_evovle_result;}
protected:
    ConfigXML        _cfg;
    vec              _center_spin_coord;
    string           _center_spin_isotope;
    char             _bath_spin_filename[500];
    char             _result_filename[500];
    double           _t0;
    double           _t1;
    int              _nTime;
    int              _cut_off_dist;
    int              _max_order;
    vec              _magB;

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
    cSPIN            create_center_spin();
    cSpinCollection  create_bath_spins();
    void             create_spin_clusters();
    void             job_distribution();
    void             run_each_clusters();
    void             DataGathering(const mat& resMat, int cce_order, int clst_num);

    virtual vec      cluster_evolution(int cce_order, int index)=0;
    void             post_treatment();
    void             cce_coherence_reduction();
    void             compuate_final_coherence();
    void             export_mat_file();
};

class EnsembleCCE:public CCE
{
public:
    EnsembleCCE(){};
    EnsembleCCE(int my_rank, int worker_num, const string& config_file);
    ~EnsembleCCE(){};
protected:

private:
    void set_parameters();
    vec cluster_evolution(int cce_order, int index);
    Hamiltonian create_spin_hamiltonian(const cSPIN& espin, const int spin_state, const vector<cSPIN>& spin_list);
    Liouvillian create_spin_liouvillian(const Hamiltonian& hami0, const Hamiltonian hami1);
    DensityOperator create_spin_density_state(const vector<cSPIN>& spin_list);
};
