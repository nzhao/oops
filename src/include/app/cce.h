#include "include/oops.h"
#include "include/app/DefectCenter.h"

extern string INPUT_PATH;
extern string OUTPUT_PATH;

namespace po = boost::program_options;

////////////////////////////////////////////////////////////////////////////////
//{{{  CCE
class CCE
{
public:
    CCE(int my_rank, int worker_num, const ConfigXML& cfg);
    CCE(int my_rank, int worker_num, const po::variables_map& para);

    void set_defect_center(DefectCenter* defect);
    void set_bath_spin(cSpinSource * source);
    void set_bath_cluster(cSpinGrouping * spin_grouping);
    void run_each_clusters();
    void post_treatment();

    cSPIN           getCenterSpin() const {return _center_spin;}
    cSpinCollection getSpinCollecion() const {return _bath_spins;}
    cSpinCluster    getSpinClusters() const {return _my_clusters;}
    vector<mat>     getResultMatrix() const {return _cce_evovle_result;}

protected:
    // I/O
    po::variables_map _para;
    ConfigXML        _cfg;
    string           _bath_spin_filename;
    string           _result_filename;
    
    // Defect Center
    DefectCenter*    _defect_center;
    string           _center_spin_name;
    cSPIN            _center_spin;
    int              _state_idx0;
    int              _state_idx1;
    pair<PureState, 
         PureState>  _state_pair;


    // Spin Bath
    double           _cut_off_dist;
    int              _max_order;
    cSpinSource*     _spin_source;
    cSpinGrouping*   _spin_grouping;
    cSpinCollection  _bath_spins;
    cSpinCluster     _spin_clusters;
    Lattice          _lattice;

    // Conditions
    vec              _magB;
    string           _pulse_name;
    int              _pulse_num;
    double           _t0;
    double           _t1;
    int              _nTime;
    vec              _time_list;

    // MPI
    int              _my_rank;
    int              _worker_num;
    cSpinCluster     _my_clusters;

    // Results
    vector<mat>      _cce_evovle_result;
    vector<mat>      _cce_evovle_result_tilder;
    mat              _final_result;
    mat              _final_result_each_order;

private:
    virtual void     set_parameters()=0;
    virtual vec      cluster_evolution(int cce_order, int index)=0;

    void             job_distribution();
    void             DataGathering(mat& resMat, int cce_order, int clst_num);
    void             cce_coherence_reduction();
    void             compuate_final_coherence();
    void             export_mat_file();
};
//}}}
////////////////////////////////////////////////////////////////////////////////


