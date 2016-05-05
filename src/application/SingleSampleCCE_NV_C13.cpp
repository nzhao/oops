#include "include/app/app.h"
#include "include/app/single_sample_cce.h"

_INITIALIZE_EASYLOGGINGPP

ConfigXML set_parameters(const string& xml_file_name);

NVCenter create_defect_center(const ConfigXML& cfg);
cSpinSourceFromFile create_spin_source(const ConfigXML& cfg);
cDepthFirstPathTracing create_spin_cluster_algrithm(const ConfigXML& cfg, const cSpinCollection& bath_spins);

int  main(int argc, char* argv[])
{
    ConfigXML cfg = set_parameters("SingleSampleCCE_NV_C13.xml");

    string log_file = LOG_PATH + cfg.getStringParameter("Data", "log_file");
    _START_EASYLOGGINGPP(argc, argv);
    easyloggingpp::Configurations confFromFile(log_file.c_str());
    easyloggingpp::Loggers::reconfigureAllLoggers(confFromFile);

    
    int worker_num(0), my_rank(0);
    int mpi_status = MPI_Init(&argc, &argv);
    assert (mpi_status == MPI_SUCCESS);

    MPI_Comm_size(MPI_COMM_WORLD, &worker_num);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    LOG(INFO) << "my_rank = " << my_rank << "  vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv Program begins vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv"; 

    SingleSampleCCE sol(my_rank, worker_num, cfg);

    // Step 1: make a defect center
    NVCenter nv = create_defect_center(cfg);  
    sol.set_defect_center(&nv);

    // Step 2: make bath spins 
    cSpinSourceFromFile spin_file = create_spin_source(cfg);
    sol.set_bath_spin(&spin_file);
    
    // Step 3: make clusters
    cSpinCollection bath_spins = sol.getSpinCollecion();
    cDepthFirstPathTracing   dfpt = create_spin_cluster_algrithm(cfg, bath_spins);
    sol.set_bath_cluster(&dfpt);

    // Step 4: run_each_cluster 
    sol.prepare_bath_state();
    sol.run_each_clusters();

    // Step 5: post_treatment
    sol.post_treatment();


    LOG(INFO) << "my_rank = " << my_rank << "  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ Program ends ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"; 

    mpi_status = MPI_Finalize();
    assert (mpi_status == MPI_SUCCESS);
}


ConfigXML set_parameters(const string& xml_file_name)
{/*{{{*/
    char *env_path = std::getenv("CCE_PROJ_PATH");
    if(env_path!=NULL)
        PROJECT_PATH = env_path;
    else
    {
        char pwd[500];
        getcwd(pwd, sizeof(pwd));
        PROJECT_PATH = pwd;
    }

    LOG_PATH    = PROJECT_PATH + "/dat/log/";
    INPUT_PATH  = PROJECT_PATH + "/dat/input/";
    OUTPUT_PATH = PROJECT_PATH + "/dat/output/";
    CONFIG_PATH = PROJECT_PATH + "/dat/config/";
    DEBUG_PATH  = PROJECT_PATH + "/dat/debug/";

    ConfigXML cfg( CONFIG_PATH+xml_file_name );
    return cfg;
}/*}}}*/

NVCenter create_defect_center(const ConfigXML& cfg)
{/*{{{*/
    vec coord = cfg.getVectorParameter("CenterSpin", "coordinate");
    vec magB = cfg.getVectorParameter("Condition", "magnetic_field");

    NVCenter nv(NVCenter::N14, coord);
    nv.set_magB(magB);
    nv.make_espin_hamiltonian();

    return nv;
}/*}}}*/

cSpinSourceFromFile create_spin_source(const ConfigXML& cfg)
{/*{{{*/
    string input_filename  = INPUT_PATH + cfg.getStringParameter("Data", "input_file");
    cSpinSourceFromFile spin_file(input_filename);
    return spin_file;
}/*}}}*/

cDepthFirstPathTracing create_spin_cluster_algrithm(const ConfigXML& cfg, const cSpinCollection& bath_spins)
{/*{{{*/
    double cut_off_dist = cfg.getDoubleParameter("SpinBath",   "cut_off_dist");
    int    max_order    = cfg.getIntParameter   ("SpinBath",   "max_order");
    sp_mat c = bath_spins.getConnectionMatrix(cut_off_dist);
    cDepthFirstPathTracing dfpt(c, max_order);
    return dfpt;
}/*}}}*/


