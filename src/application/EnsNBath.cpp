#include "include/app/app.h"
#include "include/app/ensemble_cce.h"

_INITIALIZE_EASYLOGGINGPP

namespace po = boost::program_options;
po::variables_map OPTIONS;        

void ParseCommandLineOptions(int argc, char* argv[]);
ConfigXML set_parameters(const string& xml_file_name);
NVCenter create_defect_center(const ConfigXML& cfg);
cSpinSourceUniformRandom create_spin_source(const ConfigXML& cfg);
cDepthFirstPathTracing create_spin_cluster_algrithm(const ConfigXML& cfg, const cSpinCollection& bath_spins);

int  main(int argc, char* argv[])
{
    ConfigXML cfg = set_parameters("EnsembleCCE_NV_E.xml");

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

    EnsembleCCE sol(my_rank, worker_num, cfg);

    // Step 1: make a defect center
    NVCenter nv = create_defect_center(cfg);  
    sol.set_defect_center(&nv);

    // Step 2: make bath spins 
    cSpinSourceUniformRandom spinUR = create_spin_source(cfg);
    sol.set_bath_spin(&spinUR);
    
    // Step 3: make clusters
    cSpinCollection bath_spins = sol.getSpinCollecion();
    cDepthFirstPathTracing   dfpt = create_spin_cluster_algrithm(cfg, bath_spins);
    sol.set_bath_cluster(&dfpt);

    // Step 4: run_each_cluster 
    sol.run_each_clusters();

    // Step 5: post_treatment
    sol.post_treatment();

    LOG(INFO) << "my_rank = " << my_rank << "  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ Program ends ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"; 

    mpi_status = MPI_Finalize();
    assert (mpi_status == MPI_SUCCESS);
}

void ParseCommandLineOptions(int argc, char* argv[])
{/*{{{*/
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help,h", "Print help message")
        ("input,i",          po::value<string>()->default_value("C13Bath/RoyCoord.xyz"), "Input .xyz file of bath spins")
        ("output,o",         po::value<string>()->default_value("Ensemble_NV_C13.mat"),  "Output .mat file of results")
        ("logfile,l",        po::value<string>()->default_value("EnsembleCCE.conf"),     "Config. file of logging")
        ("source,u",        po::value<string>()->default_value("file"),                  "Spin source")
        ("state0,a",         po::value<int>()->default_value(0),                         "Central spin state index - a")
        ("state1,b",         po::value<int>()->default_value(1),                         "Central spin state index - b")
        ("cce,c",            po::value<int>()->default_value(3),                         "CCE order")
        ("cutoff,d",         po::value<double>()->default_value(6.0),                    "Cut-off distance of bath spins")
        ("dephasing_rate,r", po::value<double>()->default_value(0.0),                    "dephasing rate of bath spins")
        ("dephasing_axis,x", po::value<string>()->default_value("1.0 1.0 1.0"),          "dephasing axis of bath spins")
        ("nTime,n",          po::value<int>()->default_value(101),                       "Number of time points")
        ("start,s",          po::value<double>()->default_value(0.0),                    "Start time (in unit of sec.)")
        ("finish,f",         po::value<double>()->default_value(0.002),                  "Finish time (in unit of sec.)")
        ("magnetic_field,B", po::value<string>()->default_value("0.1 0.1 0.1"),          "magnetic field vector in Tesla")
        ("pulse,p",          po::value<string>()->default_value("CPMG"),                 "Pulse name")
        ("pulse_num,m",      po::value<int>()->default_value(1),                         "Pulse number")
        ;

    po::store(parse_command_line(argc, argv, desc), OPTIONS);
    ifstream ifs("config.cfg");
    if (ifs!=NULL)
        po::store(parse_config_file(ifs,desc),OPTIONS);
    else
        cout << "No configure file found!" << endl;
    po::notify(OPTIONS);    

    if (OPTIONS.count("help")) {
        cout << desc;
        exit(0);
    }

    cout << OPTIONS["cce"].as<int>() << endl;
}/*}}}*/

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

cSpinSourceUniformRandom create_spin_source(const ConfigXML& cfg)
{/*{{{*/
    double   range = cfg.getDoubleParameter("SpinBath", "range");
    int        num = cfg.getIntParameter("SpinBath", "spin_number");
    int       seed = cfg.getIntParameter("SpinBath", "bath_seed");
    string isotope = cfg.getStringParameter("SpinBath", "isotope");

    cSpinSourceUniformRandom spinUR(range, num, isotope, seed);
    return spinUR;
}/*}}}*/

cDepthFirstPathTracing create_spin_cluster_algrithm(const ConfigXML& cfg, const cSpinCollection& bath_spins)
{/*{{{*/
    double cut_off_dist = cfg.getDoubleParameter("SpinBath",   "cut_off_dist");
    int    max_order    = cfg.getIntParameter   ("SpinBath",   "max_order");
    sp_mat c = bath_spins.getConnectionMatrix(cut_off_dist);
    cDepthFirstPathTracing dfpt(c, max_order);
    return dfpt;
}/*}}}*/
