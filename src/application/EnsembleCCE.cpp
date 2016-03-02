#include "include/app/app.h"
#include "include/app/cce.h"
#include "include/misc/xmlreader.h"
#include <cstdlib>

_INITIALIZE_EASYLOGGINGPP

ConfigXML set_parameters(const string& xml_file_name);

int  main(int argc, char* argv[])
{
    ConfigXML cfg = set_parameters("EnsembleCCE.xml");

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

    // create defect center
    double x = cfg.getDoubleParameter("CenterSpin",  "coordinate_x");
    double y = cfg.getDoubleParameter("CenterSpin",  "coordinate_y");
    double z = cfg.getDoubleParameter("CenterSpin",  "coordinate_z");
    vec coord; coord << x << y << z;
    NVCenter nv(NVCenter::N14, coord);
    
    double magBx = cfg.getDoubleParameter("Condition",  "magnetic_fieldX");
    double magBy = cfg.getDoubleParameter("Condition",  "magnetic_fieldY");
    double magBz = cfg.getDoubleParameter("Condition",  "magnetic_fieldZ");
    nv.set_magB(magBx, magBy, magBz);
    nv.make_espin_hamiltonian();

    // CCE
    EnsembleCCE sol(my_rank, worker_num, &nv, cfg);
    sol.run();

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
    DEBUG_PATH  = PROJECT_PATH = "/dat/debug/";

    ConfigXML cfg( CONFIG_PATH+xml_file_name );
    return cfg;
}/*}}}*/

