#include "include/app/cce.h"
#include "include/misc/xmlreader.h"
#include <cstdlib>

_INITIALIZE_EASYLOGGINGPP

char PROJECT_PATH[500];
cSPINDATA SPIN_DATABASE=cSPINDATA();

int  main(int argc, char* argv[])
{
    char *env_path = std::getenv("CCE_PROJ_PATH");
    if(env_path!=NULL)
        strcpy(PROJECT_PATH, env_path);
    else
        getcwd(PROJECT_PATH, sizeof(PROJECT_PATH));

    char log_path[500];
    strcpy(log_path, PROJECT_PATH);
    strcat(log_path, "/dat/log/log.conf"); 

    char cfg_path[500];
    strcpy(cfg_path, PROJECT_PATH);
    strcat(cfg_path, "/dat/config/config.xml");
    ConfigXML cfg(cfg_path);
    cfg.printParameters();

    cout << "cut_off = " << cfg.getDoubleParameter("cut_off") << endl;
    cout << "maxorder = " << cfg.getIntParameter("maxorder") << endl;
    cout << "method = " << cfg.getStringParameter("method") << endl;


    _START_EASYLOGGINGPP(argc, argv);
    easyloggingpp::Configurations confFromFile(log_path);
    easyloggingpp::Loggers::reconfigureAllLoggers(confFromFile);

    
    int worker_num(0), my_rank(0);
    int mpi_status = MPI_Init(&argc, &argv);
    assert (mpi_status == MPI_SUCCESS);

    MPI_Comm_size(MPI_COMM_WORLD, &worker_num);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    LOG(INFO) << "my_rank = " << my_rank << "  vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv Program begins vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv"; 

    EnsembleCCE sol(my_rank, worker_num);
    sol.run();

    LOG(INFO) << "my_rank = " << my_rank << "  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ Program ends ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"; 

    mpi_status = MPI_Finalize();
    assert (mpi_status == MPI_SUCCESS);
}
