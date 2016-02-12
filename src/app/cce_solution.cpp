#include "include/app/cce.h"

_INITIALIZE_EASYLOGGINGPP

cSPINDATA SPIN_DATABASE=cSPINDATA();

int  main(int argc, char* argv[])
{
    _START_EASYLOGGINGPP(argc, argv);
    easyloggingpp::Configurations confFromFile("../src/logs/log.conf");  // Load configuration from file
    easyloggingpp::Loggers::reconfigureAllLoggers(confFromFile); // Re-configures all the loggers to current configuration file

    ////////////////////////////////////////////////////////////////////////////////
    //{{{ MPI_Initialization
    int worker_num(0), my_rank(0);
    int mpi_status = MPI_Init(&argc, &argv);
    assert (mpi_status == MPI_SUCCESS);

    MPI_Comm_size(MPI_COMM_WORLD, &worker_num);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    LOG(INFO) << "my_rank = " << my_rank << "  VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV Program begins VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV"; 
    //}}}
    ////////////////////////////////////////////////////////////////////////////////

    EnsembleCCE sol(my_rank, worker_num);
    sol.run();

    if(my_rank == 0)
    {
        vector<mat> resMatList=sol.getResultMatrix();
        for(int i=0; i<resMatList.size(); ++i)
        {
            cout << "order = " << i << endl;
            cout << resMatList[i].col(1) << endl;
        }
    }
    //cSPIN espin = sol.getCenterSpin();
    //cSpinCollection sc = sol.getSpinCollecion();
    //cSpinCluster cluster =  sol.getSpinClusters();

    ////////////////////////////////////////////////////////////////////////////////
    //{{{ MPI Finalization
    LOG(INFO) << "my_rank = " << my_rank << "  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ Program ends ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"; 

    mpi_status = MPI_Finalize();
    assert (mpi_status == MPI_SUCCESS);
    //}}}
    ////////////////////////////////////////////////////////////////////////////////

}
