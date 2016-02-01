#include <stdlib.h>
#include <mat.h>
#include <vector>
#include <string>
#include <iostream>
#include "include/spin/Spin.h"
#include "include/spin/SpinData.h"
#include "include/spin/SpinCollection.h"
#include "include/spin/SpinSource.h"
#include "include/spin/SpinCluster.h"
#include "include/spin/SpinClusterAlgorithm.h"
#include "include/spin/SpinInteraction.h"
#include "include/kron/KronProd.h"

#include "include/easylogging++.h"
#include "include/misc/misc.h"
#include "include/quantum/HilbertSpaceOperator.h"
#include "include/quantum/LiouvilleSpaceOperator.h"
#include "include/quantum/MixedState.h"
#include "include/quantum/PureState.h"
#include "include/spin/SpinState.h"
#include "include/quantum/QuantumEvolutionAlgorithm.h"
#include "include/quantum/QuantumEvolution.h"

#include <mpi.h>

_INITIALIZE_EASYLOGGINGPP

using namespace std;
using namespace arma;

cSPINDATA SPIN_DATABASE=cSPINDATA();

cSPIN            create_e_spin();
cSpinCollection  create_bath_spins_from_file();
cSpinCluster     create_spin_clusters(const cSpinCollection& sc);
Hamiltonian      create_spin_hamiltonian(const cSPIN& espin, const int spin_state, const vector<cSPIN>& spin_list);
Liouvillian      create_spin_liouvillian(const Hamiltonian& hami0, const Hamiltonian hami1);
DensityOperator  create_spin_density_state(const vector<cSPIN>& spin_list);
void             post_treatment(double ** data, const cSpinCluster& spin_clusters, int nTime);


int  main(int argc, char* argv[])
{
    _START_EASYLOGGINGPP(argc, argv);
    easyloggingpp::Configurations confFromFile("../src/logs/log.conf");  // Load configuration from file
    easyloggingpp::Loggers::reconfigureAllLoggers(confFromFile); // Re-configures all the loggers to current configuration file
    LOG(INFO) << "################################################### Program begins ###################################################"; 

    // MPI head;
    int worker_num(0), my_rank(0);
    int mpi_status = MPI_Init(&argc, &argv);
    assert (mpi_status == MPI_SUCCESS);

    MPI_Comm_size(MPI_COMM_WORLD, &worker_num);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    cSpinCollection spin_collection = create_bath_spins_from_file();

    size_t maxOrder = 3;
    sp_mat c=spin_collection.getConnectionMatrix(6.0);
    cDepthFirstPathTracing dfpt(c, maxOrder);
    cSpinCluster spin_clusters(spin_collection, &dfpt);
    spin_clusters.make();

    cSPIN espin = create_e_spin();

    double ** data = NULL;
    if(my_rank == 0)
        data = new double * [maxOrder];

    int nTime = 101;
    for(int cce_order = 0; cce_order < maxOrder; ++cce_order)
    {
        size_t clst_num = spin_clusters.getClusterNum(cce_order);
        size_t job_num = clst_num % worker_num == 0 ? clst_num / worker_num : clst_num / worker_num + 1;

        mat resMat(nTime, job_num, fill::ones);
        for(int i = 0; i < job_num; ++i)
        {
            size_t index = my_rank*job_num + i;
            if( index < clst_num)
            {
                cout << "my_rank = " << my_rank << " cce_order =" << cce_order << ", index = "  << index << endl;

                vector<cSPIN> spin_list = spin_clusters.getCluster(cce_order, index);

                int spin_up = 0, spin_down = 1;
                Hamiltonian hami0 = create_spin_hamiltonian(espin, spin_up, spin_list);
                Hamiltonian hami1 = create_spin_hamiltonian(espin, spin_down, spin_list);

                Liouvillian lv = create_spin_liouvillian(hami0, hami1);

                DensityOperator ds = create_spin_density_state(spin_list);
>>>>>>> 1366a6907533d93b4162eaa075e06fbc451a1ffd

                SimpleFullMatrixVectorEvolution kernel(lv, ds);
                kernel.setTimeSequence( linspace<vec>(0.0, 0.001, nTime) );

                ClusterCoherenceEvolution dynamics(&kernel);
                dynamics.run();

                resMat.col(i) = dynamics.calc_obs();
            }
        }

        if(my_rank != 0)
            MPI_Send(resMat.memptr(), nTime*job_num, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        else
        {
            size_t blk_size = nTime*job_num;
            data[cce_order] = new double [blk_size * worker_num];
            
            memcpy(data[cce_order], resMat.memptr(), blk_size*sizeof(double));

            for(int source = 1; source < worker_num; ++source)
                MPI_Recv(data[cce_order] + source*blk_size, blk_size, MPI_DOUBLE, source, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        }

    }

    if(my_rank == 0)
        post_treatment(data, spin_clusters, nTime);

    // MPI initialization;
    mpi_status = MPI_Finalize();
    assert (mpi_status == MPI_SUCCESS);
}

////////////////////////////////////////////////////////////////////////////////
//{{{ Create an electrion spin
cSPIN create_e_spin()
{
    string isotope="E";

    double coord[] = {0.0, 0.0, 0.0};
    vector<double> coordinate(coord, coord+3);

    cSPIN espin=cSPIN(coordinate, isotope);
    return espin;
}
//}}}
////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
//{{{ Create bath spin list from xyz file
cSpinCollection create_bath_spins_from_file()
{
    cSpinSourceFromFile spin_file("../bin/RoyCoord.xyz");
    cSpinCollection sc(&spin_file);
    sc.make();
    return sc;
}
//}}}
////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
//{{{ Create spin clusters from a given spin list
cSpinCluster create_spin_clusters(const cSpinCollection& sc)
{
    sp_mat c=sc.getConnectionMatrix(6.0);

    cDepthFirstPathTracing dfpt(c, 3);
    cSpinCluster cluster(sc, &dfpt);

    cluster.make();

    return cluster;
}
//}}}
////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
//{{{ Create spin Hamiltonian for a given cluster
Hamiltonian create_spin_hamiltonian(const cSPIN& espin, const int spin_state, const vector<cSPIN>& spin_list)
{
    SpinDipolarInteraction dip(spin_list);

    vec magB; 
    magB << 0.0e-4 << 0.0e-4 << 1.0e-4;
    SpinZeemanInteraction zee(spin_list, magB);

    PureState center_spin_state(espin); 
    center_spin_state.setComponent(spin_state, 1.0);
    DipolarField hf_field(spin_list, espin, center_spin_state);

    Hamiltonian hami(spin_list);
    hami.addInteraction(dip);
    hami.addInteraction(zee);
    hami.addInteraction(hf_field);
    hami.make();
    return hami;
}
//}}}
////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
//{{{ Create Liouvillian operator from given Hamiltonians
Liouvillian create_spin_liouvillian(const Hamiltonian& hami0, const Hamiltonian hami1)
{
    Liouvillian lv0(hami0, SHARP);
    Liouvillian lv1(hami1, FLAT);

    Liouvillian lv = lv0 - lv1;
    //cx_mat lvMat = lv.getMatrix();
    return lv;
}
//}}}
////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
//{{{ Create spin density matrix
DensityOperator create_spin_density_state(const vector<cSPIN>& spin_list)
{
    vec pol = zeros<vec>(3);
    SpinPolarization p(spin_list, pol);

    DensityOperator ds(spin_list);
    ds.addStateComponent(p);
    ds.make();
    ds.makeVector();
    //cout << ds.getVector() << endl;
    //cout << ds.getMatrix() << endl;
    return ds;
}
//}}}
////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
//{{{ Post treatment
void post_treatment(double ** data, const cSpinCluster& spin_clusters, int nTime)
{
    cout << "begin post_treatement ... storing cce_data to file" << endl;

    size_t maxOrder = spin_clusters.getMaxOrder();
    for(int i=0; i<maxOrder; ++i)
    {
        char i_str [10];
        sprintf(i_str, "%d", i);
        string idx_str = i_str;
        string label = "cce_res_" + idx_str;
        string filename = label + ".mat";
        cout << "exporting " << filename << endl;

        size_t nClst = spin_clusters.getClusterNum(i);
        mxArray *pArray = mxCreateDoubleMatrix(nTime, nClst, mxREAL);

        size_t length= nTime * nClst;
        memcpy((void *)(mxGetPr(pArray)), (void *) data[i], length*sizeof(double));
    
        MATFile *mFile = matOpen(filename.c_str(), "w");
        matPutVariableAsGlobal(mFile, label.c_str(), pArray);
        matClose(mFile);

        mxDestroyArray(pArray);
    }
}
//}}}
////////////////////////////////////////////////////////////////////////////////
