#include "include/oops.h"

_INITIALIZE_EASYLOGGINGPP

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

    ////////////////////////////////////////////////////////////////////////////////
    //{{{ MPI_Initialization
    int worker_num(0), my_rank(0);
    int mpi_status = MPI_Init(&argc, &argv);
    assert (mpi_status == MPI_SUCCESS);

    MPI_Comm_size(MPI_COMM_WORLD, &worker_num);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    //}}}
    ////////////////////////////////////////////////////////////////////////////////

    size_t maxOrder = 3;
    int nTime = 101;

    uvec clstLength;     vector<umat> clstMat;
    double ** data = NULL;

    cSPIN espin = create_e_spin();
    cSpinCollection spin_collection = create_bath_spins_from_file();
    cSpinCluster spin_clusters;

    if(my_rank == 0)
    {/*{{{ cluster Generation */
        sp_mat c=spin_collection.getConnectionMatrix(6.0);
        cDepthFirstPathTracing dfpt(c, maxOrder);
        spin_clusters=cSpinCluster(spin_collection, &dfpt);
        spin_clusters.make();

        cout << spin_clusters << endl;

        data = new double * [maxOrder];
        for(int cce_order = 0; cce_order<maxOrder; ++cce_order)
            data[cce_order] = new double [nTime * spin_clusters.getClusterNum(cce_order)];
    }/*}}}*/

    ////////////////////////////////////////////////////////////////////////////////
    //{{{ Job distribution
    if(my_rank == 0)
    {
        spin_clusters.MPI_partition(worker_num);

        clstLength = spin_clusters.getMPI_ClusterLength(0);
        clstMat = spin_clusters.getMPI_Cluster(0);

        for(int i=1; i<worker_num; ++i)
        {
            uvec clstNum = spin_clusters.getMPI_ClusterLength(i);
            MPI_Send(clstNum.memptr(), maxOrder, MPI_UNSIGNED, i, 0, MPI_COMM_WORLD);

            vector<umat> clstMatList = spin_clusters.getMPI_Cluster(i);
            for(int j=0; j<maxOrder; ++j)
            {
                umat clstMat = clstMatList[j];
                MPI_Send(clstMat.memptr(), (j+1)*clstNum(j), MPI_UNSIGNED, i, j+1, MPI_COMM_WORLD);
            }
        }

    }
    else
    {
        unsigned int * clstLengthData = new unsigned int [maxOrder];
        MPI_Recv(clstLengthData, maxOrder, MPI_UNSIGNED, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        uvec tempV(clstLengthData, maxOrder);
        clstLength = tempV;
        delete [] clstLengthData;

        for(int j=0; j<maxOrder;++j)
        {
            unsigned int * clstMatData = new unsigned int [(j+1)*clstLength(j)];
            MPI_Recv(clstMatData, (j+1)*clstLength(j), MPI_UNSIGNED, 0, j+1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            umat tempM(clstMatData, clstLength(j), j+1);
            clstMat.push_back(tempM);
            delete [] clstMatData;
        }
    }
    //}}}
    ////////////////////////////////////////////////////////////////////////////////

    cSpinCluster my_clusters(spin_collection, clstLength, clstMat);

    for(int cce_order = 0; cce_order < maxOrder; ++cce_order)
    {
        cout << "my_rank = " << my_rank << ", " << "order  = " << cce_order << endl;
        size_t clst_num = my_clusters.getClusterNum(cce_order);

        mat resMat(nTime, clst_num, fill::ones);
        for(int i = 0; i < clst_num; ++i)
        {/*{{{ cluster evolution*/
            vector<cSPIN> spin_list = my_clusters.getCluster(cce_order, i);

            int spin_up = 0, spin_down = 1;
            Hamiltonian hami0 = create_spin_hamiltonian(espin, spin_up, spin_list);
            Hamiltonian hami1 = create_spin_hamiltonian(espin, spin_down, spin_list);

            Liouvillian lv = create_spin_liouvillian(hami0, hami1);

            DensityOperator ds = create_spin_density_state(spin_list);

            SimpleFullMatrixVectorEvolution kernel(lv, ds);
            kernel.setTimeSequence( linspace<vec>(0.0, 0.001, nTime) );

            ClusterCoherenceEvolution dynamics(&kernel);
            dynamics.run();

            resMat.col(i) = dynamics.calc_obs();
        }/*}}}*/

        ////////////////////////////////////////////////////////////////////////////////
        //{{{ Gathering Data
        if(my_rank != 0)
            MPI_Send(resMat.memptr(), nTime*clst_num, MPI_DOUBLE, 0, 100+my_rank, MPI_COMM_WORLD);
        else
        {
            memcpy(data[cce_order], resMat.memptr(), nTime*clst_num*sizeof(double));

            size_t prev_clst_num = clst_num;
            for(int source = 1; source < worker_num; ++source)
            {
                pair<size_t, size_t> pos = spin_clusters.getMPI_ClusterSize(cce_order, source);
                MPI_Recv(data[cce_order] + nTime*pos.first, nTime*(pos.second - pos.first), MPI_DOUBLE, source, 100+source, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }

        }
        //}}}
        ////////////////////////////////////////////////////////////////////////////////
    }

    if(my_rank == 0)
        post_treatment(data, spin_clusters, nTime);

    ////////////////////////////////////////////////////////////////////////////////
    //{{{ MPI Finalization
    mpi_status = MPI_Finalize();
    assert (mpi_status == MPI_SUCCESS);
    //}}}
    ////////////////////////////////////////////////////////////////////////////////
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
#ifdef HAS_MATLAB
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
#endif
}
//}}}
////////////////////////////////////////////////////////////////////////////////
