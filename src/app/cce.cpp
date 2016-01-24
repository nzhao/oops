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

_INITIALIZE_EASYLOGGINGPP

using namespace std;
using namespace arma;

cSPINDATA SPIN_DATABASE=cSPINDATA();

cSPIN            create_e_spin();
cSpinCollection  create_bath_spins_from_file();
cSpinCluster     create_spin_clusters(const cSpinCollection& sc);
Hamiltonian      create_spin_hamiltonian(const cSPIN& espin, const vector<cSPIN>& spin_list);

int  main(int argc, char* argv[])
{
    _START_EASYLOGGINGPP(argc, argv);
    easyloggingpp::Configurations confFromFile("../src/logs/log.conf");  // Load configuration from file
    easyloggingpp::Loggers::reconfigureAllLoggers(confFromFile); // Re-configures all the loggers to current configuration file
    LOG(INFO) << "################################################### Program begins ###################################################"; 

    cSPIN espin = create_e_spin();

    cSpinCollection spin_collection = create_bath_spins_from_file();

    cSpinCluster spin_cluster = create_spin_clusters(spin_collection);

    cClusterIndex clst = spin_cluster.getCluster(2, 4);

    vector<cSPIN> spin_list = spin_collection.getSpinList(clst);

    Hamiltonian hami = create_spin_hamiltonian(espin, spin_list);

/*
    //Liouvillian lv(hami);
    //lv.saveMatrix();

    vec pol; pol << 0 << 0 << 1;
    vector<int> idx; idx.push_back(1); vector<vec> v_pol; v_pol.push_back(pol);
    SpinPolarization p(sl, idx, v_pol);

    DensityOperator ds(sl);
    ds.addStateComponent(p);
    ds.make();

    cx_mat rho=ds.getMatrix();

    PureState psi(8);
    psi.setComponent(1 , 1);


    cx_double i = cx_double(0, 1);
    cx_mat expH=expmat(0.1*i * h); 
    cx_vec res = expH*psi.getVector();
    cout << expH << res  << endl;

    SimpleFullMatrixVectorEvolution kernel(hami, psi);
    kernel.setTimeSequence( linspace<vec>(0.0, 1.0, 101) );

    QuantumEvolution dynamics(&kernel);
    dynamics.run();

    cx_vec st; st << 1 << 0;
    cout << dipole_field(sl[0], sl[1], st); 
    cx_vec st1; st1 << 0 << 1;
    cout << dipole_field(sl[0], sl[1], st1); 
    */
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
    sp_mat c=sc.getConnectionMatrix(8.0);

    cDepthFirstPathTracing dfpt(c, 3);
    cSpinCluster cluster(&dfpt);

    cluster.make();
    return cluster;
}
//}}}
////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
//{{{ Create spin Hamiltonian for a given cluster
Hamiltonian create_spin_hamiltonian(const cSPIN& espin, const vector<cSPIN>& spin_list)
{
    SpinDipolarInteraction dip(spin_list);

    vec magB; 
    magB << 3.0e-4 << 2.0e-4 << 1.0e-4;
    SpinZeemanInteraction zee(spin_list, magB);

    PureState center_spin_state(espin); 
    center_spin_state.setComponent(0, 1.0);
    DipolarField hf_field(spin_list, espin, center_spin_state);

    Hamiltonian hami(spin_list);
    hami.addInteraction(dip);
    hami.addInteraction(zee);
    hami.addInteraction(hf_field);
    hami.make();

    cx_mat h = hami.getMatrix();
    cout << h << endl;
    hami.saveMatrix();
    return hami;
}
//}}}
////////////////////////////////////////////////////////////////////////////////
