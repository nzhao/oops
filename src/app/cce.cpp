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



int  main(int argc, char* argv[])
{
    _START_EASYLOGGINGPP(argc, argv);
    easyloggingpp::Configurations confFromFile("../src/logs/log.conf");  // Load configuration from file
    easyloggingpp::Loggers::reconfigureAllLoggers(confFromFile); // Re-configures all the loggers to current configuration file
    LOG(INFO) << "###################################################";
    LOG(INFO) << "Program begins."; 
    
    vector<double> coordinate; coordinate.push_back(1.0);coordinate.push_back(2.0);coordinate.push_back(3.0);
    string isotope="13C";

    cSPIN s1=cSPIN(coordinate, isotope);

    cout << s1.get_coordinate()[1] << "\t" << s1.get_isotope() << endl;
    cout << s1.get_multiplicity() << "\t" << s1.get_gamma() << "\t"  << s1.get_omegaQ() << "\t" << s1.get_eta() << endl;

    cSpinSourceFromFile spin_file("../bin/RoyCoord.xyz");
    cSpinCollection sc(&spin_file);

    sc.make();

    vector<cSPIN> sl=sc.getSpinList();

    for ( int i=0; i<sl.size(); ++i)
        sl[i].get_coordinate().t().print();

    mat m=sc.getDistanceMatrix();// m.print("m=:");
    sp_mat c=sc.getConnectionMatrix(8.0);

    cDepthFirstPathTracing dfpt(c, 1);
    cSpinCluster cluster(&dfpt);

    cluster.make();
    cout << cluster << endl;

    SpinDipolarInteraction dip(sl);

    vec magB; magB << 3.0e-4 << 2.0e-4 << 1.0e-4;
    SpinZeemanInteraction zee(sl, magB);

    Hamiltonian hami(sl);
    hami.addInteraction(dip);
    hami.addInteraction(zee);
    hami.make();

    cx_mat h = hami.getMatrix();
    cout << h << endl;
    hami.saveMatrix();

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
}
