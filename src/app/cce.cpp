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

INITIALIZE_EASYLOGGINGPP

using namespace std;
using namespace arma;

cSPINDATA SPIN_DATABASE=cSPINDATA();



int  main(int argc, char* argv[])
{
    START_EASYLOGGINGPP(argc, argv);
    el::Configurations conf("/Users/nzhao/code/lib/active/oops/src/logs/log.conf");
    el::Loggers::reconfigureAllLoggers(conf);

    LOG(INFO) << endl;
    LOG(INFO) << "###################################################";
    LOG(INFO) << "Program begins."; 
    vector<double> coordinate {1.0, 2.0, 3.0};
    string isotope="13C";

    cSPIN s1=cSPIN(coordinate, isotope);

    cout << s1.get_coordinate()[1] << "\t" << s1.get_isotope() << endl;
    cout << s1.get_multiplicity() << "\t" << s1.get_gamma() << "\t"  << s1.get_omegaQ() << "\t" << s1.get_eta() << endl;

    cSpinSourceFromFile spin_file("../bin/RoyCoord.xyz");
    cSpinCollection sc(&spin_file);

    sc.make();

    vector<cSPIN> sl=sc.getSpinList();

    for ( cSPIN s:sl)
    {
        s.get_coordinate().t().print();
    }

    mat m=sc.getDistanceMatrix();// m.print("m=:");
    sp_mat c=sc.getConnectionMatrix(8.0);

    cDepthFirstPathTracing dfpt(c, 1);
    cSpinCluster cluster(&dfpt);

    cluster.make();
    cout << cluster << endl;

    SpinDipolarInteraction dip(sl);

    vec magB={0.0, 0.0, 1.0};
    SpinZeemanInteraction zee(sl, magB);

    Hamiltonian hami(sl);
    hami.addInteraction(dip);
    hami.addInteraction(zee);
    hami.make();

    Liouvillian lv(hami);
    lv.saveMatrix();

    vec pol = {0 , 0, 1};
    SpinPolarization p(sl, vector<int> {1}, vector<vec> {pol});

    DensityOperator ds(sl);
    ds.addStateComponent(p);
    ds.make();

    cx_mat rho=ds.getMatrix();

    PureState psi(16);
    psi.setComponent(1 , 1);

    cx_mat h = hami.getMatrix();
    cout << h << endl;

    cx_double i = cx_double(0, 1);
    cx_mat expH=expmat(0.1*i * h); 
    cout << expH << endl;
}
