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

#include "include/easylogging++.h"
#include "include/misc/misc.h"

INITIALIZE_EASYLOGGINGPP

using namespace std;
using namespace arma;

cSPINDATA SPIN_DATABASE=cSPINDATA();

int  main(int argc, char* argv[])
{
    START_EASYLOGGINGPP(argc, argv);
    LOG(INFO) << "My first info log using default logger";

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
    dip.make();
    cout << dip << endl;

    vec magB={0.0, 0.0, 1.0};
    SpinZeemanInteraction zee(sl, magB);
    zee.make();
    cout << zee << endl;

    SumKronProd skp1=dip.getSumKronProd();
    SumKronProd skp2=zee.getSumKronProd();

    SumKronProd skp = skp1 + skp2;
//    cout << skp << endl;
}
