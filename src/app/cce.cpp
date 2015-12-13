#include <vector>
#include <string>
#include <iostream>
#include "include/spin/Spin.h"
#include "include/spin/SpinData.h"
#include "include/spin/SpinCollection.h"
#include "include/spin/SpinSource.h"
#include "include/spin/SpinCluster.h"
#include "include/spin/SpinGrouping.h"

#include "include/easylogging++.h"
#include "include/misc/misc.h"

INITIALIZE_EASYLOGGINGPP

using namespace std;

cSPINDATA SPIN_DATABASE=cSPINDATA();

int  main(int argc, char* argv[])
{
    START_EASYLOGGINGPP(argc, argv);
    LOG(INFO) << "My first info log using default logger";

    vector<double> coordinate {1.0, 2.0, 3.0};
    string isotope="13C";

    cSPIN s1=cSPIN(coordinate, isotope);

    cout << s1.get_coordinate()[1] << "\t" << s1.get_isotope() << endl;
    cout << s1.get_multiplicity() << "\t" << s1.get_gamma() << endl;

    cSpinSourceFromFile spin_file("RoyCoord.xyz");
    cSpinCollection sc(&spin_file);

    sc.make();

    vector<cSPIN> sl=sc.getSpinList();

    for ( cSPIN s:sl)
    {
        s.get_coordinate().t().print();
    }

    arma::mat m=sc.getDistanceMatrix();
    cout << m(0, 0) << "\t" << m(0, 1) <<"\t" << m(0, 2) << endl;
    cout << m(1, 0) << "\t" << m(1, 1) <<"\t" << m(1, 2) << endl;
    cout << m(2, 0) << "\t" << m(2, 1) <<"\t" << m(2, 2) << endl;

    arma::umat c=sc.getConnectionMatrix(10.0);
    cout << c(0, 0) << "\t" << c(0, 1) <<"\t" << c(0, 2) << endl;
    cout << c(1, 0) << "\t" << c(1, 1) <<"\t" << c(1, 2) << endl;
    cout << c(2, 0) << "\t" << c(2, 1) <<"\t" << c(2, 2) << endl;

    cDepthFirstPathTracing cg(c);
    cg.generate();

    vector<int> vIdx1 {1, 5, 2};
    vector<int> vIdx2 {1, 4, 3};

    cClusterIndex clst_idx1(vIdx1);
    cClusterIndex clst_idx2(vIdx2);

    cout << clst_idx1 << endl;
    cout << clst_idx2 << endl;

    cout << "compare: " <<  (clst_idx1 == clst_idx2) << endl;
    cout << "compare less: " <<  (clst_idx1 <  clst_idx2) << endl;

    cDepthFirstPathTracing dfpt(c);
    dfpt.generate();

    CLST_IDX_LIST clst_idx_lst;
    clst_idx_lst.insert(vIdx1);
    clst_idx_lst.insert(vIdx2);
    cout << "sizeof ... " << clst_idx_lst.size() << endl;

    for( auto idx:clst_idx_lst)
        cout << idx << "\t";
    cout << endl;
}
