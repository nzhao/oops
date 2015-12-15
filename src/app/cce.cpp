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
    cout << s1.get_multiplicity() << "\t" << s1.get_gamma() << endl;

    cSpinSourceFromFile spin_file("../bin/RoyCoord.xyz");
    cSpinCollection sc(&spin_file);

    sc.make();

    vector<cSPIN> sl=sc.getSpinList();

    for ( cSPIN s:sl)
    {
        s.get_coordinate().t().print();
    }

    mat m=sc.getDistanceMatrix(); m.print("m=:");
    umat c=sc.getConnectionMatrix(10.0); c.print("c=");


    cClusterIndex clst_idx1(uvec {1, 5, 2});
    cClusterIndex clst_idx2(uvec {2, 5, 1});
    cClusterIndex clst_idx3(uvec {1, 7, 9});

    cout << clst_idx1 << endl;
    cout << clst_idx2 << endl;
    cout << clst_idx3 << endl;

    CLST_IDX_LIST clst_idx_lst;
    clst_idx_lst.insert(clst_idx1);
    clst_idx_lst.insert(clst_idx2);
    clst_idx_lst.insert(clst_idx3);
    cout << "sizeof ... " << clst_idx_lst.size() << endl;

    cDepthFirstPathTracing dfpt(c, 10);
//    dfpt.insert_index_list(clst_idx_lst);
//    sp_mat XX=dfpt.get_cluster_mat(0);
//    mat XXFull(XX); XXFull.print("XX");

}
