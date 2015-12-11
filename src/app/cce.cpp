#include <vector>
#include <string>
#include <iostream>
#include "include/spin/Spin.h"
#include "include/spin/SpinData.h"
#include "include/spin/SpinCollection.h"
#include "include/spin/SpinSource.h"
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

    double x;
    x=distance( sl[1], sl[2] );
    cout << "dist = " << x  << endl;
    return 0;
}
