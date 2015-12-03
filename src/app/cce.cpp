#include <vector>
#include <string>
#include <iostream>
#include "include/spin/Spin.h"
#include "include/spin/SpinData.h"
#include "include/spin/SpinCollection.h"
#include "include/spin/SpinSource.h"

using namespace std;

cSPINDATA SPIN_DATABASE=cSPINDATA();

int  main()
{

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
        cout << s.get_isotope() << "\t" 
            << "[" <<  s.get_coordinate()[0] 
            << ","  << s.get_coordinate()[1] 
            << ","  << s.get_coordinate()[2]  << "]" << endl;
    }
    return 0;
}
