#include <vector>
#include <string>
#include <armadillo>
#include "include/spin/Spin.h"

using namespace std;

cSPIN::cSPIN()
{
    coordinate = arma::vec {0, 0, 0};
    isotope = "None";
}

cSPIN::cSPIN(arma::vec coord, string isotope_str)
{
    coordinate = coord;
    isotope = isotope_str;
    multiplicity = SPIN_DATABASE.getData(isotope_str).multiplicity;
    gamma = SPIN_DATABASE.getData(isotope_str).gamma;
}
