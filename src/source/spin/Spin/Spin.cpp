#include <vector>
#include <string>
#include <armadillo>
#include "include/spin/Spin.h"

using namespace std;
using namespace arma;

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

cx_mat cSPIN::sx()
{
    cx_mat res(2, 2);
    res(0, 0) = 0.0; res(0, 1) = 1.0;
    res(1, 0) = 1.0; res(1, 1) = 0.0;
    return res;
}

cx_mat cSPIN::sy()
{
    cx_mat res(2, 2);
    res(0, 0) = 0.0;                 res(0, 1) = cx_double(0.0, -1.0);
    res(1, 0) = cx_double(0.0, 1.0); res(1, 1) = 0.0;
    return res;
}

cx_mat cSPIN::sz()
{
    cx_mat res(2, 2);
    res(0, 0) = 1.0; res(0, 1) = 0.0;
    res(1, 0) = 0.0; res(1, 1) = -1.0;
    return res;
}
