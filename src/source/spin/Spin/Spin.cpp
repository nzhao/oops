#include <vector>
#include <string>
#include <armadillo>
#include "include/spin/Spin.h"

using namespace std;
using namespace arma;

////////////////////////////////////////////////////////////////////////////////
//{{{ cSPIN
cSPIN::cSPIN()
{
    coordinate.zeros(3);
    isotope = "None";
}

cSPIN::cSPIN(arma::vec coord, string isotope_str)
{
    coordinate = coord;
    isotope = isotope_str;
    multiplicity = SPIN_DATABASE.getData(isotope_str).multiplicity;
    gamma = SPIN_DATABASE.getData(isotope_str).gamma;
    omegaQ = SPIN_DATABASE.getData(isotope_str).omegaQ;
    eta = SPIN_DATABASE.getData(isotope_str).eta;
}

cx_mat cSPIN::sx() const
{
/// Spin operator Sx.
/// So far, only spin-1/2 is implemented.
    cx_mat res(2, 2);
    res(0, 0) = 0.0; res(0, 1) = 1.0;
    res(1, 0) = 1.0; res(1, 1) = 0.0;
    return 0.5*res;
}

cx_mat cSPIN::sy() const
{
/// Spin operator Sy.
/// So far, only spin-1/2 is implemented.
    cx_mat res(2, 2);
    res(0, 0) = 0.0;                 res(0, 1) = cx_double(0.0, -1.0);
    res(1, 0) = cx_double(0.0, 1.0); res(1, 1) = 0.0;
    return 0.5*res;
}

cx_mat cSPIN::sz() const
{
/// Spin operator Sz.
/// So far, only spin-1/2 is implemented.
    cx_mat res(2, 2);
    res(0, 0) = 1.0; res(0, 1) = 0.0;
    res(1, 0) = 0.0; res(1, 1) = -1.0;
    return 0.5*res;
}

vec cSPIN::get_spin_vector(const cx_vec& state) const
{
    vec res; 
    cx_mat x = state.t() * sx() * state;
    cx_mat y = state.t() * sy() * state;
    cx_mat z = state.t() * sz() * state;
    res << real( x[0] ) << real( y[0] ) << real( z[0] );
    return res;
}
//}}}
////////////////////////////////////////////////////////////////////////////////
