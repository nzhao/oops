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
    //cSPINDATA SPIN_DATABASE=cSPINDATA();
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
/// So far, only spin-1/2 and spin-1 are implemented.
    cx_mat res;
    if(multiplicity == 2)
    {
        cx_mat res1(2, 2);
        res1<< 0.0 << 1.0 << endr
            << 1.0 << 0.0;
        res = 0.5*res1;
    }
    else if(multiplicity == 3)
    {
        cx_mat res1(3, 3);
        res1<< 0.0 << 1.0 << 0.0 << endr
            << 1.0 << 0.0 << 1.0 << endr
            << 0.0 << 1.0 << 0.0;
        res = res1/sqrt(2.0);
    }
    return res;
}

cx_mat cSPIN::sy() const
{
/// Spin operator Sy.
/// So far, only spin-1/2 and spin-1 are implemented.
    cx_mat res; cx_double II = cx_double(0.0, 1.0);
    if(multiplicity == 2)
    {
        cx_mat res1(2, 2);
        res1<< 0.0 << -II << endr
            << II  << 0.0;
        res = 0.5*res1;
    }
    else if(multiplicity == 3)
    {
        cx_mat res1(3, 3);
        res1<< 0.0 << -II << 0.0 << endr
            << II  << 0.0 << -II  << endr
            << 0.0 << II << 0.0;
        res = res1/sqrt(2.0);
    }
    return res;
}

cx_mat cSPIN::sz() const
{
/// Spin operator Sz.
/// So far, only spin-1/2 and spin-1 are implemented.
    cx_mat res; cx_double II = cx_double(0.0, 1.0);
    if(multiplicity == 2)
    {
        cx_mat res1(2, 2);
        res1<< 1.0 << 0.0 << endr
            << 0.0  << -1.0;
        res = 0.5*res1;
    }
    else if(multiplicity == 3)
    {
        cx_mat res1(3, 3);
        res1<< 1.0 << 0.0 << 0.0 << endr
            << 0.0  << 0.0 << 0.0  << endr
            << 0.0 << 0.0 << -1.0;
        res = res1;
    }
    return res;
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
