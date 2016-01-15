#ifndef MISC_H
#define MISC_H
#include <armadillo>

using namespace arma;

extern cx_double II;

template<class T> 
double distance(T& obj1, T& obj2) {
    return norm(obj1.get_coordinate() - obj2.get_coordinate()); };

template<class T> 
vec r_vect(T& obj1, T& obj2) {
    return obj1.get_coordinate() - obj2.get_coordinate(); };

template<class T>
vec dipole(T& spin1, T& spin2)
{
    double d=distance(spin1, spin2);
    vec    r=r_vect(spin1, spin2);
    vec    n=r/d;

    double nx = n[0];
    double ny = n[1];
    double nz = n[2];

    double d0=d*1e-10;
    double g1=spin1.get_gamma();
    double g2=spin2.get_gamma();
    double prefactor = datum::h_bar * (datum::mu_0)/(4.0 * datum::pi) * (g1*g2)/(d0*d0*d0);

    vec res;
    res << 1.0-3.0*nx*nx <<     -3.0*nx*ny <<     -3.0*nx*nz
        <<    -3.0*ny*nx <<  1.0-3.0*ny*ny <<     -3.0*ny*nz
        <<    -3.0*nz*nx <<     -3.0*nz*ny <<  1.0-3.0*nz*nz;
    return prefactor*res;
};
template<class T>
vec zeeman(T&spin, const vec& magB)
{
    double bx=magB[0], by=magB[1], bz=magB[2];
    double g=spin.get_gamma();
    double q=spin.get_omegaQ();
    double e=spin.get_eta();

    vec res;
    res << -g*bx <<  -g*by <<  -g*bz <<   e/3.0 <<  -e/3.0 << q;
    return res;
};
#endif
