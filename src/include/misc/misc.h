#ifndef MISC_H
#define MISC_H
#include <armadillo>
#include "include/spin/Spin.h"

using namespace arma;

extern cx_double II;

double spin_distance(const cSPIN& spin1, const cSPIN& spin2);

vec dipole(const cSPIN& spin1, const cSPIN& spin2);

vec r_vect(const cSPIN& obj1, const cSPIN& obj2);

vec zeeman(const cSPIN&spin, const vec& magB);
#endif
