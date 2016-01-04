#ifndef SPIN_H
#define SPIN_H

#include <vector>
#include <string>
#include <armadillo>
#include "include/spin/SpinData.h"

using namespace std;
using namespace arma;

extern cSPINDATA SPIN_DATABASE;

class cSPIN
{
public:
    cSPIN(vec coord, string isotope_str);
    cSPIN();

    vec get_coordinate() { return coordinate; };
    string get_isotope() { return isotope; };
    int get_multiplicity() { return multiplicity; };
    double get_gamma() { return gamma; };
    double get_omegaQ() { return omegaQ; };
    double get_eta() { return eta; };

    void set_coordinate(vec coord) { coordinate = coord; };
    void set_isotope(string iso_str) { isotope =  iso_str; };
    
    cx_mat sx();
    cx_mat sy();
    cx_mat sz();
private:
    vec coordinate;
    string isotope;
    int multiplicity;
    double gamma;
    double omegaQ;
    double eta;
};

#endif
