#ifndef SPIN_H
#define SPIN_H

#include <vector>
#include <string>
#include <armadillo>
#include "include/spin/SpinData.h"

using namespace std;

extern cSPINDATA SPIN_DATABASE;

class cSPIN
{
public:
    cSPIN(arma::vec coord, string isotope_str);
    cSPIN();

    arma::vec get_coordinate() { return coordinate; };
    string get_isotope() { return isotope; };
    int get_multiplicity() { return multiplicity; };
    double get_gamma() { return gamma; };

    void set_coordinate(arma::vec coord) { coordinate = coord; };
    void set_isotope(string iso_str) { isotope =  iso_str; };
    
private:
    arma::vec coordinate;
    string isotope;
    int multiplicity;
    double gamma;
};

#endif
