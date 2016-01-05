#ifndef SPIN_H
#define SPIN_H

#include <vector>
#include <string>
#include <armadillo>
#include "include/spin/SpinData.h"

using namespace std;
using namespace arma;

extern cSPINDATA SPIN_DATABASE;

////////////////////////////////////////////////////////////////////////////////
//{{{ cSPIN
/// This class creates spins with given coordinate and spin name string (e.g., nulcear spin isotope "13C", "1H", etc.).
/// The name string are defined in class cSPINDATA, or the properties can be set manually.
/// To use the data stored in the cSPINDATA class, you must define create an object named SPIN_DATABASE.
class cSPIN
{
public:
    cSPIN(vec coord, string isotope_str);///< Create a spin with given coordinate and spin name.
    cSPIN();///< default constructor: create a spin at [0, 0, 0], with name "None".

    //@{ 
    vec get_coordinate() { return coordinate; };
    string get_isotope() { return isotope; };
    int get_multiplicity() { return multiplicity; };
    int get_dimension() {return multiplicity; };
    double get_gamma() { return gamma; };
    double get_omegaQ() { return omegaQ; };
    double get_eta() { return eta; };
    //@}

    //@{
    void set_coordinate(vec coord) { coordinate = coord; };
    void set_isotope(string iso_str) { isotope =  iso_str; };
    void set_multiplicity(int val) { multiplicity = val; };
    void set_gamma(double val) {gamma =  val;};
    void set_omegaQ(double val) {omegaQ =  val;};
    void set_eta(double val) {eta =  val;};
    //@}
    
    //@{
    double S() {return 0.5*(multiplicity-1);}
    double S2() {return S()*(S()+1.0);}
    cx_mat sx();
    cx_mat sy();
    cx_mat sz();
    //@}
private:
    vec coordinate;
    string isotope;
    int multiplicity;
    double gamma;
    double omegaQ;
    double eta;
};
//}}}
////////////////////////////////////////////////////////////////////////////////
#endif
