#include "include/misc/misc.h"
cx_double II = cx_double( 0.0, 1.0);
double DISTANCE_EPSILON = 1e-10;

double spin_distance(const cSPIN& spin1, const cSPIN& spin2) {
    return  norm(spin1.get_coordinate() - spin2.get_coordinate() ); }

vec r_vect(const cSPIN& obj1, const cSPIN& obj2){
    return obj1.get_coordinate() - obj2.get_coordinate(); };

vec dipole(const cSPIN& spin1, const cSPIN& spin2)
{
    double d=spin_distance(spin1, spin2);
    if(d<=DISTANCE_EPSILON)
    {
        vec res = zeros<vec>(9);
        return  res;
    }
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

vec zeeman(const cSPIN&spin, const vec& magB)
{
    double bx=magB[0], by=magB[1], bz=magB[2];
    double g=spin.get_gamma();
    double q=spin.get_omegaQ();
    double e=spin.get_eta();

    vec res;
    res << -g*bx <<  -g*by <<  -g*bz <<   e/3.0 <<  -e/3.0 << q;
    return res;
};

vec dipole_field(const cSPIN& spin, const cSPIN& source_spin, const cx_vec& source_state_vect)
{
    vec dip=dipole(spin, source_spin);
    mat dip_m = reshape(dip, 3, 3);
    vec s_vec = source_spin.get_spin_vector(source_state_vect);
    vec res = dip_m * s_vec;
    return res;
};

vector<double> Pulse_Timing(string pulsename, int n)
{
    vector<double> res;
    res.push_back( 0.0 );
    if( !strcmp(pulsename.c_str(), "CPMG") )
    {
        for(int i=0; i<n; ++i)
            res.push_back( (2.0*i+1.0) / (2.0*n) ); 
    }
    else
    {
        cout << "Pulse type not supported." << endl;
        assert(0);
    }
    res.push_back( 1.0 );
    return res;
}

vector<double> Pulse_Interval(string pulsename, int n)
{
    vector<double> res;
    vector<double> timings= Pulse_Timing(pulsename, n);
    for(int i=0; i<=n; ++i)
        res.push_back( timings[i+1] - timings[i] );
    return res;
}

vector<int> base_transform(int num, const vector<int>& base)
{
    vector<int> res;
    int max = 1;
    for(int i=0; i<base.size(); ++i)
        max *= base[i];
    if(num>max-1)
    {
        cout << "Too large input number" << endl;
        assert(0);
    }
    int i = base.size()-1, q = num;
    for(int i = base.size()-1; i>=0; --i)
    {
        res.push_back( q % base[i] );
        q = q / base[i];
    }
    reverse(res.begin(),res.end()); 
    return res;
}

int base_number(const vector<int>& num_in_base, const vector<int> base)
{
    int res=0, acc_base = 1;
    for(int i=base.size()-1; i>=0; --i)
    {
        res += num_in_base[i]*acc_base;
        acc_base *= base[i];
    }
    return res;
}
