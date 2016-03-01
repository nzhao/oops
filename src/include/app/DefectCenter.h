#ifndef DEFECTCENTER_H
#define DEFECTCENTER_H
#include "include/oops.h"

class DefectCenter
{
public:
    DefectCenter() {};
    ~DefectCenter() {};

    cSPIN get_espin() const {return _electron_spin;}
    cSPIN get_nspin() const {return _nuclear_spin;}
    cx_vec get_eigen_state(int i) const {return _eigen_vectors.col(i);}

    virtual void make_espin_hamiltonian()=0;

    void  set_magB(double bx, double by, double bz) {_magB << bx << by << bz;}
    void  set_eleE(double ex, double ey, double ez) {_eleE << ex << ey << ez;}
protected:
    cSPIN _electron_spin;
    cSPIN _nuclear_spin;
    vec   _magB;
    vec   _eleE;
    cx_mat _espin_hamiltonian;
    cx_mat _eigen_vectors;
    vec _eigen_vals;
private:
    void  add_spin_member();
};

class NVCenter:public DefectCenter
{
public:
    enum NVNuclearSpin { N14, N15 };
    
    NVCenter(NVNuclearSpin nuc);
    NVCenter(NVNuclearSpin nuc, vec coord);
    ~NVCenter() {};
    
    void  make_espin_hamiltonian();

protected:
private:
        void  add_spin_member(NVNuclearSpin nuc, vec coord);
};
#endif
