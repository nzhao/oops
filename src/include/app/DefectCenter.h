#ifndef DEFECTCENTER_H
#define DEFECTCENTER_H
#include "include/oops.h"

class DefectCenter
{
public:
    DefectCenter() {};
    ~DefectCenter() {};
protected:
private:
};

class NVCenter:public DefectCenter
{
public:
    enum NVNuclearSpin { N14, N15 };
    
    NVCenter() : _nuc(N14) { make_spin();};
    NVCenter(NVNuclearSpin nuc ) {_nuc = nuc; make_spin();};
    ~NVCenter() {};

    cSPIN get_espin() const {return _electron;}
    cSPIN get_nspin() const {return _nitrogen;}
    cx_mat get_espin_hamiltonian() const {return _espin_hamiltonian;}
    PureState get_electron_spin_eigen_state(int i) const {return PureState(_eigen_vec_list[i]); }
    vec       get_spin_vector(int i) const {return _electron.get_spin_vector(_eigen_vec_list[i]); }

    void  set_magB(vec magB) {_magB = magB;}
    void  set_eleE(vec eleE) {_eleE = eleE;}
    void  make_espin_hamiltonian();

protected:
private:
    NVNuclearSpin _nuc;
    cSPIN         _nitrogen;
    cSPIN         _electron;
    vec           _magB;
    vec           _eleE;
    cx_mat        _espin_hamiltonian;
    vector<cx_vec> _eigen_vec_list;

    void make_spin();
};

//Adapter
class NVCenterElectronSpin : public cSPIN
{
public:
    NVCenterElectronSpin(const NVCenter& nv){ _nv = nv;};
    PureState get_eigenState(int i) const {return _nv.get_electron_spin_eigen_state(i);}
protected:
private:
    NVCenter _nv;
};
#endif
