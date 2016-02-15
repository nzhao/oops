#include "include/app/DefectCenter.h"

void NVCenter::make_spin()
{
    double DIAMOND_LATTICE_CONST = 3.57; double dia4 = DIAMOND_LATTICE_CONST/4.0;
    vec nitrogen_coord; nitrogen_coord << 0.0 << 0.0 << 0.0;
    vec electron_coord; electron_coord << dia4 << dia4 << dia4;
    switch(_nuc)
    {
        case N14:
            _nitrogen = cSPIN(nitrogen_coord, "14N");
            break;
        case N15:
            _nitrogen = cSPIN(nitrogen_coord, "15N");
            break;
        default:
            assert(0);
    }
    _electron = cSPIN(electron_coord, "NVe");
    _magB << 0.0 << 0.0 << 0.0;
    _eleE << 0.0 << 0.0 << 0.0;
}

void NVCenter::make_espin_hamiltonian()
{
    cx_mat sx=_electron.sx();
    cx_mat sy=_electron.sy();
    cx_mat sz=_electron.sz();
    double omegaQ = _electron.get_omegaQ() * 2.0 * datum::pi * 1e6;
    double gamma = _electron.get_gamma();

    cx_mat s111 = 1.0/sqrt(3.0) * (sx + sy + sz);
    _espin_hamiltonian = omegaQ * s111 * s111 + gamma * (_magB(0)*sx + _magB(1)*sy + _magB(2)*sz );

    vec eigval;
    cx_mat eigvec;
    eig_sym(eigval, eigvec, _espin_hamiltonian);
    for(int i=0; i<3; ++i)
        _eigen_vec_list.push_back( eigvec.col(i) );
}
