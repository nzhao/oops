#include "include/app/DefectCenter.h"

NVCenter::NVCenter(NVNuclearSpin nuc)
{
    vec coord; coord << 0.0 << 0.0 << 0.0;
    add_spin_member(nuc, coord);
}

NVCenter::NVCenter(NVNuclearSpin nuc, vec coord)
{
    add_spin_member(nuc, coord);
}

void NVCenter::add_spin_member(NVNuclearSpin nuc, vec coord)
{
    double DIAMOND_LATTICE_CONST = 3.57; double dia4 = DIAMOND_LATTICE_CONST/4.0;
    vec nitrogen_coord; nitrogen_coord << coord[0] << coord[1] << coord[2];
    vec electron_coord; electron_coord << coord[0]+dia4 << coord[1]+dia4 << coord[2]+dia4;
    switch(nuc)
    {
        case N14:
            _nuclear_spin = cSPIN(nitrogen_coord, "14N");
            break;
        case N15:
            _nuclear_spin = cSPIN(nitrogen_coord, "15N");
            break;
        default:
            assert(0);
    }
    _electron_spin = cSPIN(electron_coord, "NVe");
}
void NVCenter::make_espin_hamiltonian()
{
    cx_mat sx=_electron_spin.sx();
    cx_mat sy=_electron_spin.sy();
    cx_mat sz=_electron_spin.sz();
    double omegaQ = _electron_spin.get_omegaQ() * 2.0 * datum::pi * 1e6;
    double gamma = _electron_spin.get_gamma();

    cx_mat s111 = 1.0/sqrt(3.0) * (sx + sy + sz);
    _espin_hamiltonian = omegaQ * s111 * s111 + gamma * (_magB(0)*sx + _magB(1)*sy + _magB(2)*sz );

    eig_sym(_eigen_vals, _eigen_vectors, _espin_hamiltonian);
}
