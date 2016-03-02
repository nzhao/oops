#include "include/app/app.h"
#include "include/math/MatExp.h"


int  main(int argc, char* argv[])
{
    cx_double II = cx_double(0.0, 1.0);

    cx_mat H; 
    H << 6.0        << 1.0 + 1.0*II << endr
      << 1.0-1.0*II << 3.0;

    cx_double prefactor = 3.0;

    cout << "input data: " << endl;
    cout << prefactor * H << endl;
    cx_mat res;
    
    MatExp expM(H, prefactor, MatExp::ArmadilloExpMat);
    expM.run();
    res = expM.getResultMatrix();
    cout << "arma output: " << endl;
    cout << res << endl;
    
    MatExp expM2(H, prefactor, MatExp::PadeApproximation);
    expM2.run();

    res = expM2.getResultMatrix();
    cout << "pade output: " << endl;
    cout << res << endl;
    
    return 0;
}
