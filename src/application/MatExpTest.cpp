#include "include/app/app.h"
#include "include/math/MatExp.h"


int  main(int argc, char* argv[])
{
    cx_double II = cx_double(0.0, 1.0);

    cx_mat H; 
    H << 6.0        << 1.0 + 1.0*II << endr
      << 1.0-1.0*II << 3.0;

    cx_double prefactor = 3.0;

    MatExp expM(H, prefactor, MatExp::ArmadilloExpMat);
    expM.run();

    cx_mat res = expM.getResultMatrix();
    cout << res << endl;
    return 0;
}
