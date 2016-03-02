#include "include/app/app.h"
#include "include/math/MatExp.h"

void test_large_mat();
void test_small_mat();

int  main(int argc, char* argv[])
{

    test_small_mat();
//    test_large_mat();
    return 0;
}

void test_small_mat()
{
    cx_double II = cx_double(0.0, 1.0);

    cx_mat H; 
    H << 6.0        << 1.0 + 1.0*II << endr
      << 1.0-1.0*II << 3.0;

    cx_double prefactor = 300.0;

    cout << "input data: " << endl;
    cout << II*prefactor * H << endl;
    cx_mat resArma, resPade;
    
    MatExp expM(H, II*prefactor, MatExp::ArmadilloExpMat);
    expM.run();
    resArma = expM.getResultMatrix();
    cout << "arma output: " << endl;
    cout << resArma << endl;
    
    MatExp expM2(H, II*prefactor, MatExp::PadeApproximation);
    expM2.run();

    resPade = expM2.getResultMatrix();
    cout << "pade output: " << endl;
    cout << resPade << endl;

    cout << "diff" << endl;
    cout << resArma-resPade << endl; 
}

void test_large_mat()
{
    //please run this application in "oops/" direcotry 
    cSpinSourceFromFile spin_file("./dat/input/magR16E.xyz");
    cSpinCollection spins(&spin_file);
    spins.make();

    vector<cSPIN> sl = spins.getSpinList();
    for(int i=0; i<sl.size(); ++i)
        cout << sl[i].get_coordinate() << endl;

    SpinDipolarInteraction dip(sl);
    Hamiltonian hami(sl);
    hami.addInteraction(dip);
    hami.make();

    SumKronProd skp = hami.getKronProdForm();

    vector<int> dim_list=skp.getDimList();
    for(int i=0; i<dim_list.size();++i)
        cout << dim_list[i] << endl;

    vector<KronProd> kp_list = skp.getKronProdList();
    for(int i=0; i<kp_list.size(); ++i)
        cout << kp_list[i] << endl;

    cout << kp_list.size() << endl;

}
