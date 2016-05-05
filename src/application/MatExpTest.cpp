#include "include/app/app.h"
#include "include/math/MatExp.h"
#include <complex>
#include "include/math/krylov_expv.h"

cx_mat MAT;
cx_vec VEC;
vec    TIME_LIST;
SumKronProd SKP;
cx_double PREFACTOR;

void prepare_data(string filename);
void test_small_mat();
cx_mat test_large_mat();
cx_mat test_large_mat_sparse();
cx_mat test_very_large_mat_CPU();
cx_mat test_very_large_mat_GPU();

int  main(int argc, char* argv[])
{
    string filename = "./dat/input/RoyCoord.xyz4";
    prepare_data(filename);

    test_small_mat();
    cx_mat res_large = test_large_mat();
    cx_mat res_large_sp = test_large_mat_sparse();
    cx_mat res_very_large_CPU = test_very_large_mat_CPU();
    cx_mat res_very_large_GPU = test_very_large_mat_GPU();

    cout << "diff 1 = " << norm(res_large_sp - res_large) << endl;
    cout << "diff 2 = " << norm(res_very_large_CPU - res_large) << endl;
    cout << "diff 3 = " << norm(res_very_large_GPU - res_large) << endl;
    cout << "diff 4 = " << norm(res_very_large_GPU - res_very_large_CPU) << endl;
    return 0;
}

void prepare_data(string filename)
{/*{{{*/
    cSpinSourceFromFile spin_file(filename);
    cSpinCollection spins(&spin_file);
    spins.make();

    vector<cSPIN> sl = spins.getSpinList();

    SpinDipolarInteraction dip(sl);
    Hamiltonian hami(sl);
    hami.addInteraction(dip);
    hami.make();

    Liouvillian lv0(hami);

    double rate = 1.0;
    vec axis; axis << 1.0 << 1.0 << 1.0;
    SpinDephasing dephasing(sl, rate, normalise(axis));
    LiouvilleSpaceOperator dephaseOperator(sl);
    dephaseOperator.addInteraction(dephasing);
    dephaseOperator.make();

    QuantumOperator lv = lv0 + dephaseOperator;

    PREFACTOR = cx_double(0.0, -1.0);
    MAT = lv.getMatrix(); 
    cout << "hamiltonian mat generated." <<endl;

    SKP =lv.getKronProdForm();

    int dim = MAT.n_cols;
    PureState psi(dim);
    psi.setComponent(0, 1.0);
    VEC = psi.getVector();
    cout << "vector generated." <<endl;

    TIME_LIST = linspace<vec>(0.01, 0.1, 10);
}/*}}}*/

void test_small_mat()
{/*{{{*/
    for(int i=0; i<TIME_LIST.size(); ++i)
    {
        MatExp expM(MAT, PREFACTOR*TIME_LIST(i), MatExp::ArmadilloExpMat);
        expM.run();
        
        MatExp expM2(MAT, PREFACTOR*TIME_LIST(i), MatExp::PadeApproximation);
        expM2.run();
        
        cx_mat resArma = expM.getResultMatrix();
        cx_mat resPade = expM2.getResultMatrix();
        cout << "t = " << TIME_LIST(i) << "; diff = " ;
        cout << norm(resArma-resPade) << endl; 
    }
}/*}}}*/

cx_mat test_large_mat()
{/*{{{*/
    cout << endl;
    cout << "begin LARGE DENSE MAT" <<  endl;

    MatExpVector expM(MAT, VEC, PREFACTOR, TIME_LIST);
    return expM.run();
}/*}}}*/

cx_mat test_large_mat_sparse()
{/*{{{*/
    cout << endl;
    cout << "begin LARGE SPARSE MAT" <<  endl;

    sp_cx_mat mat_sparse = sp_cx_mat(MAT);
    MatExpVector expM(mat_sparse, VEC, PREFACTOR, TIME_LIST);
    return expM.run();
}/*}}}*/

cx_mat test_very_large_mat_CPU()
{/*{{{*/
    cout << endl;
    cout << "Begin VERY_LARGE_MAT on CPU " <<  endl;

    MatExpVector expM(SKP, VEC, TIME_LIST, MatExpVector::Inexplicit);  
    return expM.run();
}/*}}}*/

cx_mat test_very_large_mat_GPU()
{/*{{{*/
    cout << endl;
    cout << "Begin VERY_LARGE_MAT on GPU " <<  endl;

    MatExpVector expM(SKP, VEC, TIME_LIST, MatExpVector::InexplicitGPU);  
    return expM.run();
}/*}}}*/
