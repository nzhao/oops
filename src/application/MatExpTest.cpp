#include "include/app/app.h"
#include "include/math/MatExp.h"
#include <complex>
#include "include/math/krylov_expv.h"

cx_mat MAT;
cx_vec VEC;
vec    TIME_LIST;
SumKronProd SKP;

void prepare_data(string filename);
void test_small_mat();
cx_mat test_large_mat();
cx_mat test_large_mat_sparse();
cx_mat test_very_large_mat_CPU();
cx_mat test_very_large_mat_GPU();

int  main(int argc, char* argv[])
{
    string filename = "./dat/input/RoyCoord.xyz8";
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

    cx_double II = cx_double(0.0, 1.0);
    MAT = II*hami.getMatrix(); 
    cout << "hamiltonian mat generated." <<endl;

    SKP = hami.getKronProdForm();

    int dim = MAT.n_cols;
    PureState psi(dim);
    psi.setComponent(0, 1.0);
    VEC = psi.getVector();
    cout << "vector generated." <<endl;

    TIME_LIST = linspace<vec>(0.01, 0.1, 10);
}/*}}}*/

void test_small_mat()
{/*{{{*/
    cout << TIME_LIST.size() << endl;
    for(int i=0; i<TIME_LIST.size(); ++i)
    {
        MatExp expM(MAT, TIME_LIST(i), MatExp::ArmadilloExpMat);
        expM.run();
        
        MatExp expM2(MAT, TIME_LIST(i), MatExp::PadeApproximation);
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
    cout << "begin LARGE DENSE MAT" << endl;

    int dim = MAT.n_cols;
    std::complex<double> * mat = MAT.memptr();
    std::complex<double> * vecC = VEC.memptr();

    // krylov_zgexpv and krylov_zcooexpv;
    int      tn = TIME_LIST.size();
    double * ta = TIME_LIST.memptr();
    std::complex<double> *w_seq = new std::complex<double> [tn * dim];
    
    int err = krylov_zgexpv(dim, (double _Complex *)mat, (double _Complex *)vecC, tn, &ta[0], (double _Complex *)w_seq);

    cx_mat res(&w_seq[0], dim, tn);
    return res;
    
    // for debug;
    //cout << "krylov err = " << err << endl;
    //cout << "krylov_zgexpv works, i, ta:" << endl;
    //for (int i = 0; i < 10; i++) {
      //cout << i;
      //for (int j = 0; j < tn; j++)
        //cout << ", " << w_seq[i + j * dim];
      //cout << endl;
    //}
}/*}}}*/

cx_mat test_large_mat_sparse()
{/*{{{*/
    cout << endl;
    cout << "begin LARGE SPARSE MAT" << endl;
    sp_cx_mat mat_sparse = sp_cx_mat(MAT);
    sp_cx_mat::const_iterator start = mat_sparse.begin();
    sp_cx_mat::const_iterator end   = mat_sparse.end();

    std::complex<double> * vecC = VEC.memptr();

    int dim = MAT.n_cols;
    int nz = distance(start , end);

    int * ia = new int [nz];
    int * ja = new int [nz];
    std::complex<double> * a = new std::complex<double> [nz];
    int i=0;
    for(sp_cx_mat::const_iterator it = start; it != end; ++it)
    {
        ia[i] = it.row()+1;
        ja[i] = it.col()+1;
        a[i] = (*it);
        i++;
    }
    cout << nz << " non-zero elements" << endl;

    // krylov_zgexpv and krylov_zcooexpv;
    int      tn = TIME_LIST.size();
    double * ta = TIME_LIST.memptr();
    std::complex<double> *w_seq = new std::complex<double> [tn * dim];

    
    int err = krylov_zcooexpv(dim, nz, ia, ja, (double _Complex *)a, (double _Complex *)vecC, tn, &ta[0], (double _Complex *)w_seq);

    cx_mat res(&w_seq[0], dim, tn);
    return res;
    
    // for debug;
    cout << "krylov err = " << err << endl;
    cout << "krylov_zgexpv works, i, ta:" << endl;
    for (int i = 0; i < 10; i++) {
      cout << i;
      for (int j = 0; j < tn; j++)
        cout << ", " << w_seq[i + j * dim];
      cout << endl;
    }
}/*}}}*/

cx_mat test_very_large_mat_CPU()
{/*{{{*/
    cout << endl;
    cout << "Begin VERY_LARGE_MAT on CPU " << endl;
    cx_double prefactor = 1.0;

    MatExpVector expM(SKP, prefactor, VEC, TIME_LIST);  
    expM.run();

    cx_mat res = expM.getResult();
    return res;

}/*}}}*/

cx_mat test_very_large_mat_GPU()
{/*{{{*/
    cout << endl;
    cout << "Begin VERY_LARGE_MAT on GPU " << endl;
    cx_double prefactor = 1.0;

    MatExpVector expM(SKP, prefactor, VEC, TIME_LIST);  
    expM.run_gpu();

    cx_mat res = expM.getResult();
    return res;
}/*}}}*/
