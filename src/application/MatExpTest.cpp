#include "include/app/app.h"
#include "include/math/MatExp.h"
//#include <stdlib.h>
#include <complex>
#include "include/math/krylov_expv.h"

void test_small_mat();
void test_large_mat();
void test_large_mat_sparse();
void test_very_large_mat();

int  main(int argc, char* argv[])
{

    //test_small_mat();
    //test_large_mat();
    test_large_mat_sparse();
    //test_very_large_mat();
    return 0;
}

void test_small_mat()
{/*{{{*/
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
}/*}}}*/

void test_large_mat()
{/*{{{*/
    //please run this application in "oops/" direcotry 
    cSpinSourceFromFile spin_file("./dat/input/RoyCoord.xyz8");
    cSpinCollection spins(&spin_file);
    spins.make();

    vector<cSPIN> sl = spins.getSpinList();
    cout << "spin generated." << endl;

    SpinDipolarInteraction dip(sl);
    Hamiltonian hami(sl);
    hami.addInteraction(dip);
    hami.make();
    cout << "hamiltonian made." <<endl;

    cx_double II = cx_double(0.0, 1.0);
    cx_mat H = II*hami.getMatrix(); 
    int dim = H.n_cols;
    cout << "hamiltonian mat generated." <<endl;


    PureState psi(dim);
    psi.setComponent(0, 1.0);
    cx_vec v = psi.getVector();
    cout << "vector generated." <<endl;

    std::complex<double> * mat = H.memptr();
    std::complex<double> * vecC = v.memptr();

    // krylov_zgexpv and krylov_zcooexpv;
    int                   tn = 2;
    double                ta[2] = {0.5,1.0};
    std::complex<double>  *w_seq = NULL;
    int                   err;
    w_seq = new std::complex<double> [tn * dim];
    
    cout << "begin evolving" << endl;
    err = krylov_zgexpv(dim, (double _Complex *)mat, (double _Complex *)vecC, tn, &ta[0], (double _Complex *)w_seq);
    
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
void test_large_mat_sparse()
{
    //please run this application in "oops/" direcotry 
    cSpinSourceFromFile spin_file("./dat/input/RoyCoord.xyz8");
    cSpinCollection spins(&spin_file);
    spins.make();

    vector<cSPIN> sl = spins.getSpinList();
    cout << "spin generated." << endl;

    SpinDipolarInteraction dip(sl);
    Hamiltonian hami(sl);
    hami.addInteraction(dip);
    hami.make();
    cout << "hamiltonian made." <<endl;

    cx_double II = cx_double(0.0, 1.0);
    cx_mat H = II*hami.getMatrix(); 
    int dim = H.n_cols;
    cout << "hamiltonian mat generated." <<endl;

    cout << "Converting to Sparse matrix..." <<endl;
    sp_cx_mat H_sparse = sp_cx_mat(H);
    sp_cx_mat::const_iterator start = H_sparse.begin();
    sp_cx_mat::const_iterator end   = H_sparse.end();


    int nz = 0;
    for(sp_cx_mat::const_iterator it = start; it != end; ++it)
        nz ++;

    int * ia = new int [nz];
    int * ja = new int [nz];
    std::complex<double> * a = new std::complex<double> [nz];
    int i=0;
    for(sp_cx_mat::const_iterator it = start; it != end; ++it)
    {
        ia[i] = it.row();
        ja[i] = it.col();
        a[i] = (*it);
        i++;
    }
    for(i=0; i<nz; ++i)
        cout << ia[i] << "," << ja[i] << "\t" << a[i] << endl;
    cout << nz << " non-zero elements" << endl;


    PureState psi(dim);
    psi.setComponent(0, 1.0);
    cx_vec v = psi.getVector();
    std::complex<double> * vecC = v.memptr();
    cout << "vector generated." <<endl;

    // krylov_zgexpv and krylov_zcooexpv;
    int                   tn = 2;
    double                ta[2] = {0.5,1.0};
    std::complex<double>  *w_seq = NULL;
    int                   err;
    w_seq = new std::complex<double> [tn * dim];
    
    cout << "begin evolving" << endl;
    err = krylov_zcooexpv(dim, nz, ia, ja, (double _Complex *)a, (double _Complex *)vecC, tn, &ta[0], (double _Complex *)w_seq);
    
    // for debug;
    cout << "krylov err = " << err << endl;
    cout << "krylov_zgexpv works, i, ta:" << endl;
    for (int i = 0; i < 10; i++) {
      cout << i;
      for (int j = 0; j < tn; j++)
        cout << ", " << w_seq[i + j * dim];
      cout << endl;
    }
      
}

void test_very_large_mat()
{/*{{{*/
    //please run this application in "oops/" direcotry 
    cSpinSourceFromFile spin_file("./dat/input/RoyCoord.xyz4");
    cSpinCollection spins(&spin_file);
    spins.make();

    vector<cSPIN> sl = spins.getSpinList();

    SpinDipolarInteraction dip(sl);
    Hamiltonian hami(sl);
    hami.addInteraction(dip);
    hami.make();

    SumKronProd skp = hami.getKronProdForm();

    PureState psi( skp.getDim() );
    psi.setComponent(0, 1.0);
    cout << psi.getVector() << endl;

    vec time_list = linspace<vec>(0, 1, 11);
    MatExpVector expM(hami.getKronProdForm(), -1.0*II, psi.getVector(), time_list);  
    MatExpVector expMgpu(hami.getKronProdForm(), -1.0*II, psi.getVector(), time_list);  
    expM.run();
    expMgpu.run_gpu();

    cx_mat res = expM.getResult();
    cx_mat resGPU = expMgpu.getResult();
    cout << res << endl;
    cout << res - resGPU << endl;

    cx_mat H = hami.getMatrix(); 
    MatExp expM2(H, -1.0*II, MatExp::PadeApproximation);
    expM2.run();
    cx_mat resPade = expM2.getResultMatrix();

    cout << resPade*psi.getVector() << endl;
}/*}}}*/
