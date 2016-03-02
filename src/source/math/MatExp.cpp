#include "include/math/MatExp.h" 

////////////////////////////////////////////////////////////////////////////////
//{{{  MatExp
MatExp::MatExp(const cx_mat& m, cx_double prefactor, MatExpMethod method)
{
    _matrix = m;
    _prefactor = prefactor;
    _method = method;
}

void MatExp::run()
{
    switch (_method) {
        case ArmadilloExpMat:
            _resMatrix = armadillo_exp_mat();
            break;
        case PadeApproximation:
            _resMatrix = pade_exp_mat();
            break;
        default:
            cout << "Exp method not sopport." << endl;
            assert(0);
    }
}

cx_mat MatExp::pade_exp_mat()
{
    cx_mat res;
    res = zeros<cx_mat>(10, 10);
    return res;
}
//}}}
////////////////////////////////////////////////////////////////////////////////




////////////////////////////////////////////////////////////////////////////////
//{{{  MatExpVector
MatExpVector::MatExpVector(const SumKronProd& skp, cx_double prefactor, const vec& time_list)
{
    _skp = skp;
    _prefactor = prefactor;
    _time_list = time_list;

    _nTime = time_list.n_elem;
    _dim = skp.getDim();
}

void MatExpVector::run()
{
    DIM_LIST   spin_dim = _skp.getDimList();  //spin_dim
    int      factor_num = _skp.getKronNum();  //nspin
    int           nTerm = _skp.getKronProdSize();  //nTerm

    vector<MULTIPLIER> coeff_list = _skp.getCoeffList();  //coeff_list
    vector<INDICES>  indices_list = _skp.getIndicesList(); 
    vector<TERM>        term_list = _skp.getTermList();
    vector<int>        nBody_list = _skp.getMatNumList(); //nBody_list

    INDICES full_indices = join_all(indices_list); //pos_list
    TERM       full_term = join_all(term_list);    //matC

    vector<int> dim_list; //dim_list
    for(int i=0; i<full_term.size(); ++i)
        dim_list.push_back( full_term[i].n_cols );

}
//}}}
////////////////////////////////////////////////////////////////////////////////
