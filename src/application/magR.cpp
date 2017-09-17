#include "include/app/app.h"
#include "include/math/MatExp.h"
#include "expv/include/expv.h"
#include <complex>
#include <math.h>
#include <mpi.h>
#include <omp.h>

_INITIALIZE_EASYLOGGINGPP

namespace po = boost::program_options;

int NSPIN;
vector<vec> COEFF;
cx_vec INI;
vector<uvec> OPLST;
cx_vec FIN;
cx_double PREFACTOR = cx_double(0.0,-1.0);
vec TIME_LIST;

po::variables_map ParseCommandLineOptions(int argc, char* argv[]);
void  prepare_coeff(const po::variables_map& para);
void prepare_ini();
vec liou_coeff_list(const int nspin, const vec coeff_list_l, const vec coeff_list_r);
cx_vec gene_mixvec(const int nspin, const vector< pair<int,cx_mat> > p_lst);
void save_matrix(const cx_mat result, const string filename);

int  main(int argc, char* argv[])
{
    po::variables_map para = ParseCommandLineOptions(argc, argv);

    ////////////////////////////////////////////////////////////////////////////////
    //{{{  MPI Preparation
    int worker_num(0), my_rank(0);
    int mpi_status = MPI_Init(&argc, &argv);
    assert (mpi_status == MPI_SUCCESS);

    MPI_Comm_size(MPI_COMM_WORLD, &worker_num);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    //}}}
    ////////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////////////
    //{{{ LOGGING
    string log_file = LOG_PATH +  para["logfile"].as<string>();
    _START_EASYLOGGINGPP(argc, argv);
    easyloggingpp::Configurations confFromFile(log_file.c_str());
    easyloggingpp::Loggers::reconfigureAllLoggers(confFromFile);
    //}}}
    ////////////////////////////////////////////////////////////////////////////////

    LOG(INFO) << "my_rank = " << my_rank << "  vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv Program begins vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv";

    // Step 1: generate coeffiien
    prepare_coeff(para);

    // Step 2: generate initial state
    prepare_ini();

    // Step 3: set observables
    uvec op1;
    op1 << 2 << 0 << 0 << 1 << 0;//sigma0_x sigma1_x
    OPLST.push_back(op1);
    uvec op2;
    op2 << 2 << 0 << 1 << 1 << 1;//sigma0_y sigma1_y
    OPLST.push_back(op2);
    uvec op3;
    op3 << 2 << 0 << 2 << 1 << 2;//sigma0_z sigma1_z
    OPLST.push_back(op3);

    // Step 4: generate final state
    vector< pair<int,cx_mat> > polarize_list_E;
    FIN = gene_mixvec(NSPIN,polarize_list_E);
    FIN *= pow(2,NSPIN);
    
    ///*{{{*/ Stop 5.0: check input
    //printf("[debug] NSPIN = %d\n",NSPIN);
    //for(int i = 0; i < FIN.size(); i++){
    //    if( abs(FIN[i]) > 0.0 ){
    //        printf("[debug] FIN(%d) = %e\n",i,abs(FIN(i)));}}
    //for(int i = 0; i < INI.size(); i++){
    //    if( abs(INI[i]) > 0.0 ){
    //        printf("[debug] INI(%d) = %e\n",i,abs(INI(i)));}}
    //for(int i = 0; i < OPLST.size(); i++){
    //    cout << OPLST[i] << endl;}
    //cout << TIME_LIST << endl;
    //cout << PREFACTOR << endl;
    ///*}}}*/
    COEFF[0] = zeros<vec>(3*2*NSPIN + 2*NSPIN*(2*NSPIN-1)*9/2);
    COEFF[1] = zeros<vec>(3*2*NSPIN + 2*NSPIN*(2*NSPIN-1)*9/2);
    COEFF[2] = zeros<vec>(1);
    // Step 5: Evolution
    MatExpVector expM(2*NSPIN,COEFF,INI,OPLST,FIN,PREFACTOR,TIME_LIST,MatExpVector::LargeVecExpv);
    expM.run();
    cx_mat result = expM.getResult();
    cx_mat time = conv_to<cx_mat>::from(TIME_LIST);
    result = join_horiz(time,result);

    // Step 6: Save data
    string filename = OUTPUT_PATH + para["output"].as<string>();
    //save_Matrix(result,filename);

    LOG(INFO) << "my_rank = " << my_rank << "  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ Program ends ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^";

    ////////////////////////////////////////////////////////////////////////////////
    //{{{ MPI Finializing
    mpi_status = MPI_Finalize();
    assert (mpi_status == MPI_SUCCESS);
    //}}}
    ////////////////////////////////////////////////////////////////////////////////
}

po::variables_map ParseCommandLineOptions(int argc, char* argv[])
{/*{{{*/
    ////////////////////////////////////////////////////////////////////////////////
    //{{{ record command line options
    string output_filename("magR");
    string command_opt("");
    for(int i=1; i<argc; ++i)
        command_opt += argv[i];
    cout << "Command line options recieved: " << command_opt << endl;
    output_filename += command_opt;
    //}}}
    ////////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////////////
    //{{{  Set path
    char *env_path = std::getenv("magR_PROJ_PATH");
    if(env_path!=NULL)
    {
        PROJECT_PATH = env_path;
        cout << "PROJECT_PATH got from environment varialbe: ";
    }
    else
    {
        char pwd[500];
        getcwd(pwd, sizeof(pwd));
        PROJECT_PATH = pwd;
        cout << "PROJECT_PATH set as present working directory (pwd): ";
    }
    cout << PROJECT_PATH << endl;

    LOG_PATH    = PROJECT_PATH + "/dat/log/";
    INPUT_PATH  = PROJECT_PATH + "/dat/input/";
    OUTPUT_PATH = PROJECT_PATH + "/dat/output/";
    CONFIG_PATH = PROJECT_PATH + "/dat/config/";
    DEBUG_PATH  = PROJECT_PATH + "/dat/debug/";
    //}}}
    ////////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////////////
     //{{{ Parse options
    po::variables_map para;
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help,h", "Print help message")
        ("initialize,I", "Initialize a dat folder")

        ("input,i",          po::value<string>()->default_value("magR16E.xyz"), "Input file name")
        ("output,o",         po::value<string>()->default_value(output_filename), "Output .mat file of results")
        ("logfile,l",        po::value<string>()->default_value("magR.conf"), "Config. file of logging")

        ("pair1,P1",         po::value<string>()->default_value("20.0 12.8  0.0"), "RadicalPair spin1 position(in unit Angstrom.)")
        ("isotope1,E1",       po::value<string>()->default_value("E"), "RadicalPair spin1 tye")
        ("pair2,P2",         po::value<string>()->default_value("20.0 -12.8 0.0"), "RadicalPair spin2 position(in unit Angstrom.)")
        ("isotope2,E2",       po::value<string>()->default_value("E"), "RadicalPair spin2 tye")
        ("Gaps0,Es0",        po::value<double>()->default_value(0.0e6), "Radical pair enregy gap between |s> and |0>(in unit Hz.)")
        ("Gapsp,Esp",        po::value<double>()->default_value(10.0e6), "Radical pair energy gap between |s> and |+>(in unit Hz.)")
        ("Gapsn,Esn",        po::value<double>()->default_value(10.0e6), "Radical pair energy gap between |s> and |->(in unit Hz.)")
        ("bath_numer,m",     po::value<int>()->default_value(4), "Bath spins number")

        ("B_stren,B", po::value<double>()->default_value(0.3e-4), "Magnetic field strength(in unit Tesla.)")
        ("B_theta,theta",    po::value<double>()->default_value(0), "Magnetic field direction(in unit degree.)")
        ("B_phi,phi",        po::value<double>()->default_value(0), "Magnetic field direction(in unit degree.)")
        ("dephasing_rate,r", po::value<double>()->default_value(0.0), "Dephasing rate of bath spins(in unit of Hz.)")
        ("dephasing_axis,x", po::value<string>()->default_value("1.0 1.0 1.0"), "Dephasing axis of bath spins")

        ("nTime,n",          po::value<int>()->default_value(101), "Number of time points")
        ("start,s",          po::value<double>()->default_value(0.0), "Start time (in unit of sec.)")
        ("finish,f",         po::value<double>()->default_value(1.0e-5), "Finish time (in unit of sec.)")
        ;

    po::store(parse_command_line(argc, argv, desc), para);
    ifstream ifs("config.cfg");
    if (ifs!=NULL)
        po::store(parse_config_file(ifs,desc),para);
    else
        cout << "No configure file found!" << endl;
    po::notify(para);

    PrintVariableMap(para);

    if (para.count("help")) {
        cout << desc;
        exit(0);
    }
    if (para.count("initialize")) {
        string path = getenv("OOPS_PATH");
        string cmd = "cp -r " + path + "/src/dat_example dat" ;
        cout << cmd << endl;
        system( cmd.c_str() );
        cout << "default dat folder initialized." << endl;
        exit(0);
    }
    //}}}
    ////////////////////////////////////////////////////////////////////////////////
    return para;
}/*}}}*/

void prepare_coeff(const po::variables_map& para)
{/*{{{*/
    string filename = INPUT_PATH + para["input"].as<string>();
    cSpinSourceFromFile spin_file(filename);
    cSpinCollection spins(&spin_file);
    spins.make();

    vec r1( para["pair1"].as<string>() );
    string E1 = para["isotope1"].as<string>();
    vec r2(para["pair2"].as<string>());
    string E2 = para["isotope2"].as<string>();
    cSPIN RadicalPair1(r1, E1);
    cSPIN RadicalPair2(r2, E2);
    
    vector<cSPIN> sl_bath = spins.getSpinList();
    vector<cSPIN> sl_tot;
    sl_tot.push_back( RadicalPair1 );
    sl_tot.push_back( RadicalPair2 );
    int nbath = para["bath_numer"].as<int>();
    for(int i = 0; i < nbath; i++){
        sl_tot.push_back( sl_bath[i] );}
    NSPIN = nbath + 2;

    double normB = para["B_stren"].as<double>();
    double thetaB = para["B_theta"].as<double>();
    thetaB = thetaB * M_PI / 180;
    double phiB = para["B_phi"].as<double>();
    phiB = phiB * M_PI / 180;
    vec magB;
    magB <<  sin(thetaB) * cos(phiB) << sin(thetaB) * sin(phiB) << cos(thetaB);
    magB = normB * magB;

    vec coeff;
    for(int i = 0; i < NSPIN; i++)
        for(int j = i + 1; j < NSPIN; j++)
            coeff = join_vert( coeff, dipole( sl_tot[i], sl_tot[j] ) );
    for(int i = 0; i < NSPIN; i++){
        vec temp = zeeman( sl_tot[i], magB );
        vec temp1;
        temp1 << temp(0) << temp(1) << temp(2);
        coeff = join_vert( coeff, temp1 );
    }
    double Es0 = para["Gaps0"].as<double>();
    double Esp = para["Gapsp"].as<double>();
    double Esn = para["Gapsn"].as<double>();
    double T1;
    double T2;
    //coeff[] = ;
    //COEFF = liou_coeff_list(NSPIN,coeff,coeff);
    
    for(int i = 0; i < coeff.size(); i++)
        printf("[debug] coeff[%d] = %e\n",i,coeff[i]);
    
    //add Lindblad
    //COEFF[] = ;
}/*}}}*/

void prepare_ini()
{/*{{{*/
  vector< pair<int,cx_mat> > polarize_list_E;
  INI = gene_mixvec(NSPIN,polarize_list_E);

  vector< pair<int,cx_mat> > polarize_list_X;
  cx_mat spin0;
  spin0 << 0.0 << 1.0 << endr
        << 1.0 << 0.0 << endr;
  spin0 = 0.5 * spin0;
  cx_mat spin1 = -spin0;
  polarize_list_X.push_back( make_pair(0,spin0) );
  polarize_list_X.push_back( make_pair(1,spin1) );
  INI += gene_mixvec(NSPIN,polarize_list_X);

  cx_double im_unit = cx_double(0.0,1.0);
  vector< pair<int,cx_mat> > polarize_list_Y;
  spin0 << 0.0     << -1.0 * im_unit << endr
        << im_unit << 0.0            << endr;
  spin0 = 0.5 * spin0;
  spin1 = -spin0;
  polarize_list_Y.push_back( make_pair(0,spin0) );
  polarize_list_Y.push_back( make_pair(1,spin1) );
  INI += gene_mixvec(NSPIN,polarize_list_Y);

  vector< pair<int,cx_mat> > polarize_list_Z;
  spin0 << 1.0  << 0.0  << endr
        << 0.0  << -1.0 << endr;
  spin0 = 0.5 * spin0;
  spin1 = -spin0;
  polarize_list_Z.push_back( make_pair(0,spin0) );
  polarize_list_Z.push_back( make_pair(1,spin1) );
  INI += gene_mixvec(NSPIN,polarize_list_Z);
}/*}}}*/

vec liou_coeff_list(const int nspin, const vec coeff_list_l, const vec coeff_list_r)
{/*{{{*/
    // -1.0 * transpose(H) [*] E;
    // transpose(a [*] b) = transpose(a) [*] transpose(b)
    int single_term = 3 * 2 * nspin;
    int pair_term = 9 * (2 * nspin) * (2 * nspin -1)/2;
    int nterm = single_term + pair_term;
    vec coeff_list = zeros<vec>(nterm,1);
    int idx = 0;
    int idx1 = 0;
    for (int i = 0; i < (2 * nspin); i++){
        for(int j = i+1; j< (2 * nspin) ; j++){
            for(int k = 0; k < 9; k++){
                //printf("i = %4d, j = %4d, k = %4d, idx = %4d, idx1 = %4d\n",i,j,k,idx,idx1);//for debug
                if( i < nspin && j < nspin){
                    if( k % 2 == 1)
                        coeff_list[idx1] = coeff_list_l[idx];
                    else
                        coeff_list[idx1] = -1.0 * coeff_list_l[idx];//transpose(sigma_y) = - sigma_y
                    idx++;}
                idx1++;}
    }
    }

    for ( int i = 0; i < single_term; i++){
        //printf("i = %4d, idx = %4d, idx1 = %4d\n",i,idx,idx1);//for debug
        if( i < 3 * nspin ){
            if( i % 3 == 1 )
                coeff_list[idx1] = coeff_list_l[idx];
            else
                coeff_list[idx1] = -1.0 * coeff_list_l[idx];//transpose(sigma_y) = - sigma_y;
            idx++;}
            idx1++;
    }

    //E [*] H
    idx = 0;
    idx1 = 0;
    for (int i = 0; i < (2 * nspin); i++){
        for(int j = i + 1; j< (2 * nspin) ; j++){
            for(int k = 0; k < 9; k++){
                //printf("i = %4d, j = %4d, k = %4d, idx = %4d, idx1 = %4d\n",i,j,k,idx,idx1);//for debug
                if( ( i > (nspin - 1) ) && ( j > (nspin - 1) ) ){
                        coeff_list[idx1] += coeff_list_r[idx];
                    idx++;}
                idx1++;}
    }
    }

    for ( int i = 0; i < single_term; i++){
        //printf("i = %4d, idx = %4d, idx1 = %4d\n",i,idx,idx1);//for debug
        if( i > (3 * nspin - 1) ){
            coeff_list[idx1] += coeff_list_r[idx];
            idx++;}
            idx1++;
    }
    return coeff_list;
}/*}}}*/

cx_vec gene_mixvec(const int nspin, const vector< pair<int,cx_mat> > p_lst)
{/*{{{*/
    //reverse spin sequence
    vector< pair<int,cx_mat> > polarize_list(p_lst.size());
    for(int i=0; i < polarize_list.size(); i++){
        polarize_list[i].first = nspin - 1 - p_lst[i].first;
        polarize_list[i].second = p_lst[i].second;
    }
	int rank,nprocess;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&nprocess);

    int rank_n = 0;
    int idx = nprocess;
    while( idx > 1){
        idx = idx / 2;
        rank_n++;}

    vector<cx_mat> rank_list(rank_n);
    for (int i = 0; i < rank_n; i++)
        rank_list[i] = 0.5 * eye<cx_mat>(2,2);
    for (int i = 0; i < polarize_list.size(); i++)
        if( polarize_list[i].first < rank_n )
            rank_list[ polarize_list[i].first ] = polarize_list[i].second;
    cx_mat rank_mat = eye<cx_mat>(1,1);
    for (int i = 0; i < rank_n; i++)
        rank_mat = kron(rank_mat, rank_list[i] );

    vector<cx_mat> state_list(nspin - rank_n );
    cx_mat thermal = 0.5 * eye<cx_mat>(2,2);
    for (int i = rank_n; i < nspin; i++)
        state_list[i - rank_n] = thermal;
    for (int i = 0; i < polarize_list.size(); i++)
        if( !(polarize_list[i].first < rank_n) )
            state_list[ polarize_list[i].first -rank_n ] = polarize_list[i].second;
    cx_mat state_mat = ones<cx_mat>(1,1) ;
    for (int i = rank_n; i < nspin; i++)
        state_mat = kron(state_mat, state_list[i - rank_n] );

    int n = nprocess;
    cx_mat rank_vec = reshape(rank_mat, rank_mat.n_rows * rank_mat.n_rows / n, n);
    cx_vec state_vec = reshape(kron(rank_vec.col(rank),state_mat),rank_mat.n_rows * state_mat.n_rows * state_mat.n_rows,1);

    return state_vec;
}/*}}}*/

void save_Matrix(const cx_mat result, const string name)
{/*{{{*/
    mat m_r = real(result);
    mat m_i = imag(result);
    int n_rows = result.n_rows;
    int n_cols = result.n_cols;

    mxArray *pArray = mxCreateDoubleMatrix(n_rows, n_cols, mxCOMPLEX);

    int dim = n_rows * n_cols;
    memcpy((void *)(mxGetPr(pArray)), (void *) m_r.memptr(), dim*sizeof(double));
    memcpy((void *)(mxGetPi(pArray)), (void *) m_i.memptr(), dim*sizeof(double));

    string dbg_filename = DEBUG_PATH + name + ".mat";
    MATFile *mFile = matOpen(dbg_filename.c_str(), "w");
    matPutVariableAsGlobal(mFile, name.c_str(), pArray);
    matClose(mFile);
    printf("[output] matrix export success!\n");

    mxDestroyArray(pArray);

}/*}}}*/
