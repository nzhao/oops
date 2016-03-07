#include "include/app/cce.h"
#include "include/misc/xmlreader.h"
#include <cstdlib>
#include "include/misc/lattice.h"

_INITIALIZE_EASYLOGGINGPP

string PROJECT_PATH;
string LOG_PATH;
string INPUT_PATH;
string OUTPUT_PATH;
string CONFIG_PATH;
string DEBUG_PATH;

cSPINDATA SPIN_DATABASE=cSPINDATA();
ConfigXML set_parameters(const string& xml_file_name);

int  main(int argc, char* argv[])
{
    //cout << "begin" << endl;
    //vector<int> base;
    //base.push_back(4);
    //base.push_back(3);
    //base.push_back(2);
    //vector<int> res = base_transform( 20, base);
    //for(int i=0; i<res.size(); ++i)
        //cout << res[i] << endl;
    //cout << "back = " << base_number(res, base) << endl;
    //cout << "end" << endl;
    int dim = 2;
    vec base1, base2;
    base1 << 1.0 << 0.0 << 0.0; base2 << 0.0 << 1.0 << 0.0;
    vector<vec> bases; bases.push_back(base1); bases.push_back(base2);
    int atom_num = 3;
    vec coord1, coord2, coord3;
    coord1 << 0.0 << 0.0 << 0.0;
    coord2 << 0.3 << 0.6 << 0.0;
    coord3 << 0.6 << 0.3 << 0.0;
    cout << "abc" << endl;
    vector<vec> pos; pos.push_back(coord1); pos.push_back(coord2); pos.push_back(coord3);
    vector<double> latt_const;
    latt_const.push_back(1.0); latt_const.push_back(1.0);
    Lattice latt(dim, bases, latt_const, atom_num, pos);

    umat range; range << -3 << 3 << endr << -2 << 2;
    latt.setRange(range);
    for(int kk = 0; kk <72; ++kk)
    {
        vector<int> idx = latt.getIndex(kk);
        for(int i=0; i<idx.size(); ++i)
            cout << idx[i] << ",\t";
        cout << endl << latt.getCoordinate(kk) << endl;
    }
    return 0;



    ConfigXML cfg = set_parameters("SingleSampleCCE.xml");

    string log_file = LOG_PATH + cfg.getStringParameter("Data", "log_file");
    _START_EASYLOGGINGPP(argc, argv);
    easyloggingpp::Configurations confFromFile(log_file.c_str());
    easyloggingpp::Loggers::reconfigureAllLoggers(confFromFile);

    
    int worker_num(0), my_rank(0);
    int mpi_status = MPI_Init(&argc, &argv);
    assert (mpi_status == MPI_SUCCESS);

    MPI_Comm_size(MPI_COMM_WORLD, &worker_num);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    LOG(INFO) << "my_rank = " << my_rank << "  vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv Program begins vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv"; 

    // create defect center
    double x = cfg.getDoubleParameter("CenterSpin",  "coordinate_x");
    double y = cfg.getDoubleParameter("CenterSpin",  "coordinate_y");
    double z = cfg.getDoubleParameter("CenterSpin",  "coordinate_z");
    vec coord; coord << x << y << z;
    NVCenter nv(NVCenter::N14, coord);
    
    double magBx = cfg.getDoubleParameter("Condition",  "magnetic_fieldX");
    double magBy = cfg.getDoubleParameter("Condition",  "magnetic_fieldY");
    double magBz = cfg.getDoubleParameter("Condition",  "magnetic_fieldZ");
    nv.set_magB(magBx, magBy, magBz);
    nv.make_espin_hamiltonian();

    // CCE
    SingleSampleCCE sol(my_rank, worker_num, &nv, cfg);
    sol.run();

    LOG(INFO) << "my_rank = " << my_rank << "  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ Program ends ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"; 

    mpi_status = MPI_Finalize();
    assert (mpi_status == MPI_SUCCESS);
}

ConfigXML set_parameters(const string& xml_file_name)
{/*{{{*/
    char *env_path = std::getenv("CCE_PROJ_PATH");
    if(env_path!=NULL)
        PROJECT_PATH = env_path;
    else
    {
        char pwd[500];
        getcwd(pwd, sizeof(pwd));
        PROJECT_PATH = pwd;
    }

    LOG_PATH    = PROJECT_PATH + "/dat/log/";
    INPUT_PATH  = PROJECT_PATH + "/dat/input/";
    OUTPUT_PATH = PROJECT_PATH + "/dat/output/";
    CONFIG_PATH = PROJECT_PATH + "/dat/config/";
    DEBUG_PATH  = PROJECT_PATH = "/dat/debug/";

    ConfigXML cfg( CONFIG_PATH+xml_file_name );
    return cfg;
}/*}}}*/

