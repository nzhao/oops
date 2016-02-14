#include "include/quantum/QuantumEvolutionAlgorithm.h"


////////////////////////////////////////////////////////////////////////////////
//{{{ QuantumEvolutionAlgorithm
QuantumEvolutionAlgorithm::QuantumEvolutionAlgorithm(QuantumOperator& op, QuantumState& st)
{
    _operator = op;
    _init_state = st;
    _state_dimension = st.getDimension()*st.getDimension();
}
//}}}
////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
//{{{ SimpleFullMatrixVectorEvolution
void SimpleFullMatrixVectorEvolution::perform()
{
    _matrix = _operator.getMatrix();
    _vector_list.push_back(_init_state.getVector());
    
    ////////////////////////////////////////////////////////////////////////////////
    //begin evolution
    double dt;
    double t_now     = _time_list[0];
    cx_vec state_now = _vector_list[0];
    cx_vec state_next;

    for(int i=1; i<_time_list.size(); ++i)
    {
        double t_next = _time_list[i];
        dt = t_next - t_now;
        state_next = expmat( -1.0*dt*II*_matrix ) * state_now;
        
        _vector_list.push_back( state_next );
        
        t_now = t_next;
        state_now = state_next;
    }
}
//}}}
////////////////////////////////////////////////////////////////////////////////




////////////////////////////////////////////////////////////////////////////////
//{{{ PiecewiseFullMatrixVectorEvolution
PiecewiseFullMatrixVectorEvolution::PiecewiseFullMatrixVectorEvolution(const vector<QuantumOperator>& op_list, const vector<double>& time_segment, const QuantumState& st)
{
   _op_list = op_list; 
   _time_segment = time_segment;
   _init_state = st;
    _state_dimension = st.getDimension()*st.getDimension();
}

void PiecewiseFullMatrixVectorEvolution::perform()
{
    _vector_list.push_back(_init_state.getVector());

    double dt = _time_list[1] - _time_list[0];

    cx_mat expm1 = expmat(-1.0*II*dt*0.5*_op_list[0].getMatrix() );
    cx_mat expm2 = expmat(-1.0*II*dt*0.5*_op_list[1].getMatrix() );

    cx_mat mat1 = expm1, mat2=expm2;
    for(int i=1; i<_time_list.size(); ++i)
    {
        _vector_list.push_back( mat1*mat2*_vector_list[0] );
        mat1 = mat1*expm1;
        mat2 = mat2*expm2;
    }
    ////////////////////////////////////////////////////////////////////////////////
    //begin evolution
    //for(int i=1; i<_time_list.size(); ++i)
    //{
        //double t_i = _time_list[i];
        
        //cx_vec state = _vector_list[0];
        //for(int j=0; j<_op_list.size(); ++j)
        //{
            //double t_ij = t_i * _time_segment[j];
            //state = expmat( -1.0*t_ij*II*_op_list[j].getMatrix()) * state;
        //}
        //_vector_list.push_back( state );
    //}
}
//}}}
////////////////////////////////////////////////////////////////////////////////
