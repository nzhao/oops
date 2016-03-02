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
        MatExp expM(_matrix, -1.0*dt*II, MatExp::PadeApproximation); expM.run();
        state_next = expM.getResultMatrix() * state_now;
        
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

    vector<cx_mat> expm_list, expm_list1;
    for(int j=0; j<_op_list.size(); ++j)
    {
        MatExp expM(_op_list[j].getMatrix(), -1.0*II* _time_segment[j]*dt, MatExp::PadeApproximation); expM.run();
        expm_list.push_back( expM.getResultMatrix() );
    }

    expm_list1 = expm_list;
    for(int i=1; i<_time_list.size(); ++i)
    {
        cx_vec state_i = _vector_list[0];
        for(int j=0; j<_op_list.size(); ++j)
        {
            state_i = expm_list1[j]*state_i; 
            expm_list1[j] = expm_list[j]*expm_list1[j];
        }
        _vector_list.push_back( state_i );
    }
}
//}}}
////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
//{{{ PiecewiseFullMatrixMatrixEvolution
PiecewiseFullMatrixMatrixEvolution::PiecewiseFullMatrixMatrixEvolution( const vector<QuantumOperator>& left_op_list, const vector<QuantumOperator>& right_op_list, const vector<double>& time_segment, const DensityOperator& ds)
{
   _left_op_list = left_op_list; 
   _right_op_list = right_op_list; 
   _time_segment = time_segment;
   _density_matrix = ds;
   _init_state = ds;
   _state_dimension = ds.getDimension()*ds.getDimension();
}

void PiecewiseFullMatrixMatrixEvolution::perform()
{
    _state_mat_list.push_back(_density_matrix.getMatrix());

    double dt = _time_list[1] - _time_list[0];

    vector<cx_mat> left_expm_list, right_expm_list, expm_list1, expm_list2;
    if (_left_op_list.size()!=_right_op_list.size()) assert(0);
    
    int op_num=_left_op_list.size();
    for(int j=0; j<op_num; ++j)
    {
        MatExp expM_left(_left_op_list[j].getMatrix(), -1.0*II* _time_segment[j]*dt, MatExp::PadeApproximation); expM_left.run();
        MatExp expM_right(_right_op_list[j].getMatrix(), 1.0*II* _time_segment[j]*dt, MatExp::PadeApproximation); expM_right.run();
        left_expm_list.push_back( expM_left.getResultMatrix() );
        right_expm_list.push_back( expM_right.getResultMatrix() );
    }

    expm_list1 = left_expm_list; expm_list2 = right_expm_list;
    for(int i=1; i<_time_list.size(); ++i)
    {
        cx_mat state_i = _density_matrix.getMatrix();
        for(int j=0; j<op_num; ++j)
        {
            state_i = expm_list1[op_num-1-j]*state_i*expm_list2[j];
            expm_list1[op_num-1-j] = left_expm_list[op_num-1-j]*expm_list1[op_num-1-j];
            expm_list2[j] = right_expm_list[j]*expm_list2[j];
        }
        
        _state_mat_list.push_back(state_i);
    }
}
//}}}
////////////////////////////////////////////////////////////////////////////////
