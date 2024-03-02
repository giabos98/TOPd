#pragma once

#include "CODE_HEADERS/codeHeader.h"

class PHYSICS
{
    public:

    //--- NS PROBLEM NAME
    std::string name;

    //--- DIMENSION ---
    int dim;

    //--- ABOUT TIME ---
    prec t_end;
    prec deltaT;
    prec deltaT_min;
    prec convergence_scale_factor;

    //--- CONSOLE LOG ---
    int completeLog;

    //--- PHYSICS PARAMETERS ---
    prec rho;
    prec mu;
    prec ni;

    //---  SOLVER FLAGS ---
    bool isStationary = true;
    bool flagBC; //0 penalty, 1 lifting 

    //--- V Struct -------
    int nNodes_v;
    int nElem_v;
    VECTOR h_v;
    MATRIX coord_v;
    MATRIX_INT elem_v;
    VECTOR_INT elem_geo_entities_ids_v; // note: the domains ids are scvaled with c++ notation (starting from 0)
    int max_geo_entity_id; 
    MATRIX Bloc_v; 
    MATRIX Cloc_v; 
    MATRIX Dloc_v;
    std::vector<MATRIX> Coef_v;
    VECTOR Volume_v;
    //--------------------

    //--- P struct -------
    int nElem;
    int nNodes;
    MATRIX coord;
    MATRIX_INT elem;
    MATRIX Bloc; MATRIX Cloc; MATRIX Dloc;
    VECTOR Volume;
    //Inlet Bound
    int n_inlet_bounds_elems;
    VECTOR_INT inlet_bounds;
    VECTOR_INT inlet_nodes_v;
    MATRIX_INT inlet_bounds_elems;
    VECTOR area_inlet_bounds_elems;
    //--------------------

    int nDof;

    // about boundaries
    int nBounds;
    std::string bound_info_file_path;
    std::string bound_nodes_v_file_path;

    // ABOUT SOLUTION TIMES
    VECTOR solution_times;
    VECTOR solution_deltaT;

    // PDE solution matrices for the current iteration
    MATRIX NS_solution;
    MATRIX ADJ_solution;

    VECTOR total_energy;
    MATRIX fluid_energy; //kinetic energy - pressure energy

    //------------
    prec f0Init = 1;
    prec func_normalization_factor = 1;
    prec gamma_change = 0;
    prec gamma_max = 0;
    int turn_on_gamma_acc = 0;
    prec V0;
    prec Vol;
    prec vol_fract;

    // Top. Opt. parameter filtering geometry

    PHYSICS(){};

    void print();

    void eval_solution_times();
};



//--------------------------------------------
// BASE CONSTRAINT CLASS
//--------------------------------------------
class CONSTRAINT
{
    public:

    //-----------
    // PARAMETERS
    //-----------
    int type;
    VECTOR parameters;

    // volume constraint
    prec Vr = -1.0;
    prec vol = -1.0;

    //edge size constraint
    int edge_id = -1;
    prec Sr = -1.0;
    prec side = -1.0;
    //---------------
    // END PARAMETERS
    //---------------

    //-------------
    // CONSTRUCTORS
    //-------------
    CONSTRAINT(){};

    CONSTRAINT(int constraint_type)
    {   
        int max_type = 0;
        if ((constraint_type < 0) || (constraint_type > max_type))
        {
            throw_line ("ERROR: invalid constraint type\n");
        }
        else
        {
            type = constraint_type;
        }
    }

    CONSTRAINT(int constraint_type, VECTOR constraint_parameters)
    {   
        int max_type = 0;
        if ((constraint_type < 0) || (constraint_type > max_type))
        {
            throw_line ("ERROR: invalid constraint type\n");
        }
        else
        {
            type = constraint_type;
        }
        parameters = constraint_parameters;
    }

    void initialize(int constraint_type);

    void initialize(int constraint_type, VECTOR constraint_parameters);

    void initialize_volume_constraint(int constr_type, prec volume_fraction);

    void initialize_edge_constraint(int constr_type, int edge, prec edge_fraction);
    //-----------------
    // END CONSTRUCTORS
    //-----------------
};
//--------------------------------------------
// END BASE CONSTRAINT CLASS
//--------------------------------------------

// //--------------------------------------------
// // VOLUME CONSTRAINT CLASS
// //--------------------------------------------
// class VOLUME_CONSTRAINT : public CONSTRAINT
// {
//     public:
//     //-------------
//     // CONSTRUCTORS
//     //-------------
//     VOLUME_CONSTRAINT(){};

//     VOLUME_CONSTRAINT(int constraint_type, VECTOR constraint_parameters) : CONSTRAINT(constraint_type)
//     {
//         // define maximum volume fraction
//         prec constraint_Vr = constraint_parameters[0];
//         if ((0.0 > constraint_Vr) || (constraint_Vr > 1.0))
//         {
//             std::cout << "\nVr: " << constraint_Vr << "\n";
//             throw_line ("ERROR: invalid maximum volume fraction\n");
//         }
//         else
//         {
//             Vr = constraint_Vr; // the Vr  must be inserted in the first index
//         }
//     }

//     void initialize(int constraint_type, VECTOR constraint_parameters);
//     //-----------------
//     // END CONSTRUCTORS
//     //-----------------
// };
// //--------------------------------------------
// // END VOLUME CONSTRAINT CLASS
// //--------------------------------------------

// //--------------------------------------------
// // EDGE SIZE CONSTRAINT CLASS
// //--------------------------------------------
// class EDGE_SIZE_CONSTRAINT : public CONSTRAINT
// {
//     public:
//     //-------------
//     // CONSTRUCTORS
//     //-------------
//     EDGE_SIZE_CONSTRAINT(){};

//     EDGE_SIZE_CONSTRAINT(int constraint_type, VECTOR constraint_parameters) : CONSTRAINT(constraint_type)
//     {
//         // define edge_id
//         int constraint_edge_id = int(constraint_parameters[0]);
//         edge_id = constraint_edge_id;

//         // define maximum edge size fraction
//         prec constraint_Sr = constraint_parameters[1];
//         if ((0.0 > constraint_Sr) || (constraint_Sr > 1.0))
//         {
//             throw_line ("ERROR: invalid maximum volume fraction\n");
//         }
//         else
//         {
//             Sr = constraint_Sr; // the Vr  must be inserted in the first index
//         }
//     }

//     void initialize(int constraint_type, VECTOR constraint_parameters);
//     //-----------------
//     // END CONSTRUCTORS
//     //-----------------
// };
// //--------------------------------------------
// // END EDGE SIZE CONSTRAINT CLASS
// //--------------------------------------------

//--------------------------------------------
// CONSTRAINTS HANDLER CLASS
//--------------------------------------------
class CONSTRAINTS
{
    public:

    //-------------
    // PARAMETERS
    //-------------
    int n_constr;
    int max_constraint_type = 0;
    int max_parameters_length = 2;
    VECTOR_INT types;
    std::vector<CONSTRAINT> list;
    //-------------
    // END PARAMETERS
    //-------------

    //-------------
    // CONSTRUCTORS
    //-------------
    CONSTRAINTS(){};

    CONSTRAINTS(int constraints_n_constr, VECTOR_INT constraints_types, std::vector<VECTOR> constraints_parameters)
    {
        if (constraints_n_constr < 1)
        {
            throw_line("ERROR: invalid number of constraints. Less than 1\n");
        }
        n_constr = constraints_n_constr;
        save_constraints_type(constraints_types);
        build_constraints_list(constraints_parameters);
    }

    void initialize(int constraints_n_constr, VECTOR_INT constraints_types, std::vector<VECTOR> constraints_parameters);

    void save_constraints_type(VECTOR_INT constraints_types);

    void build_constraints_list(std::vector<VECTOR> constraints_parameters);
    //-----------------
    // END CONSTRUCTORS
    //-----------------
};
//--------------------------------------------
// END CONSTRAINTS HANDLER CLASS
//--------------------------------------------








