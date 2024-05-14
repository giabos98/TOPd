#pragma once

#include "../CODE_HEADERS/codeHeader.h"
#include "../geometry.h"

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
    prec Vr    = -1.0;
    prec vol_0 = -1.0;
    prec vol   = -1.0;
    // subdomain volume constraint
    int domain_id = -1;

    //bound size constraint
    int bound_id = -1;
    prec Sr      = -1.0;
    prec surf_0  = -1.0;
    prec surf    = -1.0;

    //discretizing constraint
    prec discretization_toll = -1.0;
    prec discretization_res  = -1.0;

    // WSS constraint
    prec critical_WSS = -1.0;
    int sign = 1; // sign=-1: <, sign=1: >
    prec actual_WSS = -1.0;
    prec value_WSS = 0.0;

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

    void initialize_subdomain_volume_constraint(int constr_type, int subdomain, prec volume_fraction);

    void initialize_bound_constraint(int constr_type, int bound, prec bound_fraction);

    void initialize_discretizing_constraint(int constr_type, prec discretization_tollerance);

    void initialize_WSS_constraint(int constr_type, prec crit_WSS, int sign_value);
    //-----------------
    // END CONSTRUCTORS
    //-----------------
};
//--------------------------------------------
// END BASE CONSTRAINT CLASS
//--------------------------------------------

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
    int max_constraint_type = 4;
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


//-----------------------------------------------------
// WSS CONSTRAINT EVALUATION AND DERIVATION FACILITATOR
//-----------------------------------------------------

class WSS_CONSTRAINT
{
    public:

    //--------------
    // PARAMETERS
    //--------------
    PHYSICS* physics;
    std::vector<MATRIX> velocity_gradient;
    std::vector<VECTOR> gamma_gradient;
    int n_tan = 0;
    std::vector<std::vector<VECTOR>> gamma_tangents;
    std::vector<VECTOR> nodal_WSS_tangent_values;
    VECTOR nodal_WSS_values;

    //--------------
    // METHODS
    //--------------
    WSS_CONSTRAINT(PHYSICS &physics_ptr, std::vector<MATRIX> velocity_gradient_in, std::vector<VECTOR> gamma_gradient_in, std::vector<std::vector<VECTOR>> gamma_tangents_in)
    {
        initialize(physics_ptr, velocity_gradient_in, gamma_gradient_in, gamma_tangents_in);
    }

    void initialize(PHYSICS &physics_ptr, std::vector<MATRIX> velocity_gradient_in, std::vector<VECTOR> gamma_gradient_in, std::vector<std::vector<VECTOR>> gamma_tangents_in);
    void eval_nodal_values();


};
