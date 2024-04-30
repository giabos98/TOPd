#pragma once

#include "CODE_HEADERS/codeHeader.h"

class PHYSICS
{
    public:

    //--- NS PROBLEM NAME
    std::string name;

    //--- DIMENSION ---
    int dim;
    int isNS; // physics model (0:Stokes, 1:Navier-Stokes)

    //--- ABOUT TIME ---
    prec t_end;
    prec deltaT;
    prec deltaT_min;
    prec convergence_scale_factor;

    //--- CONSOLE LOG ---
    int completeLog;

    //--- PHYSICS PARAMETERS ---
    prec half_domain_thickness = 0.0;
    prec rho;
    prec mu;
    prec ni;
    VECTOR domains_permeability;

    //---  SOLVER FLAGS ---
    bool isStationary = true;
    bool flagBC; //0 penalty, 1 lifting 

    //--- V Struct -------
    int nNodes_v;
    int nElem_v;
    VECTOR h_v;
    MATRIX coord_v;
    MATRIX_INT elem_v;
    std::vector<VECTOR> centroids_v;
    VECTOR_INT elem_geo_entities_ids_v; // note: the domains ids are scaled with c++ notation (starting from 0)
    std::vector<VECTOR_INT> elems_in_doms;
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
    VECTOR h;
    MATRIX coord;
    MATRIX_INT elem;
    std::vector<VECTOR> centroids;
    VECTOR_INT elem_geo_entities_ids;
    MATRIX Bloc; MATRIX Cloc; MATRIX Dloc;
    std::vector<MATRIX> Coef;
    VECTOR Volume;

    //--- BOUNDS INFO ---
    std::string bound_nodes_v_path;
    std::string bound_elems_v_path;
    std::vector<VECTOR_INT> bound_nodes_v;
    std::vector<MATRIX_INT> bound_elems_v;

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
    std::vector<MATRIX_INT> bounds_elems_v;
    std::vector<VECTOR> bounds_elems_surface_v;

    // ABOUT SOLUTION TIMES
    VECTOR solution_times;
    VECTOR solution_deltaT;

    // PDE solution matrices for the current iteration
    MATRIX NS_solution;
    MATRIX ADJ_solution;
    MATRIX Darcy_solution;

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

    void initialize();
    void parse_bounds();

    void print();

    void eval_solution_times();
    void build_bounds_elems_v();
    void build_bounds_elems_surfaces_v();
    void build_centroids_v();

    void invert_sizes(std::vector<VECTOR> &vec, std::vector<VECTOR> &inverted_vec);

    void decompose_NS_solution(VECTOR &sol, MATRIX &U_sol, VECTOR &P_sol, int transpose_U = 0);
    void decompose_NS_solution_over_time(MATRIX &sol, std::vector<MATRIX> &U_sol, std::vector<VECTOR> &P_sol);

    void eval_mean_solution_over_time(MATRIX &solution_over_times, VECTOR &time_avg_sol);

    // static prec get_surface(MATRIX &matCoord, int dim);
    void eval_on_elements_v(VECTOR &value, VECTOR &element_value);
    void eval_on_elements_v(VECTOR &value, VECTOR_INT &nodesFromNodes_v, VECTOR_INT &elems, VECTOR &element_value);
    void eval_on_centroids_v(VECTOR &value, VECTOR &centroids_value);
    void smooth_between_elements(VECTOR &value, VECTOR &smoothed_value);
    void smooth_between_elements(VECTOR &value, VECTOR_INT &nodesFromNodes_v, VECTOR_INT &elems, VECTOR &smoothed_value);
    static prec get_surface(MATRIX &matCoord, int dim);
    void eval_gradient(VECTOR &value, std::vector<VECTOR> &gradient, int eval_method = 0); // eval method means: gradient L2 projection method
    void eval_gradient(MATRIX &value, std::vector<MATRIX> &gradient, int eval_method = 0);// eval method means: gradient L2 projection method
    void eval_divergence(VECTOR &value, VECTOR &divergence);
    void eval_divergence(MATRIX &value, std::vector<VECTOR> &divergence);
    void eval_gradient_norm(std::vector<VECTOR> &gradient, VECTOR &norm);
    void eval_gradient_norm(std::vector<MATRIX> &gradient, VECTOR &norm);
    void eval_WSS(MATRIX &value, VECTOR_INT &nodes, std::vector<VECTOR> &normals, std::vector<VECTOR> &WSS);
    void eval_WSS(MATRIX &value, int bound_id, VECTOR normal, std::vector<VECTOR> &WSS);
    void eval_WSS_avg_over_time(std::vector<MATRIX> &value, VECTOR_INT &nodes, std::vector<VECTOR> &normals, std::vector<VECTOR> &WSS);
    void eval_directional_gradient(VECTOR &value, VECTOR_INT &nodes, std::vector<VECTOR> &directions, VECTOR &dir_gradient);
    void eval_directional_gradient(MATRIX &value, VECTOR_INT &nodes, std::vector<VECTOR> &directions, std::vector<VECTOR> &dir_gradient);
    void build_tangents_from_normal(VECTOR &normal, std::vector<VECTOR> &tangent_vectors);
    void build_tangents_from_normals(std::vector<VECTOR> &normals, std::vector<std::vector<VECTOR>> &tangent_vectors);

    static int factorial(int n)
    {
        if (n == 0)
        {
            return 1;
        }
        else
        {
            return n * factorial(n-1);
        }
        throw_line("ERROR: not handled case in factorial function\n");
    }

    static VECTOR_INT match_indices(VECTOR_INT index)
    {
        int max = index.max();
        VECTOR_INT counter;
        counter.setZeros(max+1);
        for (int id = 0; id < index.length; id++)
        {
            counter[index[id]] += 1;
        }
        return counter;
    }

    static int eval_elemental_integral_factor_by_indices(VECTOR_INT index)
    {
        int factor = 1;
        VECTOR_INT counter = match_indices(index);
        for (int id = 0; id < counter.length; id++)
        {
            factor *= factorial(counter[id]);
        }
        return factor;
    }
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








