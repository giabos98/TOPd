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

    prec alpha_min;
    prec alpha_max;
    int use_smooth_alpha_max = 0;
    MATRIX smooth_alpha_mat;
    int n_alpha_max_steps = 0;
    VECTOR real_alpha_max_history;
    VECTOR real_alpha_max_vec;
    prec real_alpha_max;
    VECTOR_INT alpha_it;
    prec q;
    VECTOR alpha;
     
    //--- SOLUTION INFO ---
    prec curr_opt_it = 0;

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
    MATRIX coord;
    MATRIX_INT elem;
    MATRIX Bloc; MATRIX Cloc; MATRIX Dloc;
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
    prec target_vol_fract;

    // Top. Opt. parameter filtering geometry

    PHYSICS(){};

    void initialize();
    void parse_bounds();

    void set_alpha_min(VECTOR &alpha_min_vector);
    void set_alpha_max();
    void set_alpha_max_staircase();
    void update_real_alpha_max();
    void update_alpha(VECTOR &gamma, VECTOR &alpha);

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
    void print_mesh_for_mmg(std::string file_name, std::string rel_path = "");

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










