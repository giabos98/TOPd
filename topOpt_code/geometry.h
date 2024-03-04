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

    // Top. Opt. parameter filtering geometry

    PHYSICS(){};

    void initialize();

    void print();

    void eval_solution_times();
    void build_bounds_elems_v();
    void build_bounds_elems_surfaces_v();

    // static prec get_surface(MATRIX &matCoord, int dim);
    static prec get_surface(MATRIX &matCoord, int dim)
    {
        if (dim != 2 && dim != 3) throw_line("ERROR: Getting area in a space different from 2D or 3D");
        if (matCoord.nRow != dim || matCoord.nCol != dim) throw_line("ERROR: Getting area of coords incompatible with dim\n");
        prec area;
        switch (dim)
        {
            case 2:
            {
                std::shared_ptr<prec[]>  vec(new prec[dim]);
                vec[0] = matCoord[0][0] - matCoord[1][0];
                vec[1] = matCoord[0][1] - matCoord[1][1];
                area = VECTOR::norm(vec,2);
                break;
            }
            case 3:
            {
                VECTOR vec31(dim);
                vec31[0] = matCoord[2][0] - matCoord[0][0];
                vec31[1] = matCoord[2][1] - matCoord[0][1];
                vec31[2] = matCoord[2][2] - matCoord[0][2];
                VECTOR vec21(dim);
                vec21[0] = matCoord[1][0] - matCoord[0][0];
                vec21[1] = matCoord[1][1] - matCoord[0][1];
                vec21[2] = matCoord[1][2] - matCoord[0][2];
                // VECTOR tempVec; tempVec = vec31.cross(vec21);
                // vec21 = tempVec;
                VECTOR newVec;
                VECTOR::cross(vec31, vec21, newVec);
                
                area = VECTOR::norm(newVec)/2;
                break;
            }
        }
        return area;
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
    prec Vr = -1.0;
    prec vol_0 = -1.0;
    prec vol = -1.0;
    // subdomain volume constraint
    int domain_id = -1;

    //bound size constraint
    int bound_id = -1;
    prec Sr = -1.0;
    prec surf_0 = -1.0;
    prec surf = -1.0;
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
    int max_constraint_type = 2;
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








