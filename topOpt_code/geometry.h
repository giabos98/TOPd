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
    

    PHYSICS(){}

    void print()
    {
        printf(" dim: %d\n rho: %" format "\n: mu: %" format" \n ni: %" format "\n nNodes_v:%d\n", dim, rho, mu, ni, nNodes_v);
    }

    void eval_solution_times()
    {
        switch (isStationary)
        {
            case 1: // stationary solution
            {
                t_end = 1e16;
                deltaT = t_end;
                deltaT_min = 1e-5;
                solution_times.append(t_end);
                solution_deltaT.append(deltaT);
                convergence_scale_factor = 0.5;
                break;
            }
            default: // time dependent solution
            {
                solution_times.resetZeros();
                prec curr_time = 0.0;
                if (abs(t_end-curr_time) <= 1e-10) throw_line("ERROR: too small ending time (<= 1e-10)\n");
                while (abs(t_end-curr_time) > 1e-10)
                {
                    solution_times.append(curr_time);
                    curr_time += deltaT;
                }
                solution_times.append(t_end);
                solution_deltaT.initialize(solution_times.length-1);
                solution_deltaT.reset(deltaT);
                solution_deltaT[solution_deltaT.length-1] = t_end - solution_times[solution_times.length-2]; //adjust the last time step
                break;
            }
        }
    }
};





