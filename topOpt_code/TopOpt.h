#pragma once

#include "CODE_HEADERS/codeHeader.h"
#include "OPTIMIZATION/Optimizer.h"
#include "problemNSHeader.h"
#include "geometry.h"
#include "UTILITIES/Optimization/opt_utils.h"
#include "UTILITIES/mesh_scraper/scrape_meshes.h"

// #include <sys/time.h>

class TOP_OPT
{
public:
    //----------------------------
    //-----------------
    // CLASS PARAMETERS
    //-----------------
    //----------------------------
    std::string name;
    int optimization_scheme;
    std::string inputFile;
    TIME_PROFILER time_profiler;
    MESH_SCRAPER mesh_scraper;
    int use_MMG_optimization;
    prec mmg_level_set = 0.5;
    int custom_opt_mesh_sizes = 0;
    prec mmg_hmin;
    prec mmg_hmax;
    int perturb_solution;
    bool already_perturbed = false;
    int perturb_maxIt = 0;

    int funcId = 1;
    bool binPrint;
    int deltaPrint;
    int minIt;
    int maxIt;
    int complete_maxIt;
    prec change_toll;

    prec V0;
    int customFunc;
    int nFunctionals = 4;
    int onlyGrad;
    VECTOR beta;
    prec Vr;
    prec Vol;
    VECTOR alpha;
    int smooth_gamma_between_element = 0;

    //--- PHYSICS PARAMETERS ---
    PHYSICS physics;
    int n_subdomains;
    int n_opt_subdomains;
    VECTOR_INT subdomains;
    VECTOR_INT opt_subdomains;
    VECTOR_INT not_opt_nodes_in_subdomains;
    VECTOR_INT is_elem_in_dom;
    VECTOR_INT is_node_in_dom;
    VECTOR subdomains_initial_values;
    int nNodeInDom;
    int nElemInDom;
    VECTOR_INT optNodeFromGlobNode;
    VECTOR_INT nodeInDom;
    VECTOR_INT elemInDom;

    //--- EQUATIONS ---
    PROBLEM_NS NS;
    ADJOINT_NS ADJ;
    OPTIMIZER Optimizer;

    //--- DIFFUSION FILTER FOR THE TOP. OPT. PARAMETER
    int flagDomain;
    MATRIX optBox;
    prec max_opt_box_side;
    int enableDiffusionFilter;
    prec diffusionRadiusPercentage;
    prec diffusionFilterWeight;

    // INITIAL TOP. OPT. PARAMETER VALUE
    prec gamma_initial_value;

    //--- SOLUTIONS ---
    VECTOR lastSolNS;
    VECTOR lastSolADJ;

    //--- FLUID ENERGY ---
    VECTOR temp_fluid_energy;

    //------------------
    // PRINT
    //------------------
    VTK VTKWriter;
    int write_on_velocity_mesh = 0;
    int write_inital_condition = 0;

    //------------------
    // STATISTICS
    //------------------
    prec import_time;
    prec solution_time;

    //----------------------------
    //--------------
    // CLASS METHODS
    //--------------
    //----------------------------
    TOP_OPT(std::string InputFile);
    //-----------------------
    // PROBLEM INIZIALIZATION
    //-----------------------
    void importParameters(std::string inputFile);
    //-------
    // SOLVER
    //-------
    void prepareNS();
    void solveNS();
    void prepareADJ();
    void solveADJ();
    void handle_optimization_domain();
    void handle_gamma_initial_condition(VECTOR &gamma_opt, VECTOR &gamma);
    void save_gammaOpt_in_gammaFull(VECTOR &gamma_opt, VECTOR &gamma_full);

    void solve();

    //-----------------------------------
    // EXPORT OPTIMIZED GEOMETRY NODES
    //-----------------------------------
    void exportOptimizedDomain(VECTOR gamma, prec gammaMin, MATRIX_INT &optElem);

    void print_stats(prec totalTime);

    void eval_gamma_gradiend(VECTOR &gamma, MATRIX &grad_gamma, VECTOR &grad_gamma_norm);
    void eval_gamma_gradient_norm_with_filter(VECTOR &gamma, VECTOR &grad_gamma_norm);
    //-----------------------------------
    // PREPROCESS QUANTITIES TO PRINT
    //-----------------------------------
    void prepare_solution_print(MATRIX &U_print, VECTOR &P, VECTOR &U_magnitude, VECTOR P_print);
    void evaluate_U_magnitude(MATRIX &U_print, VECTOR &U_magnitude);
    void get_pressure_in_nodes_v(VECTOR &P, VECTOR &P_print);
    
    void print_optimization_results(int &nNodes_v, int &dim, int &nNodes, int &loop, int &currLoopPrint, prec &obj, prec &change, VECTOR &gamma, bool &feasible, MATRIX &funcValues, MATRIX &no_weights_funcValues, VECTOR &changes, VECTOR_INT &valid, VECTOR &grad_gamma_norm);
    
    void print_results_in_console(int &loop, prec &obj, prec &change);

    void print_results_in_vtk(int &nNodes_v, int &dim, int &nNodes, int &loop, int &currLoopPrint, VECTOR &gamma, VECTOR &grad_gamma_norm);
    
    void print_for_matlab_interface(int &loop, prec &obj, prec &change, bool &feasible, MATRIX &funcValues, MATRIX &no_weights_funcValues, VECTOR &changes, VECTOR_INT &valid);
    
    void print_mesh_for_postprocessing(VECTOR &gamma);

    void export_optimized_domain_mesh_with_mmg();

    void evaluate_total_energy();

    void perturb_gamma_solution(VECTOR &gamma);

};