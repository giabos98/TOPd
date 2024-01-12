#pragma once

#include "CODE_HEADERS/codeHeader.h"
#include "OPTIMIZATION/Optimizer.h"
#include "problemNSHeader.h"
#include "geometry.h"
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

    int funcId = 1;
    bool binPrint;
    int flagPrint;
    int deltaPrint;
    int minIt;
    int maxIt;
    int change_toll;

    prec alpha_min;
    prec alpha_max;
    prec q;
    prec V0;
    int customFunc;
    int nFunctionals = 4;
    int onlyGrad;
    VECTOR beta;
    prec Vr;
    prec Vol;
    VECTOR alpha;

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
    int write_on_velocity_mesh;

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

    //-----------------------------------
    // PREPROCESS QUANTITIES TO PRINT
    //-----------------------------------
    void prepare_solution_print(MATRIX &U_print, VECTOR &P, VECTOR &U_magnitude, VECTOR P_print);
    void evaluate_U_magnitude(MATRIX &U_print, VECTOR &U_magnitude);
    void get_pressure_in_nodes_v(VECTOR &P, VECTOR &P_print);
    
    void print_optimization_results(int &nNodes_v, int &dim, int &nNodes, int &loop, int &currLoopPrint, prec &obj, prec &change, bool &printNSSol, VECTOR &gamma, bool &feasible, MATRIX &funcValues, MATRIX &no_weights_funcValues, VECTOR &changes, VECTOR_INT &valid);
    
    void print_results_in_console(int &loop, prec &obj, prec &change);

    void print_results_in_vtk(int &nNodes_v, int &dim, int &nNodes, int &loop, int &currLoopPrint, bool &printNSSol, VECTOR &gamma);
    
    void print_for_matlab_interface(int &loop, prec &obj, prec &change, bool &feasible, MATRIX &funcValues, MATRIX &no_weights_funcValues, VECTOR &changes, VECTOR_INT &valid);
    
    void evaluate_total_energy();

};