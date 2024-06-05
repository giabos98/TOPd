#include "TopOpt.h"

namespace fs = std::filesystem;

TOP_OPT::TOP_OPT(std::string InputFile)
{
    prec startTime = omp_get_wtime();
    inputFile = InputFile;
    time_profiler.initialize();
    printf("\n-------\n--| INITIALIZING NS PROBLEM |--\n-------\n");

    PHYSICS* tempP = &physics;
    NS.initialize(tempP, inputFile, physics.alpha, false);
    alpha.length = physics.alpha.length;
    alpha.P = physics.alpha.P;

    printf("\n-------\n--| INITIALIZING PHYSICS |-- \n-------\n");
    physics.initialize();

    printf("\n-------\n--| INITIALIZING ADJOINT PROBLEM |-- \n-------\n");
    ADJ.initialize(NS, physics.alpha);


    // initialize TopOpt parameters
    printf("\n-------\n-| INITIALIZING TOP OPT PROBLEM |--\n-------\n");
    importParameters("INPUT_FILES/TopOptInput.txt");
    if (enableDiffusionFilter > 0) 
    {
        Optimizer.diffusionFilter.initialize(enableDiffusionFilter, tempP, nNodeInDom, nodeInDom, optNodeFromGlobNode ,optBox, diffusionRadiusPercentage, Optimizer.diffusion_filter_case);
        //diffusionFilter.printNeighbourhood();
    }

    ADJ.constraints = &(Optimizer.constraints);

    prec endTime = omp_get_wtime();
    import_time = endTime - startTime;
}

void TOP_OPT::importParameters(std::string inputFile)
{
    std::ifstream ParameterFile;
    ParameterFile.open(inputFile);

    MATRIX boxOpt(3,2);
    std::string line;
    std::istringstream iss;

    // PROBLEM NAME
    STREAM::getLines(ParameterFile, line, 2);
    STREAM::getValue(ParameterFile, line, iss, name);
    
    // OPTIMIZATION SCHEME
    STREAM::getLines(ParameterFile, line, 2);
    STREAM::getValue(ParameterFile, line, iss, optimization_scheme);
    if (optimization_scheme == 2)
    {
        throw_line("ERROR: GCMMA currently not working\n");
    }
    else if (optimization_scheme > 2)
    {
        throw_line("ERROR: non valid optimization scheme id\n");
    }
    
    // SUBDOMAINS
    STREAM::getLines(ParameterFile, line, 3);
    STREAM::getValue(ParameterFile, line, iss, n_subdomains);
    subdomains.initialize(n_subdomains);
    for (int idom = 0; idom < n_subdomains; idom++)
    {
        subdomains[idom] = idom;
    }
    STREAM::getLines(ParameterFile, line, 1);
    subdomains_initial_values.initialize(n_subdomains);
    STREAM::getRowVector(ParameterFile, line, iss, subdomains_initial_values);
    STREAM::getLines(ParameterFile, line, 1);
    STREAM::getValue(ParameterFile, line, iss, n_opt_subdomains);
    STREAM::getLines(ParameterFile, line, 1);
    opt_subdomains.initialize(n_opt_subdomains);
    STREAM::getRowVector(ParameterFile, line, iss, opt_subdomains);
    handle_optimization_domain();

    // BRINKMAN PENALIZATION
    STREAM::getLines(ParameterFile, line, 3);
    VECTOR alpha_min_vector(2);
    STREAM::getRowVector(ParameterFile, line, iss, alpha_min_vector);
    physics.set_alpha_min(alpha_min_vector);
    STREAM::getLines(ParameterFile, line, 1);
    STREAM::getValue(ParameterFile, line, iss, physics.alpha_max);
    STREAM::getLines(ParameterFile, line, 1);
    STREAM::getValue(ParameterFile, line, iss, physics.q);
    STREAM::getLines(ParameterFile, line, 1);
    STREAM::getValue(ParameterFile, line, iss, physics.alpha_it);

    // CONSTRAINTS
    CONSTRAINTS constraints;
    STREAM::getLines(ParameterFile, line, 3);
    int n_constr;
    STREAM::getValue(ParameterFile, line, iss, n_constr);
    STREAM::getLines(ParameterFile, line, 6);
    VECTOR_INT constraints_types_list(n_constr);
    STREAM::getRowVector(ParameterFile, line, iss, constraints_types_list);
    std::vector<VECTOR> constraints_parameters(n_constr);
    STREAM::getLines(ParameterFile, line, 1);
    for (int icons = 0; icons < n_constr; icons++)
    {
        constraints_parameters[icons].initialize(constraints.max_parameters_length);
        STREAM::getRowVector(ParameterFile, line, iss, constraints_parameters[icons]);
    }
    constraints.initialize(n_constr, constraints_types_list, constraints_parameters);

    // get total volume
    VECTOR Volume = physics.Volume_v;
    V0 = 0;
    // elemInDom.printRow("elDom");
    for (int iel = 0; iel < nElemInDom; iel++)
    {
        int globElem = elemInDom[iel];
        V0 += Volume[globElem];
    } 
    physics.V0 = V0;
    physics.Vol = V0;
    physics.vol_fract = 1.0;
    Vr = constraints.list[0].Vr;
    // std::cout << "V0: " << V0 << "\n";
    
    // USELESS BASE FUNCTINONALS
    STREAM::getLines(ParameterFile, line, 5);
    STREAM::getValue(ParameterFile, line, iss, onlyGrad);

    // ITERATIONS
    STREAM::getLines(ParameterFile, line, 3);
    STREAM::getValue(ParameterFile, line, iss, minIt);
    STREAM::getLines(ParameterFile, line, 1);
    STREAM::getValue(ParameterFile, line, iss, maxIt);

    // FUNCTIONAL DEFINITION
    STREAM::getLines(ParameterFile, line, 2);
    STREAM::getValue(ParameterFile, line, iss, customFunc);
    STREAM::getLines(ParameterFile, line, 2);
    int time_integration;
    STREAM::getValue(ParameterFile, line, iss, time_integration);
    STREAM::getLines(ParameterFile, line, 2);
    beta.initialize(nFunctionals);
    STREAM::getRowVector(ParameterFile, line, iss, beta);

    STREAM::getLines(ParameterFile, line, 2);
    STREAM::getValue(ParameterFile, line, iss, change_toll);

    STREAM::getLines(ParameterFile, line, 2);
    int opt_acceleration_case;
    STREAM::getValue(ParameterFile, line, iss, opt_acceleration_case);
    if ((opt_acceleration_case < 0) || (opt_acceleration_case > 2)) throw_line("ERROR: not handled gamma acceleration case\n");
    getline(ParameterFile, line);
    prec beta_MAX;
    STREAM::getValue(ParameterFile, line, iss, beta_MAX);
    getline(ParameterFile, line);
    prec beta_min;
    STREAM::getValue(ParameterFile, line, iss, beta_min);
    getline(ParameterFile, line);
    int beta_interpolation;
    STREAM::getValue(ParameterFile, line, iss, beta_interpolation);
    getline(ParameterFile, line);
    prec change_max;
    STREAM::getValue(ParameterFile, line, iss, change_max);
    getline(ParameterFile, line);
    prec change_min;
    STREAM::getValue(ParameterFile, line, iss, change_min);
    getline(ParameterFile, line);
    prec crit_change;
    STREAM::getValue(ParameterFile, line, iss, crit_change);
    getline(ParameterFile, line);
    prec crit_beta;
    STREAM::getValue(ParameterFile, line, iss, crit_beta);

    // check gamma projection settings
    if (opt_acceleration_case == 2)
    {
        if (beta_min <= 0)
        {
            throw_line("ERROR: beta_min has non valid value.\n");
        }
        if (beta_interpolation == 1)
        {
            if (beta_MAX <= beta_min)
            {
                throw_line("ERROR: beta_max <= beta_min.\n");
            }
            if (crit_beta <= beta_min)
            {
                throw_line("ERROR: crit_beta <= beta_min.\n");
            }
            else if (crit_beta >= beta_MAX)
            {
                throw_line("ERROR: crit_beta >= beta_max.\n");
            }
            if (change_min < change_toll)
            {
                throw_line("ERROR: change_min < change_toll.\n");
            }
            if (change_max >= 1.0)
            {
                throw_line("ERROR: change_max >= 1.\n");
            }
            if (crit_change <= change_min)
            {
                throw_line("ERROR: crit_change <= change_min.\n");
            }
            else if (crit_change >= change_max)
            {
                throw_line("ERROR: crit_change >= change_max.\n");
            }
        }
    }
    
    STREAM::getLines(ParameterFile, line, 2);
    STREAM::getValue(ParameterFile, line, iss, enableDiffusionFilter);
    getline(ParameterFile, line);
    STREAM::getValue(ParameterFile, line, iss, diffusionRadiusPercentage);
    getline(ParameterFile, line);
    STREAM::getValue(ParameterFile, line, iss, diffusionFilterWeight);

    STREAM::getLines(ParameterFile, line, 2);
    STREAM::getValue(ParameterFile, line, iss, smooth_gamma_between_element);
    if (customFunc == 1)
    {
        Optimizer.initialize(&physics, nodeInDom, elemInDom, optNodeFromGlobNode, V0, Vr, customFunc, time_integration, onlyGrad, beta, opt_acceleration_case, beta_MAX, beta_min, beta_interpolation, change_max, change_min, crit_change, crit_beta, enableDiffusionFilter, constraints, smooth_gamma_between_element);
    }
    else
    {
        throw_line("\nERROR: customFunc != 1\n");
    }
    
    //------------------------------
    // READ PRINT INFO
    //------------------------------
    STREAM::getLines(ParameterFile, line, 3);
    STREAM::getValue(ParameterFile, line, iss, binPrint);

    //
    STREAM::getLines(ParameterFile, line, 2);
    STREAM::getValue(ParameterFile, line, iss, deltaPrint);

    //
    STREAM::getLines(ParameterFile, line, 1);
    STREAM::getValue(ParameterFile, line, iss, write_on_velocity_mesh);

    //
    STREAM::getLines(ParameterFile, line, 1);
    STREAM::getValue(ParameterFile, line, iss, write_inital_condition);


    //------------------------------
    // READ OPENMP INFO
    //------------------------------
    STREAM::getLines(ParameterFile, line, 3);
    int thread_case;
    STREAM::getValue(ParameterFile, line, iss, thread_case);
    switch (thread_case)
    {
        case -1:
        {
            PARALLEL::nThread = std::thread::hardware_concurrency();
            break;
        }
        case 0:
        {
            break;
        }
        default:
        {
            PARALLEL::nThread = thread_case;
            break;
        }
    }

    // CLOSE STREAMING
    ParameterFile.close();

    std::string folderName = NS.name + "/" + name;
    if (write_on_velocity_mesh == 0)
    {
        VTKWriter.initializeForTopOpt(folderName, (physics).dim, (physics).nNodes, (physics).nElem, deltaPrint, maxIt, binPrint);
    }
    else if (write_on_velocity_mesh == 1)
    {
        VTKWriter.initializeForTopOpt(folderName, (physics).dim, (physics).nNodes_v, (physics).nElem_v, deltaPrint, maxIt, binPrint);
    }
    else
    {
        throw_line("ERROR: invalid mesh option for print\n");
    }

    fs::copy_file("INPUT_FILES/readProblemNS.txt", "results/" + folderName + "/readProblemNS.txt", fs::copy_options::overwrite_existing);
    fs::copy_file("INPUT_FILES/TopOptInput.txt", "results/" + folderName + "/TopOptInput.txt", fs::copy_options::overwrite_existing);
    
    NS.VTKWriter.binWrite = binPrint;
}

//-----------------------
// NS SOLVER
//-----------------------
void TOP_OPT::prepareNS() // prepare NS solver
{
    printf("\n-------\n-| PREPARE NS |--\n-------\n");
    NS.setBC();
    NS.prepareSolver();
    physics.NS_solution.complete_reset();
    physics.NS_solution.initialize(physics.solution_times.length, physics.nDof);
}
//---
void TOP_OPT::solveNS() // NS solver
{
    printf("\n-------\n-| NAVIER-STOKES SOLVER |--\n-------\n\n");
    if (physics.isStationary == 1)
    {
        NS.StatSolverIterative();
    } 
    else 
    {
        NS.Solver();
    }
    lastSolNS = physics.NS_solution.get_row(physics.solution_times.length-1);
}

//-----------------------
// ADJ SOLVER
//-----------------------
void TOP_OPT::prepareADJ() // prepare ADJ solver
{
    printf("\n-------\n-| PREPARE ADJOINT |--\n-------\n");
    
    if (customFunc == 1)
    {
        ADJ.fWeights.initialize(nFunctionals);
        ADJ.fWeights = beta;
    } 
    ADJ.customFunc =  customFunc;
    ADJ.onlyGrad = onlyGrad;
    ADJ.setBC(NS);
    printf("\n-------\n-| ADJOINT PREPRO |--\n-------\n");
    ADJ.prepareSolver(NS);
    physics.ADJ_solution.complete_reset();
    physics.ADJ_solution.initialize(physics.solution_times.length, physics.nDof);
}
//---
void TOP_OPT::solveADJ() // ADJ solver
{
    if (NS.completeLog < 2) printf("\n-------\n-| ADJOINT STAT FORCING |--\n-------\n");
    // MATRIX U = NS.getVelocityFromSol(lastSolNS);
    // ADJ.globIter = NS.globIter;
    ADJ.applyStatForcing();

    printf("\n-------\n-| ADJOINT SOLVER |--\n-------\n\n");

    if (physics.isStationary == 1)
    {
        ADJ.StatSolver(lastSolNS);
    } 
    else 
    {
        ADJ.Solver();
    }
    lastSolADJ = ADJ.lastSol;
}

// save the gamma value in the optimization nodes in the whole gamma vector
void TOP_OPT::save_gammaOpt_in_gammaFull(VECTOR &gamma_opt, VECTOR &gamma_full)
{
    if (gamma_opt.length != nodeInDom.length) throw_line("ERROR: wrong size of gamma_opt in copying procedure.\n");
    if (gamma_full.length != physics.nNodes_v) throw_line("ERROR: wrong size of gamma_full in copying procedure.\n");
    for (int inode = 0; inode < nodeInDom.length; inode++)
        {
            int iglob = nodeInDom[inode];
            gamma_full[iglob] = gamma_opt[inode]; 
        }
}

//-----------------------------------
// EXPORT OPTIMIZAED GEOMETRY NODES
//-----------------------------------
void TOP_OPT::exportOptimizedDomain(VECTOR gamma, prec gammaMin, MATRIX_INT &optElem)
{
    int dim = physics.dim;
    //int nNodes_v = physics.nNodes_v;
    int nElem_v = physics.nElem_v;
    // VECTOR optNodes(nNodes_v);
    MATRIX_INT elem_v = physics.elem_v;
    //MATRIX coord_v = physics.coord_v;
    //int nNodesOpt = 0;
    int nElOpt = 0;
    for (int iel = 0; iel < nElem_v; iel++)
    {
        int* tempElem_v = elem_v[iel];
        int isOptEl = 0;
        for (int inod = 0; inod < dim+1; inod++)
        {
            int iglob = tempElem_v[inod];
            if (gamma[iglob] < gammaMin)
            {
                isOptEl++;
            }
        }
        if (isOptEl < 1)
        {
            int* tempOptElem = optElem[nElOpt];
            for (int inod = 0; inod < dim+1; inod++)
            {
                tempOptElem[inod] = tempElem_v[inod];
            }
            nElOpt++;
        }
    }
    optElem.shrinkRows(nElOpt);
}
//-----------------------------------

void TOP_OPT::eval_gamma_gradiend(VECTOR &gamma, MATRIX &grad_gamma, VECTOR &grad_gamma_norm)
{
    grad_gamma.resetZeros();
    grad_gamma_norm.reset(1.0);
    int dim = physics.dim;
    int nNodes_v = physics.nNodes_v;
    int nElem_v = physics.nElem_v; 
    for (int iel = 0; iel < nElem_v; iel++)
    {
        for (int iloc = 0; iloc < (dim+1); iloc++)
        {
            int iglob = physics.elem_v[iel][iloc];
            prec temp_gamma = gamma[iglob];
            for (int icomp = 0; icomp < dim; icomp++)
            {
                prec temp_grad = physics.Coef_v[icomp][iel][iloc] * temp_gamma;
                grad_gamma[iglob][icomp] += temp_grad;
            }
        }
    }

    for (int inod = 0; inod < nNodes_v; inod++)
    {
        prec temp_value = 0;
        for (int icomp = 0; icomp < dim; icomp++)
        {
            temp_value += grad_gamma[inod][icomp]*grad_gamma[inod][icomp];
        }
        grad_gamma_norm[inod] = sqrt(temp_value);
    }
}
void TOP_OPT::eval_gamma_gradient_norm_with_filter(VECTOR &gamma, VECTOR &grad_gamma_norm)
{
    VECTOR filtered_gamma(nNodeInDom);
    Optimizer.diffusionFilter.filter_gamma(gamma, filtered_gamma, 1);
    for (int inod = 0; inod < physics.nNodes_v; inod++)
    {
        grad_gamma_norm = abs(gamma[inod] - filtered_gamma[inod]) / Optimizer.diffusionFilter.diffRadius;
    }
}


//---------------------------------------
// Print Topology Optimization Code Stats
//---------------------------------------
void TOP_OPT::print_stats(prec totalTime)
{
    std::cout << "\n   ----------------------------------------------------";
    std::cout << "\n   | TOPOLOGY OPTIMIZATION PROBLEM EXECUTED CORRECTLY |";
    std::cout << "\n   |---> Name: " << NS.name + " / " + name;
    std::cout << "\n   |---> Import Time: " << import_time;
    std::cout << "\n   |---> Solution Time: " << solution_time;
    std::cout << "\n   |---> Total Time: " << totalTime;
    std::cout << "\n   ----------------------------------------------------\n\n";
}

//-----------------------------------
// PREPROCESS QUANTITIES TO PRINT
//-----------------------------------
void TOP_OPT::prepare_solution_print(MATRIX &U_print, VECTOR &P, VECTOR &U_magnitude, VECTOR P_print)
{
    evaluate_U_magnitude(U_print, U_magnitude);
    get_pressure_in_nodes_v(P, P_print);
}

void TOP_OPT::evaluate_U_magnitude(MATRIX &U_print, VECTOR &U_magnitude)
{
    int n_nodes = U_print.nRow;
    if (U_magnitude.length != n_nodes) throw_line("ERROR: trying to save magnitude values in a vector with incompatible size with the data matrix\n");
    int dim = U_print.nCol;
    for (int inod = 0; inod < n_nodes; inod++)
    {
        prec* curr_U = U_print[inod];
        std::shared_ptr<prec[]> temp_U = 0;
        temp_U = std::shared_ptr<prec[]>(new prec[dim]);
        for (int icomp = 0; icomp < dim; icomp++)
        {
            temp_U[icomp] = curr_U[icomp];
        }
        U_magnitude[inod] = VECTOR::norm(temp_U, dim);
    }
}

void TOP_OPT::get_pressure_in_nodes_v(VECTOR &P, VECTOR &P_print)
{
    VECTOR_INT passed_nodes;
    passed_nodes.setZeros(physics.nNodes_v);
    P_print.reset(80);
    for (int inod = 0; inod < physics.nNodes; inod++)
    {
        P_print[inod] = P[inod];
        passed_nodes[inod] = 1;
    }
    int n_nodes_x_el = physics.dim + 1;
    int n_el_fact = 1; 
    for (int icount = 0; icount < physics.dim; icount++)
    {
        n_el_fact *=2;
    }
    for (int iel = 0; iel < physics.nElem; iel++)
    {
        VECTOR nodes_in_el(n_nodes_x_el);
        for (int inod = 0; inod < n_nodes_x_el; inod++)
        {
            nodes_in_el[inod] = physics.elem[iel][inod]; 
        }
        int start_el_v = n_el_fact * iel;
        for (int iloc = 0; iloc < n_nodes_x_el-1; iloc++)
        {
            int iglob = nodes_in_el[iloc];
            for (int jloc = iloc+1; jloc < n_nodes_x_el; jloc++)
            {
                int jglob = nodes_in_el[jloc];
                int glob_id_v = 0;
                switch (physics.dim)
                {
                    case 2:
                    {
                        if (iloc == 0)
                        {
                            glob_id_v = physics.elem_v[start_el_v+iloc][jloc];
                        }
                        else if (iloc == 1)
                        {
                            glob_id_v = physics.elem_v[start_el_v+iloc][jloc-1];
                        }
                        else
                        {
                            throw_line("ERROR: impossible case for 2D meshes\n");
                        }
                        break;
                    }
                    case 3:
                    {
                        glob_id_v = physics.elem_v[start_el_v+iloc][jloc];
                        break;
                    }
                }
                P_print[glob_id_v] = (P[iglob] + P[jglob]) / 2.0;
                passed_nodes[glob_id_v] = 1;
            }
        }
    }

    int not_pass = 1;
    for (int inod = 0; inod < physics.nNodes_v; inod++)
    {
        if (passed_nodes[inod] == 0)
        {
            std::cout << "inod: 0" << inod << "\n";
            not_pass = 0;
        }
    }
    if (not_pass == 0) pause();
}

void TOP_OPT::print_optimization_results(int &nNodes_v, int &dim, int &nNodes, int &loop, int &currLoopPrint, prec &obj, prec &change, VECTOR &gamma, bool &feasible, MATRIX &funcValues, MATRIX &no_weights_funcValues, VECTOR &changes, VECTOR_INT &valid, VECTOR &grad_gamma_norm)
{
    print_results_in_console(loop, obj, change);

    print_results_in_vtk(nNodes_v, dim, nNodes, loop, currLoopPrint, gamma, grad_gamma_norm);

    print_for_matlab_interface(loop, obj, change, feasible, funcValues, no_weights_funcValues, changes, valid);
}

void TOP_OPT::print_results_in_console(int &loop, prec &obj, prec &change)
{
    // prec VolPerc = Vol/V0*100;
    prec changePerc = change*100;

    VECTOR rel_func_val = Optimizer.func_val / obj * 100;

    std::cout << "\n-| It. " << loop << " <OPTIMIZATION RESULTS>| \n";
    std::cout << "--| VOLUME INFO| " << "\n";
    std::cout << "    |-----> Opt Box Vol: " << V0 << "\n"; 
    std::cout << "    |-----> Current Vol: " << Vol << "\n";
    
    std::cout << "--| FUNCTIONAL INFO|" << "\n";
    std::cout << "    |-----> Initial obj: " << Optimizer.f0Init << "\n";
    std::cout << "    |-----> Non-normalized obj: " << Optimizer.obj_abs_val << "\n";
    if (Optimizer.customFunc == 1)
    {
        std::cout << "    |-------> Inertial Part  | abs: " << Optimizer.func_val[0]  << "; rel: " << rel_func_val[0]  << "\n";
        std::cout << "    |-------> Shear Stress   | abs: " << Optimizer.func_val[1]  << "; rel: " << rel_func_val[1]  << "\n";
        std::cout << "    |-------> Vorticity      | abs: " << Optimizer.func_val[2]  << "; rel: " << rel_func_val[2]  << "\n";
        std::cout << "    |-------> Inlet Pressure | abs: " << Optimizer.func_val[3]  << "; rel: " << rel_func_val[3]  << "\n";
    }
    std::cout << "--| FUNCTIONAL   | obj:" << obj << "\n";
    // std::cout << "--| PERCENT VOL. | vol: " << VolPerc << "% \n";
    for (int iconstr = 0; iconstr < Optimizer.constraints.n_constr; iconstr++)
    {
        int type = Optimizer.constraints.list[iconstr].type;
        switch (type)
        {
            case 0:
            {
                prec vol_perc = Optimizer.constraints.list[iconstr].vol / Optimizer.constraints.list[iconstr].vol_0 * 100;
                std::cout << "--| PERCENT VOL. | vol: " << vol_perc << "% \n";
                break;
            }
            case 1:
            {
                prec vol_perc = Optimizer.constraints.list[iconstr].vol / Optimizer.constraints.list[iconstr].vol_0 * 100;
                std::cout << "--| PERCENT VOL. | vol:" << vol_perc << "% \tsubdomain: " << Optimizer.constraints.list[iconstr].domain_id << "\n";
                break;
            }
            case 2:
            {
                prec surf_perc = Optimizer.constraints.list[iconstr].surf / Optimizer.constraints.list[iconstr].surf_0 * 100;
                std::cout << "--| PERCENT SUR. | surf: " << surf_perc << "% \tbound: " << Optimizer.constraints.list[iconstr].bound_id << "\n";
                break;
            }
            case 3:
            {
                prec discratization_res = Optimizer.constraints.list[iconstr].discretization_res;
                std::cout << "--| DISCRET RES. | res: " << discratization_res << "\n";
                break;
            }
            case 4:
            {
                prec WSS = Optimizer.constraints.list[iconstr].actual_WSS;
                prec crit_WSS = Optimizer.constraints.list[iconstr].critical_WSS;
                int sign = Optimizer.constraints.list[iconstr].sign;
                std::string WSS_str = "(>";
                if (sign == -1)
                {
                    WSS_str = "(<";
                }
                WSS_str +=  std::to_string(crit_WSS) + ")";
                std::cout << "--| WSS " << WSS_str << " | WSS: " << WSS << "\n";
                break;
            }
            default:
            {
                std::cout << "\ntype: " << type << "\n";
                throw_line("ERROR: not handled constraint case\n");
                break;
            }
        }
    }
    std::cout << "--| GAMMA CHANGE | change: " << changePerc << " %\n";
    switch (physics.turn_on_gamma_acc)
    {
        case 0:
        {
            std::cout << "--| GAMMA ACC.   | state: OFF \t beta: " << Optimizer.beta_proj << "\n\n";
            break;
        }
        case 1:
        {
            std::cout << "--| GAMMA ACC.   | state: BASE\t beta: " << Optimizer.beta_proj << "\n\n";
            break;
        }
        case 2:
        {
            std::cout << "--| GAMMA ACC.   | state: ON\t beta: " << Optimizer.beta_proj << "\n\n";
            break;
        }
        case 3:
        {
            std::cout << "--| GAMMA ACC.   | state: FULL\t beta: " << Optimizer.beta_proj << "\n\n";
            break;
        }
        default:
        {
            throw_line("ERROR: invalid gamma acc case.\n");
        }
        break;
    }
    
}

void TOP_OPT::print_results_in_vtk(int &nNodes_v, int &dim, int &nNodes, int &loop, int &currLoopPrint, VECTOR &gamma_print, VECTOR &grad_gamma_norm)
{
    if (loop == currLoopPrint)
    {
        VECTOR alpha_print(gamma_print.length);
        for (int inod = 0; inod < alpha_print.length; inod++)
        {
            alpha_print[inod] = physics.alpha_min + (physics.real_alpha_max - physics.alpha_min) * physics.q * gamma_print[inod] / ( physics.q + 1 - gamma_print[inod]);
        }
        
        MATRIX U_print(nNodes_v, dim);
        VECTOR P(nNodes);
        NS.getUPFromSol_v(lastSolNS, U_print, P);
        VECTOR U_magnitude(nNodes_v);
        VECTOR P_print(nNodes_v);
        prepare_solution_print(U_print, P, U_magnitude, P_print);
        MATRIX U_print_T = U_print.transpose();

        // VECTOR WSS;
        // VECTOR normal(dim);
        // normal[0] = 0.0; normal[1] = -1.0;
        // physics.eval_WSS(U_print_T, 1, normal, WSS);
        // WSS.print();
        // std::vector<VECTOR> U_magnitude_gradient;
        // physics.eval_gradient(U_magnitude, U_magnitude_gradient);
        // std::vector<VECTOR> U_dir_gradient;
        // physics.eval_directional_gradient(U_print_T, nodes, directions, U_dir_gradient);
        // VECTOR U_x_dir_gradient = U_dir_gradient[0];
        
        switch (write_on_velocity_mesh)
        {
            case 0: // pressure nodes
            {
                //  VTKWriter.write(physics.coord_v, physics.elem_v, loop, gamma, 1, "Gamma", U_print, dim, "Velocity", P_print, 1, "Pressure", U_magnitude, 1, "||Velocity||");
                // VTKWriter.write(physics.coord, physics.elem, loop, gamma_print, 1, "Gamma", U_print, dim, "Velocity", P_print, 1, "Pressure", U_magnitude, 1, "||Velocity||", grad_gamma_norm, 1, "Grad(gamma)");
                VTKWriter.write(physics.coord, physics.elem, loop, gamma_print, 1, "Gamma", alpha_print, 1, "Inv. Perm.", U_print, dim, "Velocity", P_print, 1, "Pressure", U_magnitude, 1, "||Velocity||", grad_gamma_norm, 1, "Grad(gamma)");
                // VTKWriter.write(physics.coord_v, physics.elem_v, loop, gamma, 1, "Gamma", alpha, 1, "Alpha", U_print, dim, "Velocity", P_print, 1, "Pressure", U_magnitude, 1, "||Velocity||", grad_gamma_norm, 1, "Grad(gamma)", U_x_dir_gradient, 1, "d(u_x)/dx)");
                break;
            }
            case 1: // velocity nodes
            {
                //  VTKWriter.write(physics.coord_v, physics.elem_v, loop, gamma, 1, "Gamma", U_print, dim, "Velocity", P_print, 1, "Pressure", U_magnitude, 1, "||Velocity||");
                // VTKWriter.write(physics.coord_v, physics.elem_v, loop, gamma_print, 1, "Gamma", U_print, dim, "Velocity", P_print, 1, "Pressure", U_magnitude, 1, "||Velocity||", grad_gamma_norm, 1, "Grad(gamma)");
                VTKWriter.write(physics.coord_v, physics.elem_v, loop, gamma_print, 1, "Gamma", alpha_print, 1, "Inv. Perm.", U_print, dim, "Velocity", P_print, 1, "Pressure", U_magnitude, 1, "||Velocity||", grad_gamma_norm, 1, "Grad(gamma)");
                // VTKWriter.write(physics.coord_v, physics.elem_v, loop, gamma, 1, "Gamma", alpha, 1, "Alpha", U_print, dim, "Velocity", P_print, 1, "Pressure", U_magnitude, 1, "||Velocity||", grad_gamma_norm, 1, "Grad(gamma)", U_x_dir_gradient, 1, "d(u_x)/dx)");
                break;
            }
        }
        currLoopPrint += deltaPrint;
    }
}

void TOP_OPT::print_for_matlab_interface(int &loop, prec &obj, prec &change, bool &feasible, MATRIX &funcValues, MATRIX &no_weights_funcValues, VECTOR &changes, VECTOR_INT &valid)
{
    funcValues[loop-1][0] = obj;
    funcValues[loop-1][1] = Optimizer.func_out_box;
    for (int ifunc = 0; ifunc < Optimizer.fWeights.length; ifunc++)
    {
        funcValues[loop-1][ifunc+2] = Optimizer.func_val[ifunc];
        no_weights_funcValues[loop-1][ifunc] = Optimizer.no_weighted_func_val[ifunc];
    }        
    MATRIX temp_func_values(loop+1, Optimizer.fWeights.length+2);
    MATRIX temp_no_weights_func_values(loop, Optimizer.fWeights.length);
    temp_func_values[0][0] = Optimizer.f0Init;
    temp_func_values[0][1] = Optimizer.func_out_box_init;
    for (int ifunc = 0; ifunc < Optimizer.fWeights.length; ifunc++)
    {
        temp_func_values[0][ifunc+2] = Optimizer.func_val_init[ifunc];
    }
    for (int iloop = 1; iloop < loop+1; iloop++)
    {
        for (int ifunc = 0; ifunc < (Optimizer.fWeights.length+2); ifunc++)
        {
            temp_func_values[iloop][ifunc] = funcValues[iloop-1][ifunc];
            if (ifunc < Optimizer.fWeights.length)
            {
                temp_no_weights_func_values[iloop-1][ifunc] = no_weights_funcValues[iloop-1][ifunc];
            }
        }
    }

    evaluate_total_energy();
    physics.total_energy.append(temp_fluid_energy[0]);
    physics.fluid_energy.append_row(temp_fluid_energy);

    std::string path = "results/" + NS.name+ "/" + name + "/Matlab_interface";
    fs::create_directories(path);
    temp_func_values.print((path  + "/func.txt").c_str());
    temp_no_weights_func_values.print((path + "/no_weight_func.txt").c_str());
    physics.total_energy.print((path + "/total_energy.txt").c_str());
    physics.fluid_energy.print((path + "/fluid_energy.txt").c_str());
    changes.print((path + "/changes.txt").c_str());
    valid.print((path + "/valid.txt").c_str());
}

void TOP_OPT::evaluate_total_energy()
{
    temp_fluid_energy.resetZeros();

    int nNodes_v = physics.nNodes_v;

    VECTOR time_weights(physics.solution_times.length);
    time_weights.resetZeros();
    time_weights[0] = physics.solution_deltaT[0] / 2;
    for (int itime_step = 1; itime_step < (time_weights.length-1); itime_step++)
    {
        time_weights[itime_step] = (physics.solution_deltaT[itime_step-1] + physics.solution_deltaT[itime_step]) / 2;
    }
    time_weights[time_weights.length-1] = physics.solution_deltaT[physics.solution_deltaT.length-1] / 2;
    time_weights /= physics.t_end;

    for (int itime = 0; itime < time_weights.length; itime++)
    {
        VECTOR temp_velocity(nNodes_v);
        prec temp_velocity_energy = 0.0;
        prec temp_pressure_energy = 0.0;
        for (int icomp = 0; icomp < physics.dim; icomp++)
        {
            temp_velocity.resetZeros();
            for (int inod = 0; inod < nNodes_v; inod++)
            {
                temp_velocity[inod] = physics.NS_solution[itime][icomp*nNodes_v+inod];
            }
            temp_velocity_energy += temp_velocity.dot(NS.M * temp_velocity);
        }
        temp_fluid_energy[0] += temp_velocity_energy * time_weights[itime];

        for (int iel = 0; iel < physics.nElem; iel++)
        {
            prec temp_volume = physics.Volume[iel];
            for (int inod = 0; inod < physics.nNodes; inod++)
            {
                temp_pressure_energy += physics.NS_solution[itime][physics.dim*nNodes_v+inod] * temp_volume / (physics.dim + 1);
            }
        }
        temp_fluid_energy[1] += temp_pressure_energy * time_weights[itime];
    }
    temp_fluid_energy[0] *= (physics.rho / 2.0);
}

//----------------------------------------
// Define the elements and the node of the optimization domain region.
void TOP_OPT::handle_optimization_domain()
{
    std::cout << "\n----------\n--| SET OPTIMIZATION DOMAIN  |--\n----------\n";
    int nElem_v = physics.nElem_v;
    is_elem_in_dom.setZeros(nElem_v);
    physics.elems_in_doms.resize(n_subdomains);
    for (int idom = 0; idom < n_subdomains; idom++)
    {
        physics.elems_in_doms[idom].initialize(nElem_v);
    }
    VECTOR subd_el_count; subd_el_count.setZeros(n_subdomains);
    int nodes_per_el = physics.dim + 1;
    int nNodes_v = physics.nNodes_v;
    optNodeFromGlobNode.setZeros(nNodes_v);
    optNodeFromGlobNode += - 1;
    is_node_in_dom.setZeros(nNodes_v);
    is_node_in_dom += - 2; //set all to -2: -2=not passed, -1=not in opt domain, others=relative domain id
    not_opt_nodes_in_subdomains.setZeros(nNodes_v);
    not_opt_nodes_in_subdomains += -1;
    time_profiler.set();
    nodeInDom.setZeros(nNodes_v);
    elemInDom.setZeros(nElem_v);
    if (n_subdomains == 1) // optimizing on all the domain
    {
        for (int inod = 0; inod < nNodes_v; inod++)
        {
            nodeInDom[inod] = inod;
            optNodeFromGlobNode[inod] = inod;
        }
        for (int iel = 0; iel < nElem_v; iel++)
        {
            elemInDom[iel] = iel;
        }
        nNodeInDom = nNodes_v;
        nElemInDom = nElem_v;
    }
    else
    {
        //first check if is elem_in_dom
        for (int iel = 0; iel < nElem_v; iel++)
        {
            int temp_geo_id = physics.elem_geo_entities_ids_v[iel];
            physics.elems_in_doms[temp_geo_id][subd_el_count[temp_geo_id]] = iel;
            subd_el_count[temp_geo_id] += 1;
            if (opt_subdomains.hasIn(temp_geo_id))
            {
                is_elem_in_dom[iel] = 1;
            }
        }
        for (int idom = 0; idom < n_subdomains; idom++)
        {
            physics.elems_in_doms[idom].shrink(subd_el_count[idom]);
        }

        //check node in dom
        for (int iel = 0; iel < nElem_v; iel++)
        {
            int temp_elem_geo_id = physics.elem_geo_entities_ids_v[iel];
            int temp_elem_id_relevance;
            opt_subdomains.hasIn(temp_elem_geo_id, temp_elem_id_relevance);
            if (is_elem_in_dom[iel] == 0)
            {
                for (int iloc = 0; iloc < nodes_per_el; iloc++)
                {
                    int iglob = physics.elem_v[iel][iloc];
                    is_node_in_dom[iglob] = -1;
                    if (not_opt_nodes_in_subdomains[iglob] == -1) 
                    {
                        not_opt_nodes_in_subdomains[iglob] = temp_elem_geo_id; // the sumdomain priority of common nodes for the not_opt subdomains is given by their numbering
                    }
                    else if (temp_elem_geo_id < not_opt_nodes_in_subdomains[iglob])
                    {
                        not_opt_nodes_in_subdomains[iglob] = temp_elem_geo_id;
                    }
                }
            }
            else //  is_elem_in_dom[iel] == 1
            {
                for (int iloc = 0; iloc < nodes_per_el; iloc++)
                {
                    int iglob = physics.elem_v[iel][iloc];
                    int temp_node_geo_id = is_node_in_dom[iglob];
                    if (temp_node_geo_id == -2)
                    {
                        is_node_in_dom[iglob] = temp_elem_geo_id;
                    }
                    else if (temp_node_geo_id != -1) //not yet passed node
                    {
                        int temp_node_id_relevance;
                        opt_subdomains.hasIn(temp_node_geo_id, temp_node_id_relevance);
                        if (temp_elem_id_relevance < temp_node_id_relevance) // in this way yhe order of relevance of the subdomains belonging is the one inserted by the user
                        {
                            is_node_in_dom[iglob] = temp_elem_geo_id;
                        }
                    }
                }
            } 
        }
        int node_count = 0;

        for (int inod = 0; inod < nNodes_v; inod++)
        {
            switch (is_node_in_dom[inod])
            {
                case -2:
                {
                    throw_line("ERROR: not defined case of is_node_in_dom.\n");
                    break;
                }
                case -1:
                    break;
                default:
                {
                    optNodeFromGlobNode[inod] = node_count;
                    nodeInDom[node_count] = inod;
                    node_count++;
                    break;
                }
                
            }
        }
        nodeInDom.shrink(node_count);
        nNodeInDom = nodeInDom.length;

        // second check if is elem_in_dom: the idea is now to exclude elements in which one node belongs also to a non in dom element.
        int elem_count = 0;
        for (int iel = 0; iel < nElem_v; iel++)
        {
            int temp_is_in_dom = 1;
            for (int iloc = 0; iloc < nodes_per_el; iloc++)
            {
                int iglob = physics.elem_v[iel][iloc];
                if (is_node_in_dom[iglob] == -1) 
                {
                    temp_is_in_dom = 0;
                    break;
                }
            }
            is_elem_in_dom[iel] = temp_is_in_dom;
            if (temp_is_in_dom == 1) 
            {
                elemInDom[elem_count] = iel;
                elem_count++;
            }
        }
        elemInDom.shrink(elem_count);
        nElemInDom = elemInDom.length;
        // elemInDom.printRowMatlab("nodeInDom");
    }
    
    // finally get the optimization box to be able to use the gamma filtering
    optBox.initialize(physics.dim, 2);
    for (int idim = 0; idim < physics.dim; idim++)
    {
        optBox[idim][0] = 1e10;
        optBox[idim][1] = -1e10;
    }   
    for (int inod = 0; inod < nNodeInDom; inod ++)
    {
        int iglob = nodeInDom[inod];
        for (int icomp = 0; icomp < physics.dim; icomp++)
        {
            if (physics.coord_v[iglob][icomp] < optBox[icomp][0])
            {
                optBox[icomp][0] = physics.coord_v[iglob][icomp];
            }
            if(physics.coord_v[iglob][icomp] > optBox[icomp][1])
            {
                optBox[icomp][1] = physics.coord_v[iglob][icomp];
            }
        }
    }
}

void TOP_OPT::handle_gamma_initial_condition(VECTOR &gamma_opt, VECTOR &gamma)
{
    if (n_subdomains == 1)
    {
        prec init_value = subdomains_initial_values[0]; //only one subdomain
        gamma_opt.reset(init_value);
        gamma.reset(init_value);
    }
    else
    {
        for (int inod = 0; inod < physics.nNodes_v; inod++)
        {
            int temp_node_domain = is_node_in_dom[inod];
            if (temp_node_domain != -1)
            {
                prec temp_init_value = subdomains_initial_values[temp_node_domain];
                int temp_opt_node = optNodeFromGlobNode[inod];
                gamma_opt[temp_opt_node] = temp_init_value;
                gamma[inod] = temp_init_value;
            }
            else 
            {
                int temp_not_node_domain = not_opt_nodes_in_subdomains[inod];
                if (temp_not_node_domain == -1)
                {
                    throw_line("ERROR: not is_node_in_dom node is not classified as not_opt_nodes_in_subdomains\n");
                }
                else
                {
                    prec temp_init_value = subdomains_initial_values[temp_not_node_domain];
                    gamma[inod] = temp_init_value;
                }
            }
        }
    }
    
}

//-----------------------------------
//---| TOPOLOGY OPTIMIZATION SOLVER |
//-----------------------------------
void TOP_OPT::solve()
{
    prec startTime = omp_get_wtime();

    //------------------------------
    std::cout << "\n\n----------| NAVIER-STOKES TOPOLOGY OPTIMIZATION SOLVER |----------\n";
    std::cout << "\n----------| GEOMETRIC CONFIGURATION: ";
    std::cout << "\n-------------| dimension: " << physics.dim;
    std::cout << "\n-------------| n° nodes: " << physics.nNodes;
    std::cout << "\n-------------| n° velocity nodes: " << physics.nNodes_v;
    std::cout << "\n-------------| n° elements: " << physics.nElem;
    std::cout << "\n-------------| n° velocity elements: " << physics.nElem_v;
    std::cout << "\n-------------| n° d.o.f.: " << physics.nDof;
    std::cout << "\n----------| PARALLEL CONFIGURATION: " << PARALLEL::nThread << " threads\n\n";
    //------------------------------

    int nNodes_v = physics.nNodes_v;
    int nNodes = physics.nNodes;
    int dim = physics.dim;
    
    //---------------------
    // Preapre TopOpt
    int loop  = 0;
    prec change = change_toll + 1;
    prec obj = 1000;
    int currLoopPrint = 1;
    if (write_inital_condition == 1)
    {
        currLoopPrint = 0; // loop(=0) != currLoopPrint(=1)
        VTKWriter.opt_initial_cond_rescale = 1;
    }
    bool feasible = false;
    
    // Prepare Physics Solver
    prepareNS();
    // Prepare Adjoint Solver
    prepareADJ();

    // Prepare Gamma
    VECTOR gamma;
    gamma.setZeros(nNodes_v);
    gamma += 1;
    VECTOR gammaNew;
    VECTOR gammaOpt(nNodeInDom);
    handle_gamma_initial_condition(gammaOpt, gamma);
    gammaNew = gamma;
    MATRIX grad_gamma(nNodes_v, dim);
    VECTOR grad_gamma_norm;
    grad_gamma_norm.setZeros(nNodes_v);
    VECTOR grad_gamma_opt_norm(nNodeInDom);
    VECTOR gamma_print = (gamma-1.0)*(-1.0);

    // Print Initial Conditions
    print_results_in_vtk(nNodes_v, dim, nNodes, loop, currLoopPrint, gamma_print, grad_gamma_norm);

    MATRIX U(dim, nNodes_v);
    MATRIX Ua(dim, nNodes_v);

    MATRIX funcValues(maxIt, (Optimizer.fWeights.length + 2));
    MATRIX no_weights_funcValues(maxIt, (Optimizer.fWeights.length + 2));
    VECTOR changes;
    VECTOR_INT valid;

    temp_fluid_energy.setZeros(2);

    //UNTILL GAMMA CONVERGENCE DO:
    while (((change > change_toll || !(feasible)) && loop < maxIt) || loop < minIt)
    {
        loop++;
        physics.curr_opt_it = loop;

        //-------------------------------
        // UPDATE ALPHA GIVEN GAMMA
        //-------------------------------
        physics.update_alpha(gamma, alpha);

        // SOLVE NS
        NS.resetPrint(loop);

        prec startNStime = omp_get_wtime();
        solveNS();
        prec endNStime = omp_get_wtime();
        std::cout << "|PHYSICS| TOTAL SOLUTION TIME: " << endNStime-startNStime << "\n";

        // SOLVE ADJOINT
        ADJ.resetPrint(loop);
        prec startADJtime = omp_get_wtime();
        solveADJ();
        prec endADJtime = omp_get_wtime();
        std::cout << "|ADJOINT| TOTAL SOLUTION TIME: " << endADJtime-startADJtime << "\n";

        // if (NS.completeLog < 2) printf("\n-------\n-| OPTIMIZER UPDATE PARAMETERS |--\n-------\n");
        // Optimizer.updateSol(NS.globIter);

        printf("\n-------\n-| OPTIMIZER SOLVER |--\n-------\n");

        switch (optimization_scheme)
        {
            case 0: // MMA
            {
                Optimizer.solveMMA(gammaOpt, Vol, obj); // the filterd and projected value of gamma is saved in the optimized solution in the updateVal method
                break;
            }
            case 1: // GOC
            {
                Optimizer.solveGOC(gammaOpt, Vol, obj); // the filterd and projected value of gamma is saved in the optimized solution in the updateVal method
                break;
            }
            case 2: // GCMMA
            {
                Optimizer.solveGCMMA(gammaOpt, Vol, obj); // the filterd and projected value of gamma is saved in the optimized solution in the updateVal method
                break;
            }
        }
        
        save_gammaOpt_in_gammaFull(gammaOpt, gammaNew);

        change = (gamma-gammaNew).norm() / gamma.norm();
        physics.gamma_change = change;
        changes.append(change);
        gamma = gammaNew;
        gamma_print = (gamma-1.0)*(-1.0);
        prec vol_fract = Vol / V0;
        if (vol_fract <= Vr)
        {
            feasible = true;
            valid.append(1);
        }
        else
        {
            feasible = false;
            valid.append(0);
        }
        physics.Vol = Vol;
        physics.vol_fract = vol_fract;
        physics.gamma_max = gamma.max();

        // if (enableDiffusionFilter > 0)
        // {
        //     eval_gamma_gradient_norm_with_filter(gammaOpt, grad_gamma_opt_norm);
        // }
        // save_gammaOpt_in_gammaFull(grad_gamma_opt_norm, grad_gamma_norm);

        std::vector<VECTOR> gamma_gradient;
        physics.eval_gradient(gamma_print, gamma_gradient);
        VECTOR gamma_gradient_norm;
        physics.eval_gradient_norm(gamma_gradient, gamma_gradient_norm);
    
        print_optimization_results(nNodes_v, dim, nNodes, loop, currLoopPrint, obj, change, gamma_print, feasible, funcValues, no_weights_funcValues, changes, valid, gamma_gradient_norm);
        // pause();
    }

    // -----------------------------------
    // EXPORT OPTIMIZAED GEOMETRY NODES
    // -----------------------------------
    prec gammaMinOptMesh = 0.85; // minimum value of gamma to save node coords in the optimized mesh
    MATRIX_INT optElem(physics.nElem_v, dim + 1);
    exportOptimizedDomain(gamma, gammaMinOptMesh, optElem);
    VTKWriter.writeMesh(physics.nNodes_v, optElem.nRow, physics.coord_v, optElem);
    // VTKWriter.closeTFile();

    prec endTime = omp_get_wtime();

    solution_time = endTime - startTime;
}
//-------------------------



