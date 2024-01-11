#include "TopOpt.h"

TOP_OPT::TOP_OPT(std::string InputFile)
{
    inputFile = InputFile;
    printf("\n-------\n--| INITIALIZING NS PROBLEM |--\n-------\n");

    PHYSICS* tempP = &physics;
    NS.initialize(tempP, inputFile, alpha, false);

    printf("\n-------\n--| INITIALIZING ADJOINT PROBLEM |-- \n-------\n");
    ADJ.initialize(NS, alpha);
    //initialize TopOpt parameters
    printf("\n-------\n-| INITIALIZING TOP OPT PROBLEM |--\n-------\n");
    importParameters("INPUT_FILES/TopOptInput.txt");
    
    if (enableDiffusionFilter > 0) 
    {
        //throw_line("PAUSE0");
        Optimizer.diffusionFilter.initialize(enableDiffusionFilter, tempP, nNodeInDom, nodeInDom, optNodeFromGlobNode ,optBox, diffusionRadiusPercentage, Optimizer.diffusion_filter_case);
        //throw_line("PAUSE1");
        //diffusionFilter.printNeighbourhood();
    }
}

void TOP_OPT::importParameters(std::string inputFile)
{
    std::ifstream ParameterFile;
    ParameterFile.open(inputFile);

    MATRIX boxOpt(3,2);
    std::string line;
    std::istringstream iss;
    //
    STREAM::getLines(ParameterFile, line, 2);
    STREAM::getValue(ParameterFile, line, iss, name);
    //
    STREAM::getLines(ParameterFile, line, 2);
    bool isStat;
    STREAM::getValue(ParameterFile, line, iss, isStat);
    physics.isStationary = isStat;
    if (isStat == 1)
    {
        throw_line("ERROR: non updated solution case. Message time: 21/12/2023\n");
    }

    //
    STREAM::getLines(ParameterFile, line, 2);
    int flag_subdomains;
    STREAM::getValue(ParameterFile, line, iss, flag_subdomains);
    std::getline(ParameterFile, line);
    STREAM::getValue(ParameterFile, line, iss, n_subdomains);
    std::getline(ParameterFile, line);
    opt_subdomains.initialize(n_subdomains);
    STREAM::getRowVector(ParameterFile, line, iss, opt_subdomains);
    handle_optimization_domain();
    std::getline(ParameterFile, line);
    opt_subdomains_initial_values.initialize(n_subdomains);
    STREAM::getColVector(ParameterFile, line, iss, opt_subdomains_initial_values, n_subdomains);

    //
    STREAM::getLines(ParameterFile, line, 2);
    STREAM::getValue(ParameterFile, line, iss, q);
    //
    STREAM::getLines(ParameterFile, line, 2);

    STREAM::getValue(ParameterFile, line, iss, Vr);
    //
    STREAM::getLines(ParameterFile, line, 2);
    STREAM::getValue(ParameterFile, line, iss, alpha_min);
    //
    STREAM::getLines(ParameterFile, line, 2);
    STREAM::getValue(ParameterFile, line, iss, alpha_max);
    //
    STREAM::getLines(ParameterFile, line, 5);
    STREAM::getValue(ParameterFile, line, iss, onlyGrad);
    //
    STREAM::getLines(ParameterFile, line, 2);
    STREAM::getValue(ParameterFile, line, iss, minIt);
    //
    STREAM::getLines(ParameterFile, line, 2);
    STREAM::getValue(ParameterFile, line, iss, maxIt);
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
    // std::cout << "V0: " << V0 << "\n";
    // pause();

    STREAM::getLines(ParameterFile, line, 2);
    STREAM::getValue(ParameterFile, line, iss, customFunc);
    STREAM::getLines(ParameterFile, line, 3);
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
    if (customFunc == 1)
    {
        Optimizer.initialize(&physics, nodeInDom, elemInDom, optNodeFromGlobNode, q, alpha_min, alpha_max, V0, Vr, customFunc, onlyGrad, beta, opt_acceleration_case, beta_MAX, beta_min, beta_interpolation, change_max, change_min, crit_change, crit_beta, enableDiffusionFilter);
    }
    else
    {
        throw_line("\nERROR: customFunc != 1\n");
    }
    
    getline(ParameterFile, line);
    STREAM::getValue(ParameterFile, line, iss, diffusionRadiusPercentage);
    getline(ParameterFile, line);
    STREAM::getValue(ParameterFile, line, iss, diffusionFilterWeight);
    
    //------------------------------
    // READ PRINT INFO
    //------------------------------
    STREAM::getLines(ParameterFile, line, 3);
    STREAM::getValue(ParameterFile, line, iss, binPrint);
    //
    STREAM::getLines(ParameterFile, line, 6);
    STREAM::getValue(ParameterFile, line, iss, flagPrint);

    //
    STREAM::getLines(ParameterFile, line, 2);
    STREAM::getValue(ParameterFile, line, iss, deltaPrint);

    //
    STREAM::getLines(ParameterFile, line, 2);
    STREAM::getValue(ParameterFile, line, iss, write_on_velocity_mesh);

    // CLOSE STREAMING
    ParameterFile.close();

    std::string folderName = NS.name + "/" + name;
    if (write_on_velocity_mesh == 0)
    {
        VTKWriter.initializeForTopOpt(folderName, (physics).dim, (physics).nNodes, (physics).nElem, deltaPrint, binPrint);
    }
    else if (write_on_velocity_mesh == 1)
    {
        VTKWriter.initializeForTopOpt(folderName, (physics).dim, (physics).nNodes_v, (physics).nElem_v, deltaPrint, binPrint);
    }
    else
    {
        throw_line("ERROR: invalid mesh option for print\n");
    }
    

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
}
//---
void TOP_OPT::solveNS() // NS solver
{
    printf("\n-------\n-| NAVIER-STOKES SOLVER |--\n-------\n\n");
    if (physics.isStationary == 1) NS.StatSolver();
    else NS.Solver();
    lastSolNS = NS.lastSol;
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
}
//---
void TOP_OPT::solveADJ() // ADJ solver
{
    if (NS.completeLog < 2) printf("\n-------\n-| ADJOINT STAT FORCING |--\n-------\n");
    // MATRIX U = NS.getVelocityFromSol(lastSolNS);
    ADJ.globIter = NS.globIter;
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


void TOP_OPT::solve()
{
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
    VECTOR gamma = VECTOR::zeros(nNodes_v) + 1;
    VECTOR gammaNew;
    VECTOR gammaOpt(nNodeInDom);

    handle_gamma_initial_condition(gammaOpt);
    
    save_gammaOpt_in_gammaFull(gammaOpt, gamma);
    gammaNew = gamma;

    
    //---------------------
    int loop  = 0;
    //prec toll = 1e-5;
    prec change = change_toll + 1;
    prec obj = 1000;
    int currLoopPrint = -1;
    bool feasible = false;
    bool printNSSol = false;
    if (flagPrint > 0) currLoopPrint = 1;
    if (flagPrint == 2) printNSSol = true;
    

    prepareNS();
    prepareADJ();

    MATRIX U(dim, nNodes_v);
    MATRIX Ua(dim, nNodes_v);

    MATRIX funcValues(maxIt, (Optimizer.fWeights.length + 2));
    MATRIX no_weights_funcValues(maxIt, (Optimizer.fWeights.length + 2));
    VECTOR changes(maxIt);
    VECTOR_INT valid(maxIt);

    temp_fluid_energy.setZeros(2);

    //UNTILL GAMMA CONVERGENCE DO:
    while (((change > change_toll || !(feasible)) && loop < maxIt) || loop < minIt)
    {
        loop++;
        //-------------------------------
        // UPDATE ALPHA GIVEN GAMMA
        //-------------------------------
        prec factor = q*(alpha_max - alpha_min);
        for (int inod = 0; inod < nNodes_v; inod++)
        {
            alpha[inod] = alpha_min + factor * (1 - gamma[inod]) / ( q + gamma[inod]);
        }

        // SOLVE NS
        NS.resetPrint(loop);
        if (printNSSol) 
        {
            printNSSol = false;
        }
        solveNS();

        // SOLVE ADJOINT
        ADJ.resetPrint(loop);
        solveADJ();

        if (NS.completeLog < 2) printf("\n-------\n-| OPTIMIZER UPDATE PARAMETERS |--\n-------\n");
        Optimizer.updateSol(NS.globIter);

        printf("\n-------\n-| OPTIMIZER SOLVER |--\n-------\n");

        Optimizer.solveMMA(gammaOpt, Vol, obj);// the filterd and projected value of gamma is saved in the optimized solution in the updateVal method
        //gammaOpt = Optimizer.gamma_acc; 
       
        //Optimizer.solveGOC(gammaOpt, Vol, obj);
        //pause();
        save_gammaOpt_in_gammaFull(gammaOpt, gammaNew);

        change = (gamma-gammaNew).norm() / gamma.norm();
        physics.gamma_change = change;
        gamma = gammaNew;
        prec vol_fract = Vol / V0;
        if (vol_fract <= Vr)
        {
            feasible = true;
        }
        else
        {
            feasible = false;
        }
        physics.Vol = Vol;
        physics.vol_fract = vol_fract;
        physics.gamma_max = gamma.max();

        print_optimization_results(nNodes_v, dim, nNodes, loop, currLoopPrint, obj, change, printNSSol, gamma, feasible, funcValues, no_weights_funcValues, changes, valid);
        
    }

    if (flagPrint == 0)
    {
        VECTOR tempGamma = (gamma-1)*(-1);
        MATRIX U_print(nNodes_v, dim);
        VECTOR P_print(nNodes);
        NS.getUPFromSol(lastSolNS, U_print, P_print);
        VTKWriter.write(physics.coord_v, physics.elem_v, loop, tempGamma, 1, "Gamma", U_print, dim, "Velocity");//, P_print, 1, "Pressure");
        currLoopPrint += deltaPrint;
    }
    

    //-----------------------------------
    // EXPORT OPTIMIZAED GEOMETRY NODES
    //-----------------------------------
    prec gammaMinOptMesh = 0.85; // minimum value of gamma to save node coords in the optimized mesh
    MATRIX_INT optElem(physics.nElem_v, dim + 1);
    exportOptimizedDomain(gamma, gammaMinOptMesh, optElem);
    //std::cout << "nRows: " << optElem.nRow << "\n";
    VTKWriter.writeMesh(physics.nNodes_v, optElem.nRow, physics.coord_v, optElem);
    
    VTKWriter.closeTFile();
}
//-------------------------

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
                // std::cout << "glob_id: " << glob_id_v << "\n";
                // pause();
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

void TOP_OPT::print_optimization_results(int &nNodes_v, int &dim, int &nNodes, int &loop, int &currLoopPrint, prec &obj, prec &change, bool &printNSSol, VECTOR &gamma, bool &feasible, MATRIX &funcValues, MATRIX &no_weights_funcValues, VECTOR &changes, VECTOR_INT &valid)
{
    print_results_in_console(loop, obj, change);

    print_results_in_vtk(nNodes_v, dim, nNodes, loop, currLoopPrint, printNSSol, gamma);
    
    print_for_matlab_interface(loop, obj, change, feasible, funcValues, no_weights_funcValues, changes, valid);
}

void TOP_OPT::print_results_in_console(int &loop, prec &obj, prec &change)
{
    prec VolPerc = Vol/V0*100;
    prec changePerc = change*100;

    VECTOR rel_func_val = Optimizer.func_val / obj * 100;

    std::cout << "\n-| It. " << loop << " <OPTIMIZATION RESULTS>| \n";
    std::cout << "--| VOLUME INFO| " << "\n";
    std::cout << "  |-----> Opt Box Vol: " << V0 << "\n"; 
    std::cout << "  |-----> Current Vol: " << Vol << "\n";
    
    std::cout << "--| FUNCTIONAL INFO|" << "\n";
    std::cout << "  |-----> Initial obj: " << Optimizer.f0Init << "\n";
    std::cout << "  |-----> Non-normalized obj: " << Optimizer.obj_abs_val << "\n";
    if (Optimizer.customFunc == 1)
    {
        std::cout << "  |-------> Inertial Part  | abs: " << Optimizer.func_val[0]  << "; rel: " << rel_func_val[0]  << "\n";
        std::cout << "  |-------> Shear Stress   | abs: " << Optimizer.func_val[1]  << "; rel: " << rel_func_val[1]  << "\n";
        std::cout << "  |-------> Vorticity      | abs: " << Optimizer.func_val[2]  << "; rel: " << rel_func_val[2]  << "\n";
        std::cout << "  |-------> Inlet Pressure | abs: " << Optimizer.func_val[3]  << "; rel: " << rel_func_val[3]  << "\n";
    }
    std::cout << "--| FUNCTIONAL   | obj:" << obj << "\n";
    std::cout << "--| PERCENT VOL. | vol: " << VolPerc << "% \n";
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

void TOP_OPT::print_results_in_vtk(int &nNodes_v, int &dim, int &nNodes, int &loop, int &currLoopPrint, bool &printNSSol, VECTOR &gamma)
{
    if (loop == currLoopPrint)
        {
            if (flagPrint == 2) printNSSol = true;
            VECTOR tempGamma = (gamma-1)*(-1);
            VECTOR alpha(tempGamma.length);
            for (int inod = 0; inod < tempGamma.length; inod++)
            {
                alpha[inod] = alpha_min + (alpha_max - alpha_min) * q * tempGamma[inod] / ( q + 1 - tempGamma[inod]);
            }
            
            MATRIX U_print(nNodes_v, dim);
            VECTOR P(nNodes);
            NS.getUPFromSol_v(lastSolNS, U_print, P);
            VECTOR U_magnitude(nNodes_v);
            VECTOR P_print(nNodes_v);
            prepare_solution_print(U_print, P, U_magnitude, P_print);
            
            VTKWriter.write(physics.coord_v, physics.elem_v, loop, tempGamma, 1, "Gamma", alpha, 1, "Alpha", U_print, dim, "Velocity", P_print, 1, "Pressure", U_magnitude, 1, "||Velocity||");
            
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

    changes[loop-1] = change;
    
    if (Vol < Vr*V0) feasible = true;
    else feasible = false;
    valid[loop-1] = feasible;
    
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

    VECTOR temp_gamma_changes(loop);
    VECTOR temp_feasible_volume(loop);
    for (int iloop = 1; iloop < loop+1; iloop++)
    {
        temp_gamma_changes[iloop-1] = changes[iloop-1];
        temp_feasible_volume[iloop-1] = valid[iloop-1];
    }

    evaluate_total_energy();
    physics.total_energy.append(temp_fluid_energy[0]);
    physics.fluid_energy.append_row(temp_fluid_energy);

    std::string path = "MATLAB_INTERFACE/" + NS.name+ "/" + name;
    fs::create_directories(path);
    temp_func_values.print((path  + "/func.txt").c_str());
    temp_no_weights_func_values.print((path + "/no_weight_func.txt").c_str());
    physics.total_energy.print((path + "/total_energy.txt").c_str());
    physics.fluid_energy.print((path + "/fluid_energy.txt").c_str());
    temp_gamma_changes.print((path + "/changes.txt").c_str());
    temp_feasible_volume.print((path + "/valid.txt").c_str());
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
    //first check if is elem_in_dom
    int nElem_v = physics.nElem_v;
    is_elem_in_dom.setZeros(nElem_v);
    for (int iel = 0; iel < nElem_v; iel++)
    {
        int temp_geo_id = physics.elem_geo_entities_ids_v[iel];
        if (opt_subdomains.hasIn(temp_geo_id))
        {
            is_elem_in_dom[iel] = 1;
        }
    }
    // is_elem_in_dom.printRowMatlab("is_el");

    //check node in dom
    int nodes_per_el = physics.dim + 1;
    int nNodes_v = physics.nNodes_v;
    optNodeFromGlobNode.setZeros(nNodes_v);
    optNodeFromGlobNode += - 1;
    is_node_in_dom.setZeros(nNodes_v);
    is_node_in_dom += - 2; //set all to -2: -2=not passed, -1=not in opt domain, others=relative domain id
    for (int iel = 0; iel < nElem_v; iel++)
    {
        int temp_elem_geo_id = physics.elem_geo_entities_ids_v[iel];
        if (is_elem_in_dom[iel] == 0)
        {
            for (int iloc = 0; iloc < nodes_per_el; iloc++)
            {
                int iglob = physics.elem_v[iel][iloc];
                is_node_in_dom[iglob] = -1;
            }
        }
        else //  is_elem_in_dom[iel] == 1
        {
            for (int iloc = 0; iloc < nodes_per_el; iloc++)
            {
                int iglob = physics.elem_v[iel][iloc];
                if (is_node_in_dom[iglob] == -2) //not yet passed node
                {
                    is_node_in_dom[iglob] = temp_elem_geo_id;
                }
            }
        } 
    }
    // is_node_in_dom.printRowMatlab("is_nod");
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
                optNodeFromGlobNode[inod] = nodeInDom.length;
                nodeInDom.append(inod);
                break;
            }
            
        }
    }
    nNodeInDom = nodeInDom.length;
    // nodeInDom.printRowMatlab("nodeInDom");
    // optNodeFromGlobNode.printRowMatlab("optFrom");

    // second check if is elem_in_dom: the idea is now to exclude elements in which one node belongs also to a non in dom element.
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
            elemInDom.append(iel);
        }
    }
    nElemInDom = elemInDom.length;
    // elemInDom.printRowMatlab("nodeInDom");

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
    // MATRIX::printForMatlab(optBox, "optBox");
    // pause();
}

void TOP_OPT::handle_gamma_initial_condition(VECTOR &gamma_opt)
{
    for (int inod = 0; inod < nNodeInDom; inod++)
    {
        int iglob = nodeInDom[inod];
        int temp_node_domain = is_node_in_dom[iglob];
        int domain_id;
        if (!opt_subdomains.hasIn(temp_node_domain, domain_id)) throw_line("ERROR: node with domain ID not in the possible set");
        gamma_opt[inod] = opt_subdomains_initial_values[domain_id];
    }
}


