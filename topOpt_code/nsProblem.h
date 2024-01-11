#pragma once

#include "CODE_HEADERS/codeHeader.h"
#include "geometry.h"
//#include <sys/time.h>

class PROBLEM_NS
{
public:
    //---------
    std::string name;
    int nDof; 
    PHYSICS* physics;
    int completeLog;
    bool printRes = false;
    int save_solution;
    prec deltaT0;

    VECTOR alpha;

    //--- ABOUT TIME ---
    prec time;
    prec deltaT;
    int globIter = 0;
    VECTOR real_solution_times;

    //---  PROBLEM_NS ---
    int flagForcing;
    std::vector<std::string> statForcing;
    std::vector<std::string> timeForcing;

    VECTOR currForcingAtNodes;
    std::string statG;
    std::string timeG;
    int flagG;
    
    //--- SYSMAT ---------
    CSRMAT SYSMAT_base; CSRMAT SYSMAT; 
    int rescale;
    int nTerms; //terms of H
    int nEvals; //evaluation (before shrinking) of matrix with same H pattern
    int dimTerms; // number of terms for each components in the SYSMAT matrix = H.nTerm + B[i].nTerm
    VECTOR_INT iSparse; VECTOR_INT jSparse; VECTOR_INT realPos; 
    std::vector<std::vector<VECTOR_INT>> realPosSysmat; // std vector containing in each comp the pattern of H coeff in sysmatNS
    // CSRMAT H; CSRMAT M; std::vector<CSRMAT> B;
    CSRMAT H; 
    CSRMAT M;
    VECTOR N; 
    VECTOR Ma;
    std::vector<CSRMAT> B;
    int nTermB;
    // std::vector<std::vector<int>> AcPath; std::vector<std::vector<int>> AgPath;
    VECTOR rhs_statForcing;
    VECTOR rhs_timeForcing;
    VECTOR rhs;



    //--------------------------------------------------------------
    //BOUNDARY CONDITIONS
    //--------------------------------------------------------------
    int flagWall;
    bool statBCAlreadyApplied = false;
    //---- BOUND INFORMATIONS ---- ALL WITH REFERENCE TO THE VELOCITY MESH
    int nBoundNodes; //# of total boundary nodes
    int nBoundElems; //# of total boundary elements
    int nBound;      //# of boundary IDs 
    std::shared_ptr<int[]> nBoundIdNodes; //# of boundary nodes of a specif  ID

    std::shared_ptr<int[]> neuBoundIdElems_buff_buff;
    std::shared_ptr<int*[]> neuBoundIdElems_buff;
    std::shared_ptr<int**[]> neuBoundIdElems; // [entryID, specific el with such ID, nodes] 
    std::shared_ptr<int[]> nNeuBoundIdElems; //# of boundary elems of a specif  ID
    std::shared_ptr<int[]> boundIdNodes_buff;
    std::shared_ptr<int*[]> boundIdNodes; // [entryID, specific el with such ID] 
    //????double* areaboundEl; // area of each boundary element

    std::shared_ptr<int[]> neuTimeBoundIdElems_buff_buff;
    std::shared_ptr<int*[]> neuTimeBoundIdElems_buff;
    std::shared_ptr<int**[]> neuTimeBoundIdElems;
    std::shared_ptr<int[]> nNeuTimeBoundIdElems;

    std::shared_ptr<int[]> symmBoundIdElems_buff_buff;
    std::shared_ptr<int*[]> symmBoundIdElems_buff;
    std::shared_ptr<int**[]> symmBoundIdElems;
    std::shared_ptr<int[]> nSymmBoundIdElems;
    std::shared_ptr<int[]> symmPressureBoundIdElems_buff_buff;
    std::shared_ptr<int*[]> symmPressureBoundIdElems_buff;
    std::shared_ptr<int**[]> symmPressureBoundIdElems;

    std::shared_ptr<int[]> dirPressureBoundIdElems_buff_buff;
    std::shared_ptr<int*[]> dirPressureBoundIdElems_buff;
    std::shared_ptr<int**[]> dirPressureBoundIdElems;
    VECTOR_INT dirBoundPressureNEl;

    //- Inner Walls ---
    int nInnerBound;
    VECTOR_INT innerBound;

    //- Symmetries --- (priority -1) SYMMETRIES
    int nSymmBound;
    int    nSymm;
    VECTOR_INT   symmNod;
    VECTOR_INT symmRealIglob;
    MATRIX symmVal;

    VECTOR_INT symmPosToZero;
    VECTOR_INT symmPosToCoef;
    VECTOR_INT symmBound;
    CSRMAT J;
    
    //- Wall BC --- (priority: 5)
    int        nWallBound;
    VECTOR_INT wallBound;
    VECTOR_INT wallNod; 
    VECTOR_INT wallIdCount;
    int  nWall;

    VECTOR_INT wallPosToZero; VECTOR_INT wallPosToOne;

    //- Temporal Dirichlet BC --- (priority: 4)
    int        nDirTimeBound;
    VECTOR_INT dirTimeBound;
    std::vector<PIECEWISE_FUNCTION> dirTimeFunc;
    VECTOR_INT nIdDirTimeCases;

    VECTOR_INT dirTimeNod; 
    VECTOR_INT dirTimeIdCount;
    MATRIX dirTimeVal;
    int  nTimeDir;

    VECTOR_INT dirTimePosToZero; VECTOR_INT dirTimeIglobToZero; VECTOR_INT dirTimeCountToZero; VECTOR_INT dirTimePosToOne;

    //- Static Dirichlet BC --- (priority: 3)
    int        nDirBound;
    VECTOR_INT dirBound;
    std::vector<std::vector<std::string>> dirFunc;
    
    VECTOR_INT dirNod; 
    VECTOR_INT dirIdCount;
    MATRIX dirVal;
    int  nDir;
    //MATRIX dirCoord;

    VECTOR_INT dirPosToZero; VECTOR_INT dirIglobToZero; VECTOR_INT dirCountToZero; VECTOR_INT dirPosToOne;

    //- Temporal Neumann BC --- (priority: 2)
    int        nNeuTimeBound;
    VECTOR_INT neuTimeBound;

    std::vector<PIECEWISE_FUNCTION> neuTimeMeanP;
    VECTOR_INT nIdNeuTimeCases;

    VECTOR_INT neuTimeNod; 
    VECTOR_INT neuTimeIdCount;
    MATRIX neuTimeVal;
    int  nTimeNeu;
    VECTOR neuTimeBaseFlux; // baseFlux = areaEl/dim
    MATRIX neuTimeNormal;
    
    //- Static Neumann BC --- (priority: 1)
    int        nNeuBound;
    VECTOR_INT neuBound;
    VECTOR neuMeanP;

    VECTOR_INT   neuNod;
    VECTOR_INT neuIdCount;
    MATRIX neuVal;
    int    nNeu;

    //- Normal Velocity Open Boundary --- (priority: -2)
    int nNVOBBound;
    VECTOR_INT NVOBBound;
    
    VECTOR_INT NVOBNod; 
    VECTOR_INT NVOBIdCount;
    int  nNVOB;

    VECTOR_INT NVOBPosToZero; VECTOR_INT NVOBIglobToZero; VECTOR_INT NVOBCountToZero; VECTOR_INT NVOBPosToOne;

    // previous solution
    VECTOR lastSol;

    MATRIX simulation_times_solution;


    //------------------
    // PRINT parameters
    //------------------
    VTK VTKWriter;

    //----------------------------------------------------------
    // CONSTRUCTOR
    //----------------------------------------------------------
    PROBLEM_NS() {}
    PROBLEM_NS(PHYSICS* &Physics, std::string probRefFile)
    {
        printRes = true;
        physics = Physics;
        importParameters(probRefFile);
        checkImportParameters();
        importPREPRO();
        if (abs(time) < 1e-16) time = 0;
        localBasis();
    }
    //---
    PROBLEM_NS(int dimension)
    {
        (*physics).dim = dimension;
    }
    //---
    void initialize(PHYSICS *&Physics, std::string probRefFile, VECTOR &alphaIn, bool print = true)
    {
        printRes = print;
        physics = Physics;
        importParameters(probRefFile);
        checkImportParameters();
        importPREPRO();
        alphaIn.initialize((*physics).nNodes_v);
        alpha.length = alphaIn.length;
        alpha.P = alphaIn.P;
        if (abs(time) < 1e-16) time = 0;
        localBasis();
    }

    void printFriendlyError()
    {
        // std::cout << "ERROR: PCU NOT WORKING CORRECTLY...\n RESTARTING THE SYSTEM...\n";for (int j = 0; j < 2e9; j++) j+16;
        // std::cout << "segmentation fault\n"; for (int j = 0; j < 1e9; j++) j+16;
        // std::cout << "SYSTEM FAILED\n ERASING CRITICAL FILES.... WAIT FOR THE END OF THE PROCEDURE\n"; for (int j = 0; j < 2e9; j++) j+16;
        // printRandomLines(45);
        // throw_line("cancella pure quì :)\n");
    }
    void printRandomLines(int nRep)
    {
        // for (int i = 0; i < nRep; i++) 
        // {
        //     for (int j = 0; j < 1e8; j++) j+16;
        //     std::cout << "....................................................................................................................\n";
        // }
        // std::cout << "BUONA FORTUNA\n";
        // std::cout << "          -----            -----        \n";
        // std::cout << "          -----            -----        \n";
        // std::cout << "     |                                  | \n";
        // std::cout << "     |                                  | \n";
        // std::cout << "                     |                  \n";
        // std::cout << "     __                                __\n";
        // std::cout << "       __                           __\n";
        // std::cout << "         ___                    ___\n";
        // std::cout << "             ___            ___\n";
        // std::cout << "                 __________\n";
    }
    //-----------------------------------------
    // IMPORT PRE-PRO
    //-----------------------------------------
    void importPREPRO();

    //-----------------------------------------
    // INITIALIZE PROBLEM_NS PARAMETERS
    //-----------------------------------------
    void importParameters(std::string readFile);

    //------------------------------------------
    // INIT CONDITIONS
    //------------------------------------------
    void setInitCond();
    //------------------------------------------
    // SET BC
    //------------------------------------------    
    void setBC();
    //---------------------------------------------------
    // PREPARE LINEAR SYSTEM (ASSEMBLY OF THE MATRICES)
    //---------------------------------------------------
    void prepareSolver();

    //-------------------------------------------
    void localBasis();
    //-------------------------------------------
    // ASSEMBLE NS
    //-------------------------------------------
    void assemble();
    void assembleMa();
    void resizeCoef(std::shared_ptr<prec[]> &coef, VECTOR &resized);
    void addToSysmat(std::shared_ptr<prec[]> &coef, prec factor = 1.0);
    void addToSysmat(std::shared_ptr<prec[]> &coef, CSRMAT &mat, prec factor = 1.0);
    
    //-------------------------------------------
    // FORCING
    //-------------------------------------------
    void applyStatForcing();
    //-------------------------------------------
    // IMPOSE BOUNDARY CONDITIONS
    //-------------------------------------------
    void imposeBC(CSRMAT &SYSMAT, VECTOR &rhs);
    //---------------------------------------------------
    // SOLVER
    //---------------------------------------------------
    void Solver();

    //-------------------------
    void StatSolver()
    {
        time = 0;
        // prec t_end = 1e16;//(*physics).t_end;

        deltaT = 1e16;
        // while (time < t_end)
        // {
            globIter = -1;
            oneStepSolver(1e-2, 200);
            // if (time+deltaT > t_end) deltaT = t_end-time;
        // }
        VTKWriter.closeTFile();
    }
    //----
    void StatSolverIterative()
    {
        time = 0;
        prec t_end = (*physics).t_end;

        while (time < t_end)
        {
            globIter = -1;
            oneStepSolver(1e-2, 20);
            if (time+deltaT > t_end) deltaT = t_end-time;
        }
        VTKWriter.closeTFile();
    }
    //----
    void oneStepSolver(prec toll = 1e-2, int itMax = 20);

    //----
    void evaluate_solution_on_requested_time_steps();
    //--------------------------------------
    // ERROR EVALUATION
    //--------------------------------------
    //---------------------------
    // ERROR EVALUATION
    //---------------------------
    prec evalErrRelV(VECTOR &solh, VECTOR &solToCompare);

    prec evalErrRelV(std::shared_ptr<prec[]> &solh, std::shared_ptr<prec[]> &solToCompare);

    prec evalErrRelV(prec* &solh, prec* &solToCompare);
    
    prec evalErr(VECTOR &solh, VECTOR &solToCompare, prec &solNorm)
    {
        if (solToCompare.length != (*physics).nNodes) throw_line("ERROR: Comparing solution with another one defined in a different number of nodes\n");
        solNorm = norm(solToCompare);
        return evalErr(solh.P, solToCompare.P);
    }
    //----------------------------
    prec evalErr(VECTOR &solh, VECTOR &solToCompare)
    {
        if (solh.length != solToCompare.length) throw_line("ERROR: Comparing solution with another one defined in a different number of nodes\n");
        if (solh.length == (*physics).nNodes) return evalErr(solh.P, solToCompare.P);
        else if (solh.length == (*physics).nNodes_v) return evalErrV(solh.P, solToCompare.P);
        else throw_line("ERROR: Comparing solution with another one defined in a different number of nodes\n"); 
    }
    //----------------------------
    prec evalErr(std::shared_ptr<prec[]> &solh, std::shared_ptr<prec[]> &solToCompare)
    {
        prec* solhP = &(solh[0]);
        prec* solToCompareP = &(solToCompare[0]);
        return evalErr(solhP, solToCompareP);
    }
    //---
    prec evalErr(prec* &solh, prec* &solToCompare)
    {
        prec err = 0;
        std::shared_ptr<prec[]> Volume = (*physics).Volume.P;
        int dim = (*physics).dim;
        int nElem = (*physics).nElem;
        std::shared_ptr<int*[]> elem = (*physics).elem.PP;
        switch (dim)
        {
            case 2:
            {
                for (int iel = 0; iel < nElem; iel++)
                {
                    prec tempNodalVol = Volume[iel]/3; // 3 = nNodesxEl
                    for (int iloc = 0; iloc < dim + 1; iloc++)
                    {
                        int iglob = elem[iel][iloc];
                        prec nodalDiff = solh[iglob] - solToCompare[iglob];
                        nodalDiff *= nodalDiff;
                        err += nodalDiff * tempNodalVol; 
                    }

                }
                err = sqrt(err);
                break;
            }
            //---
            case 3:
            {
                for (int iel = 0; iel < nElem; iel++)
                {
                    prec tempNodalVol = Volume[iel]/4; // 4 = nNodesxEl
                    for (int iloc = 0; iloc < dim + 1; iloc++)
                    {
                        int iglob = elem[iel][iloc];
                        prec nodalDiff = solh[iglob] - solToCompare[iglob];
                        nodalDiff *= nodalDiff;
                        err += nodalDiff * tempNodalVol; 
                    }

                }
                err = sqrt(err);
                break;
            }
        }
        return err;
    }
    //----------------------------
    prec evalErrV(std::shared_ptr<prec[]> &solh, std::shared_ptr<prec[]> &solToCompare)
    {
        prec err = 0;
        std::shared_ptr<prec[]>  Volume = (*physics).Volume_v.P;
        int dim = (*physics).dim;
        int nElem = (*physics).nElem_v;
        std::shared_ptr<int*[]>  elem = (*physics).elem_v.PP;
        switch (dim)
        {
            case 2:
            {
                for (int iel = 0; iel < nElem; iel++)
                {
                    prec tempNodalVol = Volume[iel]/3; // 3 = nNodesxEl
                    for (int iloc = 0; iloc < dim + 1; iloc++)
                    {
                        int iglob = elem[iel][iloc];
                        prec nodalDiff = solh[iglob] - solToCompare[iglob];
                        nodalDiff *= nodalDiff;
                        err += nodalDiff * tempNodalVol; 
                    }

                }
                err = sqrt(err);
                break;
            }
            //---
            case 3:
            {
                for (int iel = 0; iel < nElem; iel++)
                {
                    prec tempNodalVol = Volume[iel]/4; // 4 = nNodesxEl
                    for (int iloc = 0; iloc < dim + 1; iloc++)
                    {
                        int iglob = elem[iel][iloc];
                        prec nodalDiff = solh[iglob] - solToCompare[iglob];
                        nodalDiff *= nodalDiff;
                        err += nodalDiff * tempNodalVol; 
                    }

                }
                err = sqrt(err);
                break;
            }
        }
        return err;
    }
    //----------------------------
    prec norm (VECTOR &solh);
    prec norm (std::shared_ptr<prec[]> &solh);
    prec norm (prec* &solh);
    //----------------------------

    void testNormal()
    {
        int i0; int i1; int i2;
        VECTOR v10(3); VECTOR v21(3); VECTOR normal(3); VECTOR P(3); VECTOR extNormal(3);
        double res;
        std::shared_ptr<prec*[]>  coord = (*physics).coord.PP;
        for (int i = 0; i < nBound; i++)
        {
            printf("BOUND id: %d\n-----------------\n", i);
            for(int iel = 0; iel < nNeuBoundIdElems[i]; iel++)
            {
                i0 = neuBoundIdElems[i][iel][0];
                i1 = neuBoundIdElems[i][iel][1];
                i2 = neuBoundIdElems[i][iel][2];
                for (int c = 0; c < 3; c++)
                {
                    P[c]   = (coord[i0][c] + coord[i1][c]+ coord[i2][c])/3;
                    v10[c] = coord[i1][c] - coord[i0][c];
                    v21[c] = coord[i2][c] - coord[i1][c];
                }
                VECTOR::cross(v10, v21, normal);
                 
                res = extNormal.dot(normal.P);
                // if (res == 0)
                // {
                    //VECTOR::print2(v10.P, v21.P, 3, "v10", "v21"); 
                    //printf("\nboundId: %d, nBoundEl: %d\n", i, iel);
                    //VECTOR::print2(normal.P, ext, 3, "normal", "extNormal");         
                    if (res < 0) printf("\nres: %f, INTERNAL :( \n", res);
                    else printf("\nres: %f, EXTERNAL :) \n", res);                    
                // }           
            }            
        }
    }

    //------------------------------------------
    // CLEAR 
    //------------------------------------------
    // void dlt()
    // {

    // }

    //------------------------------
    // RESET PRINT
    //------------------------------
    void resetPrint(int iter);

    //-------------------------------
    // PRINT ONE STEP SOLUTION IN VTK
    //-------------------------------
    void print_one_step_sol_in_VTK(int dim, int nNodes, int nNodes_v, prec time);
    void print_sol_in_VTK(MATRIX &requested_sol);
    //------------------------------
    // CHECK CORRECT PROBLEM SETTING
    //------------------------------
    void checkImportParameters()
    {
        std::cout << "\n+++++++++++++++| NS CHECK IMPORT-PARAMETERS |+++++++++++++++\n\n";

        // PROBLEM NAME & DIMENSION
        std::cout << " Problem Name: " << name << "\n" ;
        std::cout << " Dimension: " << (*physics).dim << "\n";

        // TIME PARAMETERS
        std::cout << "\n.-=| TIME PARAMETERS |=-. \n";
        std::cout << " Ending time:  " << (*physics).t_end << "\n" ;
        std::cout << " Default time step: " << deltaT << "\n" ;
        std::cout << " Minimal time step: " << (*physics).deltaT_min << "\n" ;

        // PHYSICS PARAMETERS
        std::cout << "\n.-=| PHYSICS PARAMETERS |=-. \n";
        std::cout << " rho: "<< (*physics).rho << "\n";
        std::cout << " mu: " << (*physics).mu << "\n" ;
        
        // FORCING & COMPRESSIBILITY
        std::cout << "\n.-=| FORCING & COMPRESSIBILITY |=-. \n";
        // forcing
        CUSTOM::printRowStd(statForcing, " STATIC Forcing: ");
        if (flagForcing == 1) CUSTOM::printRowStd(timeForcing, " TIME DEPENDENT Forcing: ");
        // compressibility G
        std::cout << " STATIC Compressibility G: " << statG << "\n";
        if (flagG == 1) std::cout << " TIME DEPENDENT Compressibility G: " << timeG << "\n";

        // BOUNDARY CONDITIONS
        std::cout << "\n.-=| BOUNDARY CONDITIONS |=-. \n";
        std::cout << " Flag BC: " << (*physics).flagBC << "\n";
        std::cout << " Flag Wall: " << flagWall << "\n";
        // Wall
        std::cout << "\n--| WALL:  ";
        std::cout << "\n # Wall Bound ID: " << nWallBound << "\n";
        wallBound.printRowMatlab(" Wall Bound ID");
        // Inner 
        std::cout << "\n--| INNER:  ";
        std::cout << "\n # Inner Bound ID: " << nInnerBound << "\n";
        wallBound.printRowMatlab(" Inner Bound ID");
        // Symmetry
        std::cout << "\n--| SYMMETRY:  ";
        std::cout << "\n # Symmetry Bound ID: " << nSymmBound << "\n";
        symmBound.printRowMatlab(" Symmetry Bound ID");
        // Static Dirichlet
        std::cout << "\n--| STATIC DIRICHLET:  ";
        std::cout << "\n # STATIC Dirichlet Bound ID: " << nDirBound << "\n";
        dirBound.printRowMatlab(" STATIC Dirichlet Bound ID");
        for (int ibound = 0; ibound < nDirBound; ibound++)
        {
            std::string tempName = " " + std::to_string(ibound+1) + "th STATIC Dirichlet function: ";
            CUSTOM::printRowStd(dirFunc[ibound], tempName);
        }
        // Static Neumann
        std::cout << "\n--| STATIC NEUMANN:  ";
        std::cout << "\n # STATIC Neumann Bound ID: " << nNeuBound << "\n";
        neuBound.printRowMatlab(" STATIC Neumann Bound ID");
        neuMeanP.printRowMatlab(" STATIC Neumann flux");

        // Time Dependent Dirichlet
        std::cout << "\n--| TIME DEPENDENT DIRICHLET:  ";
        std::cout << "\n # TIME DEPENDENT Dirichlet Bound ID: " << nDirTimeBound << "\n";
        dirTimeBound.printRowMatlab(" TIME DEPENDENDT Dirichlet Bound ID");
        for (int i = 0; i < nDirTimeBound; i++)
        {
            std::string tempName = std::to_string(i+1) + "th TIME DEPENDENT Dirichlet function";
            dirTimeFunc[i].print(tempName);
        }
        // Time Dependent Neumann
        std::cout << "\n--| TIME DEPENDENT NEUMANN:  ";
        std::cout << "\n # TIME DEPENDENT Neumann Bound ID: " << nNeuTimeBound << "\n";
        neuTimeBound.printRowMatlab(" TIME DEPENDENT Neumann Bound ID");
        for (int ibound = 0; ibound < nNeuTimeBound; ibound++)
        {
            std::string tempName = std::to_string(ibound+1) + "th TIME DEPENDENT Neumann function";
            neuTimeMeanP[ibound].print(tempName);
        }
        // NVOB
        std::cout << "\n--| NVOB:  ";
        std::cout << "\n # NVOB Bound ID: " << nNVOBBound << "\n";
        NVOBBound.printRowMatlab(" NVOB Bound ID");
        std::cout << "\n----------| END CHECK IMPORT PARAMETERS |----------|\n\n";
    }

    //---------------------------
    // CHECK PREPRO
    //---------------------------
    void checkPrepro()
    {
        // std::cout << "\n+++++++++++++++| NS CHECK IMPORT PREPRO |+++++++++++++++\n\n";
        // std::cout << " nNodes: " << (*physics).nNodes << "\n";
        // std::cout << " nElem: " << (*physics).nElem << "\n";
        // std::cout << " nBoundID: " << nBound << "\n";
        // std::cout << " nBoundNodes: " << nBoundNodes << "\n";
        // std::cout << " nBoundElem: " << nBoundElems << "\n";
        // for (int i = 0; i < nBound; i++)
        // {
            
        // }
    }

    //-------------------------
    // GET SOLUTION FROM COMSOL
    //-------------------------
    std::vector<MATRIX> getComSol(std::string fileSol)
    {
        std::ifstream solFile;
        solFile.open(fileSol);
        std::string line;
        std::istringstream iss;

        prec discard;
        int nthSol = (*physics).t_end/deltaT + 1;

        STREAM::getLines(solFile, line, 9);
        int nSol = (*physics).dim + 1; // dim for velocity, 1 for pressure

        std::vector<MATRIX> comSol(nthSol);
        for (int ithSol = 0; ithSol < nthSol; ithSol++) comSol[ithSol].initialize((*physics).nNodes, nSol); // initialize comSol;

        for (int inod = 0; inod < (*physics).nNodes; inod++)
        {
            getline(solFile, line);
            iss.str(line);
            iss >> discard; iss >> discard;
            for  (int ithSol = 0; ithSol < nthSol; ithSol++) //ithSol referes for ith-step solution
            {
                for (int jsol = 0; jsol < nSol; jsol++) // jsol refers to velocity components and pressure 
                {
                    iss >> comSol[ithSol][inod][jsol];
                }
            }
        }

        return comSol;
    }
    //---
    void compareStatSol(VECTOR &sol, VECTOR &solToCompare, VECTOR &error);
    VECTOR analyticStatError(VECTOR &sol, std::vector<std::string> &correctSol);

    MATRIX getVelocityFromSol_v(VECTOR &sol)
    {
        int nNodes_v = (*physics).nNodes_v;
        int dim = (*physics).dim;
        // GET VELOCITY
        MATRIX velocity(dim, nNodes_v);
        for (int icomp = 0; icomp < dim; icomp++)
        {
            for (int inode = 0; inode < nNodes_v; inode++)
            {
                int iglob = inode + nNodes_v*icomp;
                velocity[icomp][inode] = lastSol[iglob];
            }
        }
        return velocity;
    }
    //-------------------------------------------------------
    MATRIX getVelocityFromSol(VECTOR &sol)
    {
        int nNodes = (*physics).nNodes;
        int nNodes_v = (*physics).nNodes_v;
        int dim = (*physics).dim;
        // GET VELOCITY
        MATRIX velocity(nNodes, dim);
        for (int icomp = 0; icomp < dim; icomp++)
        {
            for (int inode = 0; inode < nNodes; inode++)
            {
                int iglob = inode + nNodes_v*icomp;
                velocity[inode][icomp] = sol[iglob];
            }
        }
        return velocity;
    }
     //-------------------------------------------------------
    void getUPFromSol_v(VECTOR &sol, MATRIX &velocity, VECTOR &P)
    {
        int nNodes = (*physics).nNodes;
        int nNodes_v = (*physics).nNodes_v;
        int dim = (*physics).dim;
        // GET VELOCITY
        for (int icomp = 0; icomp < dim; icomp++)
        {
            for (int inode = 0; inode < nNodes_v; inode++)
            {
                int iglob = inode + nNodes_v*icomp;
                velocity[inode][icomp] = lastSol[iglob];
            }
        }
        for (int inode = 0; inode < nNodes; inode++)
        {
            int iglob = inode + nNodes_v*dim;
            P[inode] = lastSol[iglob];
        }
    }

    void getUPFromSol(VECTOR &sol, MATRIX &velocity, VECTOR &P)
    {
        int nNodes = (*physics).nNodes;
        int nNodes_v = (*physics).nNodes_v;
        int dim = (*physics).dim;
        // GET VELOCITY
        for (int icomp = 0; icomp < dim; icomp++)
        {
            for (int inode = 0; inode < nNodes; inode++)
            {
                int iglob = inode + nNodes_v*icomp;
                velocity[inode][icomp] = lastSol[iglob];
            }
        }
        for (int inode = 0; inode < nNodes; inode++)
        {
            int iglob = inode + nNodes_v*dim;
            P[inode] = lastSol[iglob];
        }
    }

    //-------------------------
    // GET SOLUTION FROM COMSOL
    //-------------------------
    std::vector<MATRIX> getComSol(std::string fileSol, bool isStat);
    std::vector<MATRIX> prepSolToComSol(std::vector<VECTOR> &solToPrep, bool isStat);
    VECTOR compareStatComSol(MATRIX &sol, MATRIX &comSol);
    
    void createSysmatBase(CSRMAT &H, std::vector<CSRMAT> &B);

    //------------------------------------
    // GET LENGTH OR AREA
    //------------------------------------
    prec getArea(MATRIX &matCoord, int dim)
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

    prec static getSurface(MATRIX &matCoord, int dim)
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

    protected:
    //---------------------------
    // SET BOUNDARY CONDITIONS
    //---------------------------
    //---
    void setTimeDirBC(prec time);
    //---
    void setStatDirBC();
    //---
    void setTimeNeuBC(prec time);
    //---
    void setStatNeuBC(MATRIX_INT &boundInfoMat);
    //---
    void setSymmBC(MATRIX_INT &boundInfoMat);
    //---------------------------
    // ASSEMBLE
    //---------------------------
    void assembleN(VECTOR &U_r, CSRMAT& SYSMAT_final, VECTOR &rhs_final); // (CSRMAT &SYSMAT, VECTOR &rhs_final);
    //---------------------------
    // IMPOSE BOUNDARY CONDITIONS
    //---------------------------
    //---
    void imposeStaticBC(CSRMAT &SYSMAT, VECTOR& rhs);
    //---
    void imposeTimeBC(CSRMAT &SYSMAT, VECTOR& rhs);
    //---
    void updateBC(prec time);
    //---
    void updateSYSMAT();
    //---
    void updateTimeForcing();
    //---
    void updateRHS();
    //---

    

    //----------------------------------------------
    // COMPUTE EXTERNAL NORMAL TO A BOUNDARY ELEMENT
    //----------------------------------------------
    void getNormal(MATRIX &matCoord, VECTOR &normal)
    {
        switch ((*physics).dim)
        {
            case 2:
            {
                prec* v0 = matCoord[0]; prec* v1 = matCoord[1];
                VECTOR vec; VECTOR zVec;
                vec.setZeros(3); zVec.setZeros(3);
                for (int i = 0; i < 2; i++)
                {
                    vec[i] = v1[i] - v0[i];
                }
                zVec[2] = 1;
                VECTOR::cross(zVec, vec, normal);
                normal /= normal.norm();
                normal.shrink(2);
                normal.length = 2;
                break;
            }
            case 3:
            {
                prec* v0 = matCoord[0]; prec* v1 = matCoord[1]; prec* v2 = matCoord[2];
                VECTOR v0Vec(3); v0Vec = v0; VECTOR v1Vec(3); v1Vec = v1; VECTOR v2Vec(3); v2Vec = v2;
                v1Vec -= v0Vec; v2Vec -= v0Vec;
                VECTOR::cross(v1Vec, v2Vec, normal);
                prec norm = normal.norm();
                normal /= norm;
                break;
            }
        }
    } 
    // PRECONDITIONER

    void precondJacobiSolver(CSRMAT &SYS,VECTOR &rhs, VECTOR &sol);

    void precondSchurSolver(CSRMAT &SYS,VECTOR &rhs, VECTOR &sol);

    void SIMPLESolver(CSRMAT &SYS,VECTOR &rhs, VECTOR &sol, prec &finalRes);

    void SCHURgmres(CSRMAT &SYS,VECTOR &rhs, VECTOR &sol, prec &finalRes);

    void updatePrecond(CSRMAT &SYS, CSRMAT &A, VECTOR &M1, CSRMAT &B, CSRMAT &BT, CSRMAT &S, CSRMAT &leftPrecond, CSRMAT &rightPrecond, CSRMAT &invRightPrecond, CSRMAT &M1BT);
    // OTHER METHODS
    void getAFromSYSMAT(CSRMAT &SYS, CSRMAT &A)
    {
        int dim = (*physics).dim;
        int nNodes_v = (*physics).nNodes_v;
        int nTermA = (SYS.nTerm-2*nTermB*dim);

        A.initialize(dim*nNodes_v, dim*nNodes_v, nTermA);
        std::shared_ptr<int[]> iatA = A.iat; std::shared_ptr<int[]> jaA = A.ja; std::shared_ptr<prec[]> coefA = A.coef;

        std::shared_ptr<int[]> iat = SYS.iat; std::shared_ptr<int[]> ja = SYS.ja; std::shared_ptr<prec[]> coef = SYS.coef;
        int lastCol = dim*nNodes_v;

        int pos = 0;
        iatA[0] = 0;
        int irowA = 0;
        for (int irow = 0; irow < dim*nNodes_v; irow++)
        {
            for (int tempPos = iat[irow]; tempPos < iat[irow+1]; tempPos++)
            {
                int tempCol = ja[tempPos];
                if (tempCol >= lastCol) break;
                jaA[pos] = ja[tempPos];
                coefA[pos]  = coef[tempPos];
                pos++;
            }
            iatA[irowA+1] = pos;
            irowA++;
        }
    } 

    //-----------------

    void getBFromSYSMAT(CSRMAT &SYS, CSRMAT &B)
    {
        int dim = (*physics).dim;
        int nNodes = (*physics).nNodes;
        int nNodes_v = (*physics).nNodes_v;
        int nTerm = nTermB*dim;

        B.initialize(nNodes, dim*nNodes_v, nTerm);

        std::shared_ptr<int[]> iatB = B.iat; std::shared_ptr<int[]> jaB = B.ja; std::shared_ptr<prec[]> coefB = B.coef;
        std::shared_ptr<int[]> iat = SYS.iat; std::shared_ptr<int[]> ja = SYS.ja; std::shared_ptr<prec[]> coef = SYS.coef;
        int lastRow = dim*nNodes_v+nNodes;

        int pos = 0;
        iatB[0] = 0;
        int irowB = 0;
        
        for (int irow = dim*nNodes_v; irow < lastRow; irow++)
        {
            for (int tempPos = iat[irow]; tempPos < iat[irow+1]; tempPos++)
            {
                jaB[pos] = ja[tempPos];
                coefB[pos]  = coef[tempPos];
                pos++;
            }
            iatB[irowB+1] = pos;
            irowB++;
        }
    } 
};