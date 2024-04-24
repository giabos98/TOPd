#pragma once

#include "CODE_HEADERS/codeHeader.h"
#include "geometry.h"
//#include <sys/time.h>

class PROBLEM_DARCY
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
    int n_domains;

    VECTOR alpha;

    //--- ABOUT TIME ---
    prec time;
    prec deltaT;
    int globIter = 0;
    VECTOR real_solution_times;

    //---  PROBLEM_DARCY ---
    VECTOR domains_permeability;
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
    // int dimTerms; // number of terms for each components in the SYSMAT matrix = H.nTerm + B[i].nTerm
    VECTOR_INT iSparse; VECTOR_INT jSparse; VECTOR_INT realPos; 
    std::vector<std::vector<VECTOR_INT>> realPosSysmat; // std vector containing in each comp the pattern of H coeff in sysmatNS
    // CSRMAT H; CSRMAT M; std::vector<CSRMAT> B;
    CSRMAT H; 
    CSRMAT M;
    VECTOR N; 
    VECTOR Ma;
    // std::vector<CSRMAT> B;
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
    std::shared_ptr<int[]> nBoundIdNodes; //# of boundary nodes of a specific ID

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
    PROBLEM_DARCY() {}
    PROBLEM_DARCY(PHYSICS* &Physics, std::string probRefFile)
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
    PROBLEM_DARCY(int dimension)
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

    //-----------------------------------------
    // IMPORT PRE-PRO
    //-----------------------------------------
    void importPREPRO();

    //-----------------------------------------
    // INITIALIZE PROBLEM_DARCY PARAMETERS
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
    void StatSolver();
    //----
    void StatSolverIterative();
    //----
    // void oneStepSolverStokes();
    // void oneStepSolverNavierStokes(prec toll = 1e-2, int itMax = 5);
    void oneStepSolver(prec toll = 1e-2, int itMax = 5);

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
    void checkImportParameters();

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
};