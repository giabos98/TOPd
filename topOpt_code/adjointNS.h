#pragma once

#include "CODE_HEADERS/codeHeader.h"
#include "nsProblem.h"

class ADJOINT_NS : public PROBLEM_NS
{
public: 
    //--- BC ---
    std::shared_ptr<int[]> normalBoundIdElems_buff_buff;
    std::shared_ptr<int*[]> normalBoundIdElems_buff;
    std::shared_ptr<int**[]> normalBoundIdElems;
    std::shared_ptr<int[]> nNormalBoundIdElems;

    std::shared_ptr<int[]> NVOBBoundIdElems_buff_buff;
    std::shared_ptr<int*[]> NVOBBoundIdElems_buff;
    std::shared_ptr<int**[]> NVOBBoundIdElems;
    std::shared_ptr<int[]> nNVOBBoundIdElems;
    VECTOR_INT dirPosToZero; VECTOR_INT dirIglobToZero; VECTOR_INT dirCountToZero; VECTOR_INT dirPosToOne;
    MATRIX NVOBVal;

    //--- SYSMAT ---------
    VECTOR Ac; 
    std::vector<std::vector<VECTOR>> Ag;

    int customFunc;
    int onlyGrad;
    VECTOR fWeights;

    MATRIX dirNormals;

    CSRMAT dAdGu;

    //time iteration count
    int curr_iter = 0;

    ADJOINT_NS() : PROBLEM_NS(){};

    //-----------------------------
    // INITIALIZE FROM NS PROBLEM
    //-----------------------------
    void initialize(PROBLEM_NS &NS, VECTOR &alphaIn);

    void setBC(PROBLEM_NS &NS);

    void setStatNeuBC();
    void setStatDirBC();

    void assembleA(VECTOR &velocity);

    // void createSysmatBase(CSRMAT H, std::vector<CSRMAT> B, VECTOR P_bc);
    // void addToSysmat(std::vector<std::vector<VECTOR>> &coef, int selectPath, prec factor);
    void imposeBC(CSRMAT &SYSMAT_final, VECTOR &rhs);
    void imposeStaticBC(CSRMAT &SYSMAT_NS, VECTOR &rhs);
    void imposeTimeBC(CSRMAT &SYSMAT_NS, VECTOR &rhs);

    void prepareSolver(PROBLEM_NS &NS);
    void updateRHS(VECTOR &nsSol);
    void addToSysmat(std::vector<std::vector<VECTOR>> &coef, int selectPath = 0, prec factor = 1);
    void applyStatForcing();

    void createSysmatBase(CSRMAT &H, std::vector<CSRMAT> &B);

    void updateSYSMAT(VECTOR &sol);

    void assembledAdgu();
    //---------------------------------------------------
    // SOLVER
    //---------------------------------------------------
    // void Solver(VECTOR solNS)
    // {
    //     // pos dal primo BC
    //     prec t_end = (*physics).t_end;
    //     int nTimeSteps = t_end/deltaT;

    //     while ((t_end-time) > (1e-15 * (nTimeSteps+1)))
    //     {
    //         oneStepSolver(solNS);
    //     }
    //     VTKWriter.closeTFile();
    // }
    void StatSolver(VECTOR &nsSol);
    void Solver();

    void oneStepSolverAdjoint(VECTOR &nsSol);

};