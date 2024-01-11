#pragma once

#include "../../CODE_HEADERS/codeHeader.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

/* PARDISO prototype. */
extern "C" void pardisoinit (void   *, int    *,   int *, int *, double *, int *);
extern "C" void pardiso     (void   *, int    *,   int *, int *,    int *, int *, 
                  double *, int    *,    int *, int *,   int *, int *,
                     int *, double *, double *, int *, double *);
extern "C" void pardiso_chkmatrix  (int *, int *, double *, int *, int *, int *);
extern "C" void pardiso_chkvec     (int *, int *, double *, int *);
extern "C" void pardiso_printstats (int *, int *, double *, int *, int *, int *, double *, int *);

class SOLVERLS
{
public:
    //----------------------------------------------
    // STANDARD PARDISO
    //----------------------------------------------
    static void launchPardiso(CSRMAT &mat, std::shared_ptr<prec[]> &b, std::shared_ptr<prec[]> &sol, int type = 1)
    {
        std::string licMsg = "PARDISOLICMESSAGE=1";
        putenv(&licMsg[0]);
        /* Matrix data. */
        int      n = mat.nRow;
        std::shared_ptr<int[]>   iat = mat.iat;
        std::shared_ptr<int[]>    ja = mat.ja;
        std::shared_ptr<prec[]> coef = mat.coef;

        int      nnz = iat[n];
        int      mtype = 11;        /* Real unsymmetric matrix */

        /* RHS and solution vectors. */
        int      nrhs = 1;          /* Number of right hand sides. */
        if (sol == 0) std::shared_ptr<prec[]> sol(new prec[n]);

        /* Internal solver memory pointer pt,                  */
        /* 32-bit: int pt[64]; 64-bit: long int pt[64]         */
        /* or void *pt[64] should be OK on both architectures  */ 
        void    *pt[64];

        /* Pardiso control parameters. */
        int      iparm[64];
        double   dparm[64];
        int      solver;
        int      maxfct, mnum, phase, error, msglvl;

        /* Number of processors. */
        int      num_procs;

        /* Auxiliary variables. */
        int      i;

        double   ddum;              /* Double dummy */
        int      idum;              /* Integer dummy. */

    /* -------------------------------------------------------------------- */
    /* ..  Setup Pardiso control parameters and initialize the solvers      */
    /*     internal adress pointers. This is only necessary for the FIRST   */
    /*     call of the PARDISO solver.                                      */
    /* ---------------------------------------------------------------------*/
        error = 0;
        solver = type; /* 0 use sparse direct solver */
        pardisoinit (pt,  &mtype, &solver, iparm, dparm, &error);
        if (error != 0)
        {
            if (error == -10 )
            throw_line("ERROR(PARDISO) No license file found \n");
            if (error == -11 )
            throw_line("ERROR(PARDISO) License is expired \n");
            if (error == -12 )
            throw_line("ERROR(PARDISO) Wrong username or hostname \n");
        }
        
        num_procs = omp_get_max_threads();
        iparm[2]  = num_procs;
    
        
        maxfct = num_procs;         /* Maximum number of numerical factorizations.  */
        mnum   = 1;         /* Which factorization to use. */
        
        msglvl = 0;         /* Print statistical information  */
        error  = 0;         /* Initialize error flag */


    /* -------------------------------------------------------------------- */    
    /* ..  Convert matrix from 0-based C-notation to Fortran 1-based        */
    /*     notation.                                                        */
    /* -------------------------------------------------------------------- */ 
        for (i = 0; i < n+1; i++) {
            iat[i] += 1;
        }
        for (i = 0; i < nnz; i++) {
            ja[i] += 1;
        }
    /* -------------------------------------------------------------------- */
    /*  .. pardiso_chk_matrix(...)                                          */
    /*     Checks the consistency of the given matrix.                      */
    /*     Use this functionality only for debugging purposes               */
    /* -------------------------------------------------------------------- */
        
        // pardiso_chkmatrix  (&mtype, &n, coef, iat, ja, &error);
        // if (error != 0) {
        //     printf("\nERROR in consistency of matrix: %d", error);
        //     exit(1);
        // }

    /* -------------------------------------------------------------------- */
    /* ..  pardiso_chkvec(...)                                              */
    /*     Checks the given vectors for infinite and NaN values             */
    /*     Input parameters (see PARDISO user manual for a description):    */
    /*     Use this functionality only for debugging purposes               */
    /* -------------------------------------------------------------------- */

        pardiso_chkvec (&n, &nrhs, &(b[0]), &error);
        if (error != 0) {
            printf("\nERROR  in right hand side: %d", error);
            exit(1);
        }
    /* -------------------------------------------------------------------- */
    /* .. pardiso_printstats(...)                                           */
    /*    prints information on the matrix to STDOUT.                       */
    /*    Use this functionality only for debugging purposes                */
    /* -------------------------------------------------------------------- */

        // pardiso_printstats (&mtype, &n, coef, iat, ja, &nrhs, b, &error);
        // if (error != 0) {
        //     printf("\nERROR right hand side: %d", error);
        //     exit(1);
        // }
    
    /* -------------------------------------------------------------------- */    
    /* ..  Reordering and Symbolic Factorization.  This step also allocates */
    /*     all memory that is necessary for the factorization.              */
    /* -------------------------------------------------------------------- */ 
        phase = 11; 

        pardiso (pt, &maxfct, &mnum, &mtype, &phase,
                &n, &(coef[0]), &(iat[0]), &(ja[0]), &idum, &nrhs,
                iparm, &msglvl, &ddum, &ddum, &error,  dparm);
        if (error != 0) {
            printf("\nERROR during symbolic factorization: %d", error);
            exit(1);
        }
        // printf("\nReordering completed ... ");
        // printf("\nNumber of nonzeros in factors  = %d", iparm[17]);
        // printf("\nNumber of factorization MFLOPS = %d", iparm[18]);
    
    /* -------------------------------------------------------------------- */    
    /* ..  Numerical factorization.                                         */
    /* -------------------------------------------------------------------- */    
        phase = 22;

        pardiso (pt, &maxfct, &mnum, &mtype, &phase,
                &n, &(coef[0]), &(iat[0]), &(ja[0]), &idum, &nrhs,
                iparm, &msglvl, &ddum, &ddum, &error, dparm);
        if (error != 0) {
            printf("\nERROR during numerical factorization: %d", error);
            exit(2);
        }
        printf("\nFactorization completed ...\n ");

    /* -------------------------------------------------------------------- */    
    /* ..  Back substitution and iterative refinement.                      */
    /* -------------------------------------------------------------------- */    
        phase = 33;

        iparm[7] = 2;       /* Max numbers of iterative refinement steps. */

    
        pardiso (pt, &maxfct, &mnum, &mtype, &phase,
                &n, &(coef[0]), &(iat[0]), &(ja[0]), &idum, &nrhs,
                iparm, &msglvl, &(b[0]), &(sol[0]), &error,  dparm);
    
        if (error != 0) {
            printf("\nERROR during solution: %d", error);
            exit(3);
        }

    /* -------------------------------------------------------------------- */
    /* ..  Back substitution with tranposed matrix A^t x=b                   */
    /* -------------------------------------------------------------------- */
        phase = 33;

        iparm[7]  = 2;       /* Max numbers of iterative refinement steps. */
        iparm[11] = 1;       /* Solving with transpose matrix. */

    
        pardiso (pt, &maxfct, &mnum, &mtype, &phase,
                &n, &(coef[0]), &(iat[0]), &(ja[0]), &idum, &nrhs,
                iparm, &msglvl, &(b[0]), &(sol[0]), &error,  dparm);
    
        if (error != 0) {
            printf("\nERROR during solution: %d", error);
            exit(3);
        }

    /* -------------------------------------------------------------------- */    
    /* ..  Convert matrix back to 0-based C-notation.                       */
    /* -------------------------------------------------------------------- */ 
        for (i = 0; i < n+1; i++) {
            iat[i] -= 1;
        }
        for (i = 0; i < nnz; i++) {
            ja[i] -= 1;
        }

    /* -------------------------------------------------------------------- */    
    /* ..  Termination and release of memory.                               */
    /* -------------------------------------------------------------------- */ 
        phase = -1;                 /* Release internal memory. */

        pardiso (pt, &maxfct, &mnum, &mtype, &phase,
                &n, &ddum, &(iat[0]), &(ja[0]), &idum, &nrhs,
                iparm, &msglvl, &ddum, &ddum, &error,  dparm);
    }
    //---------------------------------------------------------------------------------------------------------------
    // GMRES 
    //---------------------------------------------------------------------------------------------------------------
    static void gmres_p(CSRMAT &mat, std::shared_ptr<prec[]> &b, std::shared_ptr<prec[]> &solOut, std::shared_ptr<prec[]> &x_0, prec bNorm, prec toll = 1e-16, int p = 100, int maxIt = 1e8)
    {
        int N = mat.nRow;

        // check
        if (maxIt > N) maxIt = N;
        if (p > N) p = N; 
        if (N != mat.nCol) throw_line("ERROR, in matrix dimension, non square matrix\n");

        //double bNorm = VECTOR::norm(b,N);
        VECTOR sol; sol.length = N; sol.P = solOut;
        sol = x_0;
        //---------------------------------------------------
        // alloc of everything useful
        //---------------------------------------------------
        VECTOR w(N); // new vector linear independent
        VECTOR d(p+1); // Q matrix (in QR decomposition) first row
        VECTOR hv_j(N); // store temporary product h_{i,j}*V[j];
        VECTOR temp(N); // store temporary vector result
        VECTOR z(p); // store results of Rz = d

        HESSEMBERG_BY_COL H(p+1);    
        MATRIX_BY_COL V(N,p+1);

        // alloc for QR factorization
        UPPER_TRIANGULAR R(p+1);
        HESSEMBERG_BY_COL Q(p+1);
        
        //----------------------------------------
        // INITIALIZE PARAMETERS
        //----------------------------------------
        std::shared_ptr<prec[]> lastV;  //pointer to last vector of V
        std::shared_ptr<prec[]> Ax_old = mat.prod(x_0, N);
        std::shared_ptr<prec[]> r_old  = 0;
        VECTOR::diff(b, Ax_old, N, r_old);
        prec  beta   = VECTOR::norm(r_old, N);

        int it = 0;
        int k = 0;
        prec res;
        //------------------------------
        while (beta/bNorm > toll && it < maxIt)
        {
            prec* colAddress = V[0];
            VECTOR::divideByScalar(r_old, beta, N, colAddress);
            maxIt += k;
            k = 0;
            res = toll + 1;
            //---------------------------------------
            while (res > toll && k < p && it < maxIt)
            {
                prec* wP = &(w.P[0]);
                
                mat.prod(V[k], N, wP);
                for (int j = 0; j < k+1; j++)
                {
                    prec* tempP = V[j];
                    prec h = VECTOR::dot(w.P, tempP, N); // scalar products
                    VECTOR::multiplyByScalar(tempP, h, N, hv_j);
                    VECTOR::diff(w.P, hv_j.P, N, w.P);

                    H[k][j] = h;
                }
                prec hh = VECTOR::norm(w.P, N); // h_{k+1, k}
                prec* tempP = V[k+1];
                VECTOR::divideByScalar(w.P, hh, N, tempP); //store V last vector
                H[k][k+1] = hh;
                
                res = QR_DECOMPOSITION::forGMRES(H, R, Q, d.P, k+1);

                //-----
                if (res < 0) 
                {
                    printf("it: %d res: %" format_e "  WARNING: NEGATIVE RESIDUAL, PROBABLY DUE TO MALCONDITIONATE MATRIX.\nRESTARTED GMRES APPLIED.\n", it, res);
                    res = toll + 1;
                    break;
                }
                //-----
                res = beta*res/bNorm;
                k++;
                it++;
            }
            //-----------------------------------------------
            VECTOR::multiplyByScalar(d.P, beta, k, temp); // temp: known term of Rz = temp (d)
            R.solveLS(temp.P, z.P, k); // find z sol

            // compute solution
            for (int i = 0; i < k; i++)
            {
                prec* tempP = V[i];
                VECTOR::multiplyByScalar(tempP, z[i], N, temp.P);
                VECTOR::sum(sol.P, temp.P, N, sol.P);
                // sol.print();
            }
            std::shared_ptr<prec[]> Asol = mat.prod(sol.P, N);
            r_old = nullptr;
            VECTOR::diff(b, Asol, N, r_old);
            beta = VECTOR::norm(r_old, N); // actual residuals

            if (beta/bNorm > toll) VECTOR::copy(sol.P, x_0, N); 
        }
        // //-------------------------------------------
        printf("final it: %d relative res: %" format_e "\n", it, beta/bNorm);

    }
    //---
    //---------------------------------------------
    // JACOBI PRECONDITIONED GMRES
    //---------------------------------------------
    static std::shared_ptr<prec[]> precGmres_p(CSRMAT &mat, std::shared_ptr<prec[]> b, std::shared_ptr<prec[]> x_0, prec toll = 1e-20, int p = 1e2, int maxIt = 1e8)
    {
        int N = mat.nRow;

        // check
        if (maxIt > N) maxIt = N;
        if (p > N) p = N; 
        if (N != mat.nCol) throw_line("ERROR, in matrix dimension, non square matrix\n");

        //double bNorm = VECTOR::norm(b,N);

        VECTOR sol; sol = x_0;
        //---------------------------------------------------
        // alloc of everything useful
        //---------------------------------------------------
        VECTOR w(N); // new vector linear independent
        VECTOR d(p+1); // Q matrix (in QR decomposition) first row
        VECTOR hv_j(N); // store temporary product h_{i,j}*V[j];
        VECTOR temp(N); // store temporary vector result
        VECTOR z(p); // store results of Rz = d

        HESSEMBERG_BY_COL H(p+1);    
        MATRIX_BY_COL V(N,p+1);

        // alloc for QR factorization
        UPPER_TRIANGULAR R(p+1);
        HESSEMBERG_BY_COL Q(p+1);
        
        //----------------------------------------
        // PRECOND MATRIX 
        //----------------------------------------
        std::shared_ptr<prec[]> diag = mat.diag();
        mat.jacobiPrecond(diag);
        //----------------------------------------
        // INITIALIZE PARAMETERS
        //----------------------------------------
        std::shared_ptr<prec[]> lastV;  //pointer to last vector of V
        std::shared_ptr<prec[]> Ax_old = mat.prod(x_0, N);
        std::shared_ptr<prec[]> r_old  = nullptr;
        VECTOR::diff(b, Ax_old, N, r_old);
        prec  beta   = VECTOR::norm(r_old, N);
        prec bNorm   = VECTOR::norm(b, N);

        int it = 0;
        int k = 0;
        prec res;
        //------------------------------
        while (beta/bNorm > toll && it < maxIt)
        {
            prec* colAddress = V[0];
            VECTOR::divideByScalar(r_old, beta, N, colAddress);
            maxIt += k;
            k = 0;
            res = toll + 1;
            //---------------------------------------
            while (res > toll && k < p && it < maxIt)
            {
                prec* p1 = V[k];
                prec* p2 = &(w[0]);
                mat.prod(p1, N, p2);
                for (int j = 0; j < k+1; j++)
                {
                    prec* tempP2 = V[j];
                    std::shared_ptr<prec[]> tempVec2 (new(tempP2) prec[N], [](prec* tempP2){free(tempP2);});
                    prec h = VECTOR::dot(w.P, tempVec2, N); // scalar products
                    VECTOR::multiplyByScalar(tempVec2, h, N, hv_j);
                    VECTOR::diff(w.P, hv_j.P, N, w.P);

                    H[k][j] = h;
                }
                prec hh = VECTOR::norm(w.P, N); // h_{k+1, k}
                prec* tempP = V[k+1];
                lastV = std::shared_ptr<prec[]> (new(tempP) prec[N], [](prec* tempP){free(tempP);});
                VECTOR::divideByScalar(w.P, hh, N, lastV); //store V last vector
                H[k][k+1] = hh;
                
                res = QR_DECOMPOSITION::forGMRES(H, R, Q, d.P, k+1);
                //-----
                if (res < 0) 
                {
                    printf("it: %d res: %" format_e "  WARNING: NEGATIVE RESIDUAL, PROBABLY DUE TO MALCONDITIONATE MATRIX.\nRESTARTED GMRES APPLIED.\n", it, res);
                    res = toll + 1;
                    break;
                }
                //-----
                res = beta*res/bNorm;
                k++;
                it++;
            }
            //-----------------------------------------------
            VECTOR::multiplyByScalar(d.P, beta, k, temp.P); // temp: known term of Rz = temp (d)
            R.solveLS(temp.P, z.P, k); // find z sol

            // compute solution
            for (int i = 0; i < k; i++)
            {
                prec* tempP = V[i];
                std::shared_ptr<prec[]> tempVec (new(tempP) prec[N], [](prec* tempP){free(tempP);});
                VECTOR::multiplyByScalar(tempVec, z[i], N, temp);
                VECTOR::sum(sol.P, temp.P, N, sol.P);
            }
            std::shared_ptr<prec[]> Asol = mat.prod(sol.P, N);
            r_old = nullptr;
            VECTOR::diff(b, Asol, N, r_old);
            beta = VECTOR::norm(r_old, N); // actual residuals
            if (beta/bNorm > toll) VECTOR::copy(sol.P, x_0, N); 
        }
        // //-------------------------------------------
        VECTOR::pointdiv(sol.P, diag, N, sol.P);
        // printf("final it: %d res: %" format_e "\n", it, res);
        printf("final it: %d relative res: %" format_e "\n", it, beta/bNorm);
        return sol.P;
    }
    //---
    //----------------------------------------------
    // JACOBI FOR NAVIER-STOKES PRECONDITIONED GMRES
    //----------------------------------------------
    //---------------------------------------------------
    // GMRES PRECONDITIONED WITH JACOBI FOR NAVIER-STOKES
    //---------------------------------------------------
    static void precJacobiNSGmres_p(CSRMAT &mat, std::shared_ptr<prec[]> &b, VECTOR &sol, std::shared_ptr<prec[]> &x_0, prec bNorm, prec toll = 1e-16, int p = 100, int maxIt = 1e8)
    {
        int N = mat.nRow;

        // check
        if (maxIt > N) maxIt = N;
        if (p > N) p = N; 
        if (N != mat.nCol) throw_line("ERROR, in matrix dimension, non square matrix\n");

        //double bNorm = VECTOR::norm(b,N);

        sol = x_0;

        //----------------------------------------
        // PRECOND MATRIX 
        //----------------------------------------
        VECTOR diag; diag.length = N;
        diag.P = mat.precondNSDiag();
        mat.jacobiNSPrecond(diag.P);
        //---------------------------------------------------
        // alloc of everything useful
        //---------------------------------------------------
        VECTOR w(N); // new vector linear independent
        VECTOR d(p+1); // Q matrix (in QR decomposition) first row
        VECTOR hv_j(N); // store temporary product h_{i,j}*V[j];
        VECTOR temp(N); // store temporary vector result
        VECTOR z(p); // store results of Rz = d

        HESSEMBERG_BY_COL H(p+1);    
        MATRIX_BY_COL V(N,p+1);

        // alloc for QR factorization
        UPPER_TRIANGULAR R(p+1);
        HESSEMBERG_BY_COL Q(p+1);
        
        //----------------------------------------
        // INITIALIZE PARAMETERS
        //----------------------------------------
        std::shared_ptr<prec[]> lastV;  //pointer to last vector of V
        std::shared_ptr<prec[]> Ax_old = mat.prod(x_0, N);
        std::shared_ptr<prec[]> r_old  = 0;
        VECTOR::diff(b, Ax_old, N, r_old);
        prec  beta   = VECTOR::norm(r_old, N);

        int it = 0;
        int k = 0;
        prec res;
        //------------------------------
        while (beta/bNorm > toll && it < maxIt)
        {
            prec* colAddress = V[0];
            VECTOR::divideByScalar(r_old, beta, N, colAddress);
            maxIt += k;
            k = 0;
            res = toll + 1;
            //---------------------------------------
            while (res > toll && k < p && it < maxIt)
            {
                prec* wP = &(w.P[0]);
                
                mat.prod(V[k], N, wP);
                for (int j = 0; j < k+1; j++)
                {
                    prec* tempP = V[j];
                    prec h = VECTOR::dot(w.P, tempP, N); // scalar products
                    VECTOR::multiplyByScalar(tempP, h, N, hv_j);
                    VECTOR::diff(w.P, hv_j.P, N, w.P);

                    H[k][j] = h;
                }
                prec hh = VECTOR::norm(w.P, N); // h_{k+1, k}
                prec* tempP = V[k+1];
                VECTOR::divideByScalar(w.P, hh, N, tempP); //store V last vector
                H[k][k+1] = hh;
                
                res = QR_DECOMPOSITION::forGMRES(H, R, Q, d.P, k+1);

                //-----
                if (res < 0) 
                {
                    printf("it: %d res: %" format_e "  WARNING: NEGATIVE RESIDUAL, PROBABLY DUE TO MALCONDITIONATE MATRIX.\nRESTARTED GMRES APPLIED.\n", it, res);
                    res = toll + 1;
                    break;
                }
                //-----
                res = beta*res/bNorm;
                k++;
                it++;
            }
            //-----------------------------------------------
            VECTOR::multiplyByScalar(d.P, beta, k, temp); // temp: known term of Rz = temp (d)
            R.solveLS(temp.P, z.P, k); // find z sol

            // compute solution
            for (int i = 0; i < k; i++)
            {
                prec* tempP = V[i];
                VECTOR::multiplyByScalar(tempP, z[i], N, temp.P);
                VECTOR::sum(sol.P, temp.P, N, sol.P);
                // sol.print();
            }
            std::shared_ptr<prec[]> Asol = mat.prod(sol.P, N);
            r_old = nullptr;
            VECTOR::diff(b, Asol, N, r_old);
            beta = VECTOR::norm(r_old, N); // actual residuals
            if (beta/bNorm > toll) VECTOR::copy(sol.P, x_0, N); 
        }
        // //-------------------------------------------
        printf("final it: %d relative res: %" format_e "\n", it, beta/bNorm);

        sol.pointdiv(diag, sol);

        // free pointers
        // VECTOR::print(b, N);
    }
};

