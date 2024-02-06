#pragma once

#include "../CODE_HEADERS/codeHeader.h"

//------------------------------------------------
//---
// PARALLEL COMPUTING
//---
//------------------------------------------------

class PARALLEL
{
public:
    static int nThread;

    //------------------------------------------------
    // VECTOR_INT
    //------------------------------------------------
    static void minus(VECTOR &vec);
    static void sum(VECTOR_INT &v1, VECTOR_INT &v2, VECTOR_INT &res);
    static void diff(VECTOR_INT &v1, VECTOR_INT &v2, VECTOR_INT &res);
    static void copy(VECTOR_INT &vecIn, VECTOR_INT &vecOut);
    static void copy(int *&vecIn, int &N, int *&vecOut);
    static void resetZeros(VECTOR_INT &vec);

    //------------------------------------------------
    // VECTOR
    //------------------------------------------------
    static void set(VECTOR &vec, prec coef);
    static void sum(VECTOR &v1, VECTOR &v2, VECTOR &res);  static void sum(prec* &v1, prec* &v2, int &N, prec* &res);
    static void diff(VECTOR &v1, VECTOR &v2, VECTOR &res); static void diff(prec* &v1, prec* &v2, int &N, prec* &res);
    static void pointProd(VECTOR &v1, VECTOR &v2, VECTOR &res); static void pointProd(prec* &v1, prec* &v2, int &N, prec* &res);
    static void pointDiv(VECTOR &v1, VECTOR &v2, VECTOR &res); static void pointDiv(prec* &v1, prec* &v2, int &N, prec* &res);
    static void multiplyByScalar(VECTOR &v1, prec &coef, VECTOR &res); static void multiplyByScalar(prec* &v1, prec &coef, int &N, prec* &res);
    static void divideByScalar(VECTOR &v1, prec &coef, VECTOR &res); static void divideByScalar(prec* &v1, prec &coef, int &N, prec* &res);
    static prec dot(VECTOR &v1, VECTOR &v2); static prec dot(prec* &v1, prec* &v2, int &N);
    static prec norm(VECTOR &vec); static prec norm(prec* &vec, int &N);
    static void copy(VECTOR &vecIn, VECTOR &vecOut); static void copy(prec* &vecIn, int &N, prec* &vecOut);
    static void resetZeros(VECTOR &vec); static void resetZeros(prec* &vec, int &N);

    //------------------------------------------------
    // CSRMAT
    //------------------------------------------------
    static void minus(CSRMAT &mat);
    static void prod(CSRMAT &mat, prec &coef, CSRMAT &res); 
    static void prod(CSRMAT &mat, VECTOR &vec, VECTOR &res); static void prod(CSRMAT &mat, prec* &vec, int &N, prec* &res);   
    static void prod(CSRMAT &mat1, CSRMAT &mat2, CSRMAT &res);
    static void diag(CSRMAT &mat, VECTOR & res);
    static void RprodByDiag(CSRMAT &mat, VECTOR &diag, CSRMAT &res);
    static void LprodByDiag(CSRMAT &mat, VECTOR &diag, CSRMAT &res);
    static void getTranspose(CSRMAT &mat, CSRMAT &res);
    static void copyPatt(CSRMAT &matIn, CSRMAT &matOut);
    static void concatInRow(CSRMAT &mat1, CSRMAT &mat2, CSRMAT &res);
    static void concatInCol(CSRMAT &mat1, CSRMAT &mat2, CSRMAT &res, int delay = 0);
    static void createBlockMatrix(CSRMAT &refMat, std::vector<std::vector<VECTOR*>> &coef, CSRMAT &matRes);
    static void createBlockMatrix(std::vector<std::vector<CSRMAT*>> &refMat, CSRMAT &matRes);


    // SOLVER
    static void solveLS(UPPER_TRIANGULAR &mat, VECTOR &b, VECTOR &sol, int &N);
    static void gmres_p(CSRMAT &mat, VECTOR &b, VECTOR &sol, VECTOR &x_0, prec &bNorm, prec toll = 1e-16, int p = 100, int maxIt = 1e8);

    static void SIMPLEgmres_p(CSRMAT &mat, CSRMAT &A, CSRMAT &B, CSRMAT &BT, CSRMAT &S, CSRMAT &M1BT,
                              VECTOR &b, VECTOR &sol, VECTOR &x_0, prec &bNorm, prec& finalRes, prec toll = 1e-16, int p = 100, int maxIt = 1e3);

    static void SCHURgmres_p(CSRMAT &mat, CSRMAT &M, CSRMAT &MR, VECTOR &b, VECTOR &sol, VECTOR &x_0, prec &bNorm, prec &finalRes, prec toll = 1e-16, int p = 100, int maxIt = 1e3);

    //QR FACTORIZATION
    static prec QRforGMRES(HESSEMBERG_BY_COL &mat, UPPER_TRIANGULAR &R, HESSEMBERG_BY_COL& Q, prec* &d, int nCols);

};