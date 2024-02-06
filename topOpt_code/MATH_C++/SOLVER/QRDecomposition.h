#pragma once

#include "../../CODE_HEADERS/codeHeader.h"

class QR_DECOMPOSITION
{
    public:
    static void reduced(MATRIX &mat, MATRIX &Q, UPPER_TRIANGULAR &R)
    {
        int nCols = mat.nCol;
        int nRows = mat.nRow;
        if (nCols > nRows) throw_line("ERROR in QR factorization.\nMore columns than rows.\n");

        //-----------------------------
        // alloc something
        //-----------------------------
        VECTOR hv_j(nRows);
        VECTOR w(nRows); 
        //--------------------------------------------------
        for (int k = 0; k < nCols; k++)
        {
            std::shared_ptr<prec[]> temp = mat(k);
            w = temp;
            for (int j = 0; j < k; j++)
            {
                prec* tempP = &(Q(j)[0]);
                VECTOR pointColj(Q.nRow); pointColj = tempP;
                prec* pointColjP = &(pointColj[0]);

                // std::shared_ptr<prec[]> pointColj = Q(j);
                prec h = w.dot(pointColj);
                VECTOR::multiplyByScalar(pointColjP, h, nRows, hv_j);
                w -= hv_j;
                R.el(j,k,h);
            }
            prec hh = w.norm();
            R[k][0] = hh;
            w /= hh;
            
            for (int i = 0; i < nRows; i++)
            {
                Q[i][k] = w[i];
            }
        }
    }
    //---
    static void reduced(MATRIX_BY_COL &mat, MATRIX_BY_COL &Q, UPPER_TRIANGULAR &R)
    {
        int nCols = mat.nCol;
        int nRows = mat.nRow;
        if (nCols > nRows) throw_line("ERROR in QR factorization.\nMore columns than rows.\n");

        //-----------------------------
        // alloc something
        //-----------------------------
        VECTOR hv_j(nRows);
        VECTOR w(nRows); 
        //--------------------------------------------------
        for (int k = 0; k < nCols; k++)
        {
            prec* tempP2 = mat[k];
            w = tempP2;
            for (int j = 0; j < k; j++)
            {
                tempP2 = Q[j];
                VECTOR pointColj; pointColj = tempP2;
                prec h = w.dot(pointColj);
                VECTOR::multiplyByScalar(pointColj.P, h, nRows, hv_j);
                w -= hv_j;
                w.print();
                R.el(j,k,h);
            }
            prec hh = w.norm();
            R[k][0] = hh;
            w /= hh;
            prec* pointColk = Q[k];
            std::shared_ptr<prec[]> tempP(new(pointColk) prec[nRows], [](prec* pointColk){free(pointColk);});
            VECTOR::copy(w.P, tempP, nRows);
        }
    }
    //---
    static prec reduced(HESSEMBERG_BY_COL &mat, MATRIX_BY_COL &Q, UPPER_TRIANGULAR &R)
    {
        int nCols = mat.m;
        int nRows = mat.m + 1;

        prec res;

        //-----------------------------
        // alloc something
        //-----------------------------
        VECTOR hv_j(nRows);
        VECTOR w(nRows); 
        //--------------------------------------------------
        for (int k = 0; k < nCols; k++)
        {
            int maxRowID = 2+k;
            prec* tempRawP = mat[k];
            std::shared_ptr<prec[]> tempP(new(tempRawP) prec[maxRowID], [](prec* tempRawP){free(tempRawP);});
            VECTOR::copy(tempP, w.P, maxRowID);
            for (int j = 0; j < k; j++)
            {
                prec* tempP = Q[j];
                VECTOR pointColj; pointColj = tempP;
                prec h = VECTOR::dot(w.P, pointColj.P, maxRowID);
                VECTOR::multiplyByScalar(pointColj.P, h, nRows, hv_j);
                VECTOR::diff(w.P, hv_j.P, maxRowID, w.P);
                R.el(j,k,h);
            }
            prec hh = VECTOR::norm(w.P, maxRowID);
            R[k][0] = hh;
            VECTOR::divideByScalar(w.P, hh, maxRowID, w.P);
            prec* pointColk = Q[k];
            tempP.reset();
            tempP = std::shared_ptr<prec[]>(new(pointColk) prec[maxRowID], [](prec* pointColk){free(pointColk);});
            VECTOR::copy(w.P, tempP, maxRowID);
        }

        res = 0;
        for (int i = 0; i < nCols; i++) res += Q[i][0]*Q[i][0];
        res = sqrt(1-res);
        return res;
    }
    //---
    static prec forGMRES(HESSEMBERG_BY_COL &mat, UPPER_TRIANGULAR &R, prec* &d, int nCols)
    {
        int nRows = nCols + 1;
        MATRIX_BY_COL Q(nRows, nCols);
        if (d == nullptr)
        {
             d = (prec*) malloc(nCols * sizeof(prec));
            checkPtr(d);
        }

        //-----------------------------
        // alloc something
        //-----------------------------
        VECTOR hv_j(nRows);
        VECTOR w(nRows); 
        //--------------------------------------------------
        for (int k = 0; k < nCols; k++)
        {
            int maxRowID = 2+k;
            prec* tempRawP = mat[k];
            std::shared_ptr<prec[]> tempP(new(tempRawP) prec[maxRowID], [](prec* tempRawP){free(tempRawP);});
            VECTOR::copy(tempP, w.P, maxRowID);
            for (int j = 0; j < k; j++)
            {
                prec* tempP = Q[j];
                VECTOR pointColj; pointColj = tempP;
                prec h = VECTOR::dot(w.P, pointColj.P, maxRowID);
                VECTOR::multiplyByScalar(pointColj.P, h, nRows, hv_j);
                VECTOR::diff(w.P, hv_j.P, maxRowID, w.P);
                R.el(j,k,h);
            }
            prec hh = VECTOR::norm(w.P, maxRowID);
            R.el(k,k,hh);
            VECTOR::divideByScalar(w.P, hh, maxRowID, w.P);
            prec* pointColk = Q[k];
            tempP.reset();
            tempP = std::shared_ptr<prec[]>(new(pointColk) prec[maxRowID], [](prec* pointColk){free(pointColk);});
            VECTOR::copy(w.P, tempP, nRows);
            d[k] = Q[k][0];
        }
        prec res = 0.0;
        for (int i = 0; i < nCols; i++) res += d[i]*d[i];
        res = (1 - res);
        
        res = abs(res);
        if (abs(res) < 1e-15) res = 0; 
        res = sqrt(res);
        // R.print();
        return res;
    }
    //---

    static prec forGMRES(HESSEMBERG_BY_COL &mat, UPPER_TRIANGULAR &R, MATRIX_BY_COL& Q, std::shared_ptr<prec[]> &d, int nCols)
    {
        // int nRows = nCols + 1;

        //--------------------------------------------------
        int k = nCols - 1;       
        int maxRowID = 2+k;

        //-----------------------------
        // alloc something
        //-----------------------------
        VECTOR hv_j(maxRowID);
        VECTOR w(maxRowID); 
            
        prec* tempRawP = mat[k];
        std::shared_ptr<prec[]> tempP(new(tempRawP) prec[maxRowID], [](prec* tempRawP){free(tempRawP);});
        VECTOR::copy(tempP, w.P, maxRowID);
        for (int j = 0; j < k; j++)
        {
            int locMax = 2+j;
            prec* tempP = Q[j];
            VECTOR pointColj; pointColj = tempP;
            prec h = VECTOR::dot(w.P, pointColj.P, locMax);
            VECTOR::multiplyByScalar(pointColj.P, h, locMax, hv_j);
            VECTOR::diff(w.P, hv_j.P, locMax, w.P);
            R.el(j,k,h);
        }
        prec hh = VECTOR::norm(w.P, maxRowID);
        R.el(k,k,hh);
        VECTOR::divideByScalar(w.P, hh, maxRowID, w.P);
        prec* pointColk = Q[k];
        tempP.reset();
        tempP = std::shared_ptr<prec[]>(new(pointColk) prec[maxRowID], [](prec* pointColk){free(pointColk);});
        VECTOR::copy(w.P, tempP, maxRowID);
        d[k] = pointColk[0];
        prec res = 0.0;
        for (int i = 0; i < nCols; i++) res += d[i]*d[i];

        res = (1 - res);
        //res = abs(res);
        if (res < 0 && res > -1e-16) res = std::abs(res);
        if (res < 0)
        {
            return res;
        }
        res = sqrt(res);
        return res;
    }
    //---
    //---

    static prec forGMRES(HESSEMBERG_BY_COL &mat, UPPER_TRIANGULAR &R, HESSEMBERG_BY_COL& Q, std::shared_ptr<prec[]> &d, int nCols)
    {
        //--------------------------------------------------
        int k = nCols - 1;       
        int maxRowID = 2+k;

        //-----------------------------
        // alloc something
        //-----------------------------
        VECTOR hv_j(maxRowID);
        VECTOR w(maxRowID); 
        
        prec* tempRawP = mat[k];
        VECTOR tempP(maxRowID); tempP = tempRawP;
        VECTOR::copy(tempP.P, w.P, maxRowID);
            for (int j = 0; j < k; j++)
            {
                int locMax = 2+j;
                prec* tempP = Q[j];
                VECTOR pointColj(locMax); pointColj = tempP;
                prec h = VECTOR::dot(w.P, pointColj.P, locMax);
                VECTOR::multiplyByScalar(pointColj.P, h, locMax, hv_j);
                VECTOR::diff(w.P, hv_j.P, locMax, w.P);
                R.el(j,k,h);
            }
            prec hh = VECTOR::norm(w.P, maxRowID);
            R.el(k,k,hh);
            VECTOR::divideByScalar(w.P, hh, maxRowID, w.P);
            prec* pointColk = Q[k];
            VECTOR::copy(w.P, pointColk, maxRowID);
            d[k] = pointColk[0];
        prec res = 0.0;
        for (int i = 0; i < nCols; i++) res += d[i]*d[i];
        
        res = (1 - res);
        
        //res = abs(res);
        //if (res < 0 && res > -1e-16) res = std::abs(res);
        if (res < 0)
        {
            return res;
        }
        res = sqrt(res);
        return res;
    }
    //-------------------------------------------------------------------------------
    // HOUSEHOLDER PROCESS
    //-------------------------------------------------------------------------------
};
