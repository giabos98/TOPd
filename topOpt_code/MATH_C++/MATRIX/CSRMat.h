#pragma once

#include "../../CODE_HEADERS/codeHeader.h"

class CSRMAT
{
    public: 
    std::shared_ptr<int[]> iat = 0;
    std::shared_ptr<int[]> ja = 0;
    std::shared_ptr<prec[]> coef = 0;
    int nRow = 0;      
    int nCol = 0;   
    int nTerm = 0;  
    //------------------------------------------------
    // CONSTRUCTOR
    //------------------------------------------------
    CSRMAT()
    {
        // int* iat = nullptr;
        // int* ja = nullptr;
        // prec* coef = nullptr;
    }
    //---
    CSRMAT(int nRows, int nCols, int N)
    {
        nRow  = nRows;
        nCol  = nCols;
        nTerm = N;
        if (N == 0) throw_line("ERROR: Defining sparse matrix with zero non null entries.\n");
        iat  = VECTOR_INT::makePointer(nRows+1);
        ja   = VECTOR_INT::makePointer(N);
        coef = VECTOR::makePointer(N);
    }
    //---
    CSRMAT(int nRows, int nCols, int N, std::shared_ptr<int[]> &iatP, std::shared_ptr<int[]> &jaP)
    {
        nRow  = nRows;
        nCol = nCols;
        nTerm = N;
        if (N == 0) throw_line("ERROR: Defining sparse matrix with zero non null entries.\n");
        
        iat  = iatP;
        ja = jaP;
        coef = nullptr; 
        coef = VECTOR::makePointer(N);
    }
    //---
    CSRMAT(int nRows, int nCols, std::shared_ptr<int[]> &iatP, std::shared_ptr<int[]> &jaP, std::shared_ptr<prec[]> &coefP, int N)
    {
        nRow  = nRows;
        nCol = nCols;
        nTerm = N;
        if (N == 0) throw_line("ERROR: Defining sparse matrix with zero non null entries.\n");
        
        iat  = iatP;
        ja   = jaP;
        coef = coefP;
    }
    
    //------------------------------------------------------------------------------
    // CONTRUCTOR FOR THE STIFFNESS MATRIX ASSEMBLY
    //------------------------------------------------------------------------------
    CSRMAT(int nRows, int nCols, int N, std::shared_ptr<int[]> &rowID, std::shared_ptr<int[]> &colID, std::shared_ptr<prec[]> &coefP) 
    {
        initialize(nRows, nCols, N, rowID, colID, coefP);
    }
    //---
    CSRMAT(int &nRows, int &nCols, int &N,  VECTOR_INT &rowIDV, VECTOR_INT &colIDV, VECTOR &coefV) 
    {
        std::shared_ptr<int[]> rowID = rowIDV.P;
        std::shared_ptr<int[]> colID = colIDV.P;
        std::shared_ptr<prec[]> coef  = coefV.P;
        CSRMAT(nRows, nCols, N, rowID, colID, coef);
    }
    //-----------------------------------------------------------
    // INITIALIZE
    //-----------------------------------------------------------
    void initialize(int nRows, int nCols, int N)
    {
        nRow  = nRows;
        nCol  = nCols;
        nTerm = N;
        if (N == 0) throw_line("ERROR: Defining sparse matrix with zero non null entries.\n");
        iat  = VECTOR_INT::makePointer(nRows+1);
        ja   = VECTOR_INT::makePointer(N);
        coef = VECTOR::makePointer(N);
    }
    //---------------------------------------------------------------------------------
    void initialize(int nRows, int nCols, int N, std::shared_ptr<int[]> &iatP, std::shared_ptr<int[]> &jaP)
    {
        nRow  = nRows;
        nCol = nCols;
        nTerm = N;
        if (N == 0) throw_line("ERROR: Defining sparse matrix with zero non null entries.\n");
        
        iat  = iatP;
        ja = jaP;
        coef = nullptr; 
        coef = VECTOR::makePointer(N);
    }
    //---
    void initialize(int &nRows, int &nCols, std::shared_ptr<int[]> &iatP, std::shared_ptr<int[]> &jaP, std::shared_ptr<prec[]> &coefP, int &N)
    {
        nRow  = nRows;
        nCol = nCols;
        nTerm = N;
        if (N == 0) throw_line("ERROR: Defining sparse matrix with zero non null entries.\n");
        
        iat  = iatP;
        ja   = jaP;
        coef = coefP;
    }
    //---
    void initialize(int &nRows, int &nCols, int &N,  std::shared_ptr<int[]> &rowID, std::shared_ptr<int[]> &colID, std::shared_ptr<prec[]> &coefP) 
    {
        // inizialize parameters
        nRow  = nRows;
        nCol  = nCols;
        if (N == 0) throw_line("ERROR: Defining sparse matrix with zero non null entries.\n");

        // count the number of columns in each row
        VECTOR_INT nColInRow = VECTOR_INT::zeros(nRows);
        for (int i = 0; i < N; i++)
        {
            nColInRow[rowID[i]]++;
        }

        // allocate pointer of pointers for the columns & coeffs of each row
        std::shared_ptr<int[]> rowBuffCol = VECTOR_INT::makePointer(N);
        std::shared_ptr<int*[]> rowPPCol  =  VECTOR_INT::makeDoublePointer(nRows);      
        std::shared_ptr<prec[]> rowBuffCoeff = VECTOR::makePointer(N);
        std::shared_ptr<prec*[]> rowPPCoeff  = VECTOR::makeDoublePointer(nRows); 
        
        // fulfill *pointers with the **pointers
        int count = 0;
        for (int i = 0; i < nRows; i++)
        {
            rowPPCol[i]   = &rowBuffCol[count];
            rowPPCoeff[i] = &rowBuffCoeff[count];
            count += nColInRow[i];
        }
         
        // fulfil the *pointers
        // free(nColInRow);
        //int caso;
        nColInRow.setZeros(nRows);

        for (int i = 0; i < N; i++)
        {
            rowPPCol[rowID[i]][nColInRow[rowID[i]]]   = colID[i] ;           
            rowPPCoeff[rowID[i]][nColInRow[rowID[i]]] = coefP[i] ;
            nColInRow[rowID[i]]++;
        }
        
        // sort each *pointer of the **pointer for columns and copy changes to **pointer for coeffs, taking all the values once
        count = 0;
        std::shared_ptr<int[]> repeatedColInRow = VECTOR_INT::makePointer(nRows);
        for (int i = 0; i < nRows; i++) 
        {
            repeatedColInRow[i] = sort(rowPPCol[i], rowPPCoeff[i], nColInRow[i]);
            count += repeatedColInRow[i];
        }

        nTerm = N - count;

        //CREATING JA % IAT
        std::shared_ptr<int[]> irow = VECTOR_INT::makePointer(nTerm);
        iat       = VECTOR_INT::makePointer(nRow+1);
        ja        = VECTOR_INT::makePointer(nTerm);
        coef      = VECTOR::makePointer(nTerm);

        int idRow = 0;
        int countCol = -1;
        int incumbentNnzInRow = nColInRow[0]-repeatedColInRow[0];
        for (int i = 0; i < nTerm; i++)
        {
            //*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
            while (incumbentNnzInRow == 0)
            {
                idRow++;
                if (idRow == nRow) break;
                incumbentNnzInRow = nColInRow[idRow]-repeatedColInRow[idRow];
            }
            countCol++;
            irow[i] = idRow;
            ja[i]   = rowPPCol[idRow][countCol];
            coef[i] = rowPPCoeff[idRow][countCol];
            if (countCol == incumbentNnzInRow - 1)
            {
                idRow++;
                if (idRow == nRow) break;
                incumbentNnzInRow = nColInRow[idRow]-repeatedColInRow[idRow];
                countCol = -1;
            } 
        }
        irow2iat(nRow, nTerm, irow, iat);
    }
    //-----------------------------------
    // REFERENCED INITIALIZE
    //-----------------------------------
    void refInitialize(CSRMAT &matIn)
    {
        nRow = matIn.nRow;
        nCol = matIn.nCol;
        nTerm = matIn.nTerm;
        iat = matIn.iat;
        ja = matIn.ja;
        coef = matIn.coef;
    }
    //-----------------------------------
    // ZEROS
    //-----------------------------------
     void defineZero(int nRows, int nCols)
    {
        nRow = nRows;
        nCol = nCols;
        prec* coefP = (prec*) calloc(nTerm, sizeof(prec));
        coef = std::shared_ptr<prec[]>(new(coefP) prec[nTerm], [](prec* vecPointer){free(vecPointer);});
    }
    //-----
    void defineZero(int nRows, int nCols, int N, std::shared_ptr<int[]> &jaIn, std::shared_ptr<int[]> &iatIn)
    {
        nRow = nRows;
        nCol = nCols;
        nTerm = N;
        ja = jaIn;
        iat = iatIn;
        prec* coefP = (prec*) calloc(nTerm, sizeof(prec));
        coef = std::shared_ptr<prec[]>(new(coefP) prec[nTerm], [](prec* vecPointer){free(vecPointer);});
    }
    //-----------------------------------
    // CREATE IAT FROM ROW IDS
    //-----------------------------------
    static void irow2iat(int nRows, int N, std::shared_ptr<int[]> &irow, std::shared_ptr<int[]> &iat)
    {
        int irow_old = -1;
        int irow_new;
        for (int k = 0; k < N; k++)
        {
            irow_new = irow[k];
            if (irow_new > irow_old)
            {
                for (int j = irow_old + 1; j < irow_new + 1; j++) iat[j] = k;
                irow_old = irow_new;
            }
        }
        int k = N;
        for (int j = irow_old + 1; j < nRows + 1; j++) iat[j] = k;
    }
    //--------------------------------------------------------------
    // ACCESS TO ALREADY !=0 ENTRY WITH STANDARD COORDINATES
    //--------------------------------------------------------------
    int getPos (int iglob, int jglob)
    {
        int startRow = iat[iglob];
        int endRow   = iat[iglob+1];
        for (int pos = startRow; pos < endRow; pos++)
        {
            if (ja[pos] == jglob)
            {
                return pos;
            } 
        }
        throw_line("WARNING: trying to access a non allocated position of SPARSE matrix\n");
    }
    //---
    prec getEl(int iglob, int jglob)
    {
        int startRow = iat[iglob];
        int endRow   = iat[iglob+1];
        prec value;
        for (int pos = startRow; pos < endRow; pos++)
        {
            if (ja[pos] == jglob)
            {
                value = coef[pos];
                return value;
            } 
        }
        return 0;
        // throw_line("WARNING: trying to access a non allocated position of SPARSE matrix\n");
    }
    //---
    int operator () (int iglob, int jglob, prec value)
    {
        int startRow = iat[iglob];
        int endRow   = iat[iglob+1];
        for (int pos = startRow; pos < endRow; pos++)
        {
            if (ja[pos] == jglob)
            {
                coef[pos] = value;
                return pos;
            } 
        }
        throw_line("WARNING: trying to access a non allocated position of SPARSE matrix\n");
    }
    //---
    prec& operator () (int iglob, int jglob)
    {
        int startRow = iat[iglob];
        int endRow   = iat[iglob+1];
        for (int pos = startRow; pos < endRow; pos++)
        {
            if (ja[pos] == jglob) return coef[pos];
        }
        throw_line("WARNING: trying to access a non allocated position of SPARSE matrix\n");
    }
    
    //--------------------------------------------
    // SET ROW TO ZERO
    //--------------------------------------------
    void setRowToZero(VECTOR_INT &rows)
    {
        for (int irow = 0; irow < rows.length; irow++) setRowToZero(rows[irow]);
    }
    //---
    void setRowToZero(int iRow)
    {
        int startPos = iat[iRow];
        int endPos = iat[iRow+1];
        for (int pos = startPos; pos < endPos; pos++)
        {
            coef[pos] = 0;
        }
    }
    //---
    VECTOR_INT setColInRowToZero(int iRow, VECTOR_INT &posToZero, VECTOR_INT &iglobToZero, int &countPosToZero, VECTOR_INT &posToOne, int &countPosToOne)
    {
        int count = 0;
        int startPos = iat[iRow];
        int endPos = iat[iRow+1];
        VECTOR_INT colInRowToZero(endPos-startPos - 1);
        // printf("startPos: %d, end pos: %d\n", startPos, endPos);
        for (int pos = startPos; pos < endPos; pos++)
        {
            int col = ja[pos];

            if (col == iRow) 
            {
                coef[pos] = 1.0; 
                posToOne[countPosToOne] = pos; countPosToOne++;
                continue;
            }
            colInRowToZero[count] = col;
            coef[pos] = 0.0;
            posToZero[countPosToZero] = pos; 
            iglobToZero[countPosToZero] = iRow;
            countPosToZero++;
            count++;
        }
        return colInRowToZero;
    }
    //---
    VECTOR_INT setColInRowToZero(int iRow, VECTOR_INT &posToZero, int &countPosToZero, VECTOR_INT &posToOne, int &countPosToOne)
    {
        int count = 0;
        int startPos = iat[iRow];
        int endPos = iat[iRow+1];
        VECTOR_INT colInRowToZero(endPos-startPos - 1);
        // printf("startPos: %d, end pos: %d\n", startPos, endPos);
        for (int pos = startPos; pos < endPos; pos++)
        {
            int col = ja[pos];

            if (col == iRow) 
            {
                coef[pos] = 1.0; 
                posToOne[countPosToOne] = pos; countPosToOne++;
                continue;
            }
            colInRowToZero[count] = col;
            coef[pos] = 0.0;
            posToZero[countPosToZero] = pos; countPosToZero++;
            count++;
        }
        return colInRowToZero;
    }
    //---
    VECTOR_INT setColInRowToZero(int iRow)
    {
        int count = 0;
        int startPos = iat[iRow];
        int endPos = iat[iRow+1];
        VECTOR_INT colInRowToZero(endPos-startPos - 1);
        // printf("startPos: %d, end pos: %d\n", startPos, endPos);
        for (int pos = startPos; pos < endPos; pos++)
        {
            int col = ja[pos];

            if (col == iRow) 
            {
                coef[pos] = 1.0; 
                continue;
            }
            colInRowToZero[count] = col;
            coef[pos] = 0.0;
            count++;
        }
        return colInRowToZero;
    }
    //---
    //--------------------------------------------------------------
    // OPERATIONS OVERLOAD
    //--------------------------------------------------------------
    // SUM
    //-----------------------------------
    // ONLY FOR IDENTICAL PATTERNS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    void sumSamePatt(CSRMAT &matIn)
    {
        if (nRow != matIn.nRow || nCol != matIn.nCol) throw_line("ERROR: matrix sizes are not the same.\n");
        
        std::shared_ptr<prec[]> coefIn = matIn.coef;
        for (int pos = 0; pos < nTerm; pos++)
        {
            coef[pos] += coefIn[pos];
        }
    }
    //---
    CSRMAT operator + (CSRMAT &matIn)
    {
        if (nRow != matIn.nRow || nCol != matIn.nCol) throw_line("ERROR: matrix sizes are not the same.\n");
        std::shared_ptr<int[]> jaIn  = matIn.ja;
        std::shared_ptr<int[]> iatIn = matIn.iat; 
        std::shared_ptr<prec[]> coefIn = matIn.coef;

        // std::shared_ptr<int[]> jaRes  = VECTOR_INT::makePointer(nTerm+matIn.nTerm);
        std::shared_ptr<int[]> iatRes = VECTOR_INT::makePointer(nRow+1);
        // std::shared_ptr<prec[]> coefRes = VECTOR::makePointer(nTerm+matIn.nTerm);

        VECTOR_INT jaRes(nTerm+matIn.nTerm);
        VECTOR coefRes(nTerm+matIn.nTerm);

        int count = 0;
        for (int i = 0; i < nRow; i++)
        {
            int nColInRow    = iat[i+1]   - iat[i];
            int nColInRow_in = iatIn[i+1] - iatIn[i];
            int nColTot = nColInRow + nColInRow_in;
            // merge values to appendCol
            std::shared_ptr<int[]> appendCol = VECTOR_INT::makePointer(nColInRow+nColInRow_in);
            for (int j = 0; j < nColInRow; j++) appendCol[j] = ja[iat[i]+j];
            for (int j = 0; j < nColInRow_in; j++) appendCol[nColInRow+j] = jaIn[iatIn[i]+j];

            //----
            // merge values to appendCoef
            std::shared_ptr<prec[]> appendCoef = VECTOR::makePointer(nColInRow+nColInRow_in);
            for (int j = 0; j < nColInRow; j++) appendCoef[j] = coef[iat[i]+j];
            for (int j = 0; j < nColInRow_in; j++) appendCoef[nColInRow+j] = coefIn[iatIn[i]+j];

            int repetedCol = sort(appendCol, appendCoef, nColTot);
            
            int nNZ = nColTot-repetedCol;
            iatRes[i] = count;
        

            VECTOR_INT::copy(appendCol, jaRes.P, nNZ, count);
            VECTOR::copy(appendCoef, coefRes.P, nNZ, count);

            count += nNZ; 
        }
        iatRes[nRow] = count;

        jaRes.shrink(count); coefRes.shrink(count);

        CSRMAT matOut(nRow, nCol, iatRes, jaRes.P, coefRes.P, count);
        return matOut;
    }
    //------
    // DIFF
    //------
    CSRMAT operator - (CSRMAT &matIn)
    {
        CSRMAT tempMat = matIn*(-1);
        return *this + tempMat;
    }
    
    void minus()
    {
        for (int pos = 0; pos < nTerm; pos++) coef[pos] = -coef[pos];
    }
    //-----------------------------------
    // MATRIX BY VECTOR PRODUCT
    //-----------------------------------
    VECTOR operator * (VECTOR &vec)
    {
        if (vec.length != nCol) 
        {
            throw_line("\nERROR, wrong matrix by product dimensions\n");

        }
        std::shared_ptr<prec[]> vecP = vec.P;
        // local variables
        int endRow,beginRow;
        VECTOR res(nRow);

        //--------------------

        for (int k = 0; k < nRow; k++)
        {
            endRow   = iat[k];
            beginRow = iat[k+1];
            res[k]   = coef[endRow] * vec[ja[endRow]];
            for (int i = endRow+1; i < beginRow; i++)
            {
                res[k] = res[k] + coef[i]*vecP[ja[i]];
            }
        }
        return res;
    }
    //---
    VECTOR prod(VECTOR &vec)
    {
        if (vec.length != nCol) 
        {
            throw_line("\nERROR, wrong matrix by product dimensions\n");

        }
        std::shared_ptr<prec[]> vecP = vec.P;
        // local variables
        int endRow,beginRow;
        VECTOR res(nRow);

        //--------------------

        for (int k = 0; k < nRow; k++)
        {
            endRow   = iat[k];
            beginRow = iat[k+1];
            res[k]   = coef[endRow] * vec[ja[endRow]];
            for (int i = endRow+1; i < beginRow; i++)
            {
                res[k] = res[k] + coef[i]*vecP[ja[i]];
            }
        }
        return res;
    }
    //---
    void prod(VECTOR &vec, VECTOR& res)
    {
        if (vec.length != nCol) 
        {
            throw_line("\nERROR, wrong matrix by product dimensions\n");

        }
        std::shared_ptr<prec[]> vecP = vec.P;
        // local variables
        int endRow,beginRow;

        //--------------------

        for (int k = 0; k < nRow; k++)
        {
            endRow   = iat[k];
            beginRow = iat[k+1];
            res[k]   = coef[endRow] * vec[ja[endRow]];
            for (int i = endRow+1; i < beginRow; i++)
            {
                res[k] = res[k] + coef[i]*vecP[ja[i]];
            }
        }
    }
    
    //---
    prec* prod(prec* &vec, int N)
    {

        if (N != nCol) 
        {
            throw_line("\nERROR, wrong matrix by product dimensions\n");

        }
        // local variables
        int endRow,beginRow;
        prec* res = (prec*) malloc(nRow * sizeof(prec));

        //--------------------

        for (int k = 0; k < nRow; k++)
        {
            endRow   = iat[k];
            beginRow = iat[k+1];
            res[k]   = coef[endRow] * vec[ja[endRow]];
            for (int i = endRow+1; i < beginRow; i++)
            {
                res[k] = res[k] + coef[i]*vec[ja[i]];
            }
        }
        return res;
    }
    //---
    std::shared_ptr<prec[]> prod(std::shared_ptr<prec[]> &vec, int N)
    {

        if (N != nCol) 
        {
            throw_line("\nERROR, wrong matrix by product dimensions\n");

        }
        // local variables
        int endRow,beginRow;
        std::shared_ptr<prec[]> res(new prec[nRow]);

        //--------------------

        for (int k = 0; k < nRow; k++)
        {
            endRow   = iat[k];
            beginRow = iat[k+1];
            res[k]   = coef[endRow] * vec[ja[endRow]];
            for (int i = endRow+1; i < beginRow; i++)
            {
                res[k] = res[k] + coef[i]*vec[ja[i]];
            }
        }
        return res;
    }
    //---
    void prod(prec* vec, int N, prec* &res)
    {

        if (N != nCol) 
        {
            throw_line("\nERROR, wrong matrix by product dimensions\n");

        }
        // local variables
        int endRow,beginRow;
        if(res == nullptr) res = (prec*) malloc(nRow * sizeof(prec));

        //--------------------

        for (int k = 0; k < nRow; k++)
        {
            endRow   = iat[k];
            beginRow = iat[k+1];
            res[k]   = coef[endRow] * vec[ja[endRow]];
            for (int i = endRow+1; i < beginRow; i++)
            {
                res[k] = res[k] + coef[i]*vec[ja[i]];
            }
        }
    }

    //-------------------------------------
    // MULTIPLY BY SCALAR
    //-------------------------------------
    CSRMAT operator * (prec val)
    {
        CSRMAT resCSR(nRow, nCol, nTerm);
        copyPattTo(resCSR);
        std::shared_ptr<prec[]> res = resCSR.coef;
        VECTOR::multiplyByScalar(coef, val, nTerm, res);
        return resCSR;
    }
    //---
    CSRMAT operator * (int val)
    {
        CSRMAT resCSR(nRow, nCol, nTerm);
        copyPattTo(resCSR);
        std::shared_ptr<prec[]> res = resCSR.coef;
        VECTOR::multiplyByScalar(coef, val, nTerm, res);
        return resCSR;
    }

    // MATRIX BY MATRIX MULTIPLACTION
    CSRMAT operator * (CSRMAT &mat) 
    {
        if (nCol != mat.nRow) throw_line("ERROR: SPARSE MATRIX PRODUCT, inconstistend dimensions.\n")
        // int nRowIn = mat.nRow;
        int nColIn = mat.nCol;
        std::shared_ptr<int[]> iatIn = mat.iat;
        std::shared_ptr<int[]> jaIn = mat.ja; 
        int nTermIn = mat.nTerm;
        std::shared_ptr<prec[]> coefIn = mat.coef;
        //----------------------
        int incumbMaxStor = 2*nTerm;
        std::shared_ptr<int[]> iatRes(new int[nRow+1]);
        VECTOR_INT jaRes(incumbMaxStor);
        VECTOR coefRes(incumbMaxStor);

        int pos = 0;
        int posRes = 0;
        iatRes[0] = 0;
        for (int irow = 0; irow < nRow; irow++)
        {
            int nElInRow = iat[irow+1] - iat[irow];
            VECTOR_INT tempColInRow(nElInRow);
            VECTOR tempCoefInRow(nElInRow);
            for (int el = 0; el < nElInRow; el++)
            {
                tempColInRow[el]  = ja[pos];
                tempCoefInRow[el] = coef[pos];
                pos++;
            }

            VECTOR_INT colInResRow(nElInRow*nTermIn);
            VECTOR coefInResRow(nElInRow*nTermIn);
            int tempCount = 0;
            //------------------------------------
            for (int iter = 0; iter < nElInRow; iter++)
            {
                int jrow = tempColInRow[iter];
                for (int posIn = iatIn[jrow]; posIn < iatIn[jrow+1]; posIn++)
                {
                    int tempColIn = jaIn[posIn];
                    colInResRow[tempCount]  = tempColIn;
                    coefInResRow[tempCount] = tempCoefInRow[iter] * coefIn[posIn];
                    tempCount++;
                }
            }
            int repitedCol = sort(colInResRow.P, coefInResRow.P, tempCount);
            int tempNnz = tempCount - repitedCol;
            for (int tempPos = 0; tempPos < tempNnz; tempPos++)
            {
                jaRes[posRes]   = colInResRow[tempPos];
                coefRes[posRes] = coefInResRow[tempPos];
                posRes++;
                if (posRes == incumbMaxStor)
                {
                    prec perc = (prec) irow/nRow;
                    incumbMaxStor += incumbMaxStor*(1-perc)/perc;
                    jaRes.enlarge(incumbMaxStor);
                    coefRes.enlarge(incumbMaxStor);
                }
            }
            iatRes[irow+1] = posRes;
            // free
            //-----
        }
        jaRes.shrink(posRes); coefRes.shrink(posRes);
        CSRMAT resMat(nRow, nColIn, iatRes, jaRes.P, coefRes.P, posRes);
        return resMat;
    }
    //-------------------------------------------------------------------
    static void prod(CSRMAT &mat1, CSRMAT &mat2, CSRMAT &resMat) 
    {
        int nCol1 = mat1.nCol; 
        int nRow1 = mat1.nRow;
        int nTerm1 = mat1.nTerm;
        int nCol2 = mat2.nCol;
        int nRow2 = mat2.nRow;
        if (nCol1 != nRow2)
            throw_line("ERROR: SPARSE MATRIX PRODUCT, inconstistend dimensions.\n")
                // int nRowIn = mat.nRow;

        std::shared_ptr<int[]> iat1 = mat1.iat;
        std::shared_ptr<int[]> ja1 = mat1.ja;
        std::shared_ptr<prec[]> coef1 = mat1.coef;

        std::shared_ptr<int[]> iat2 = mat2.iat;
        std::shared_ptr<int[]> ja2 = mat2.ja;
        std::shared_ptr<prec[]> coef2 = mat2.coef;
        //----------------------
        int incumbMaxStor = 2*nTerm1;
        std::shared_ptr<int[]> iatRes(new int[nRow1+1]);
        VECTOR_INT jaRes(incumbMaxStor);
        VECTOR coefRes(incumbMaxStor);

        int pos = 0;
        int posRes = 0;
        iatRes[0] = 0;
        for (int irow = 0; irow < nRow1; irow++)
        {
            int nElInRow = iat1[irow+1] - iat1[irow];
            VECTOR_INT tempColInRow(nElInRow);
            VECTOR tempCoefInRow(nElInRow);
            for (int el = 0; el < nElInRow; el++)
            {
                tempColInRow[el]  = ja1[pos];
                tempCoefInRow[el] = coef1[pos];
                pos++;
            }

            VECTOR_INT colInResRow(nElInRow*nCol2);
            VECTOR coefInResRow(nElInRow*nCol2);
            int tempCount = 0;
            //------------------------------------
            for (int iter = 0; iter < nElInRow; iter++)
            {
                int jrow = tempColInRow[iter];
                for (int posIn = iat2[jrow]; posIn < iat2[jrow+1]; posIn++)
                {
                    int tempCol2 = ja2[posIn];
                    colInResRow[tempCount]  = tempCol2;
                    coefInResRow[tempCount] = tempCoefInRow[iter] * coef2[posIn];
                    tempCount++;
                }
            }
            int repitedCol = sort(colInResRow.P, coefInResRow.P, tempCount);
            int tempNnz = tempCount - repitedCol;
            for (int tempPos = 0; tempPos < tempNnz; tempPos++)
            {
                jaRes[posRes]   = colInResRow[tempPos];
                coefRes[posRes] = coefInResRow[tempPos];
                posRes++;
                if (posRes == incumbMaxStor)
                {
                    prec perc = (prec) irow/nRow1;
                    incumbMaxStor += incumbMaxStor*(1-perc)/perc;
                    jaRes.enlarge(incumbMaxStor);
                    coefRes.enlarge(incumbMaxStor);
                }
            }
            iatRes[irow+1] = posRes;
            // free
            //-----
        }
        jaRes.shrink(posRes); coefRes.shrink(posRes);
        resMat.initialize(nRow1, nCol2, iatRes, jaRes.P, coefRes.P, posRes);
    }
    //------------------------------
    // CSRMAT operator * (CSRMAT mat) //ONLY IF RESULTING MAT HAS ONLY NONZEROS ON THE MAIN DIAGONAL
    // {
    //     if (nCol != mat.nRow) throw_line("ERROR: SPARSE MATRIX PRODUCT, inconstistend dimensions.\n");
    //     CSRMAT tempMat = mat.getTranspose();
    //     int nRowIn = tempMat.nRow;
    //     std::shared_ptr<int[]> iatIn = tempMat.iat;
    //     std::shared_ptr<int[]> jaIn = tempMat.ja; 
    //     int nTermIn = mat.nTerm;
    //     std::shared_ptr<prec[]> coefIn = tempMat.coef;
    //     //----------------------
    //     int rescale =  2*nTerm;
    //     int incumbMaxStor = 2*nTerm;
    //     std::shared_ptr<int[]> iatRes(new int[nRow+1]);
    //     VECTOR_INT jaRes(incumbMaxStor);
    //     VECTOR coefRes(incumbMaxStor);

    //     int pos = 0;
    //     int posRes = 0;
    //     iatRes[0] = 0;
    //     for (int irow = 0; irow < nRow; irow++)
    //     {
    //         int nElInRow = iat[irow+1] - iat[irow];
    //         VECTOR_INT tempColInRow(nElInRow);
    //         VECTOR tempCoefInRow(nElInRow);
    //         for (int el = 0; el < nElInRow; el++)
    //         {
    //             tempColInRow[el]  = ja[pos];
    //             tempCoefInRow[el] = coef[pos];
    //             pos++;
    //         }
    //         //------------------------------------
    //         for (int jrow = 0; jrow < nRowIn; jrow++)
    //         {
    //             prec tempVal = 0;
    //             bool found = false;
    //             for (int posIn = iatIn[jrow]; posIn < iatIn[jrow+1]; posIn++)
    //             {
    //                 int tempId;
    //                 int tempColIn = jaIn[posIn];
    //                 if (tempColInRow.hasIn(tempColIn, tempId))
    //                 {
    //                     found = true;
    //                     tempVal += tempCoefInRow[tempId] * coefIn[posIn];
    //                 }
    //             }
    //             if (found)
    //             {
    //                 jaRes[posRes] = jrow;
    //                 coefRes[posRes] = tempVal;
    //                 posRes++;
    //             }
    //             found = false;
    //         }
    //         iatRes[irow+1] = posRes;
    //         // free
    //         //-----
    //     }
    //     jaRes.shrink(posRes); coefRes.shrink(posRes);
    //     CSRMAT resMat(nRow, nRowIn, iatRes, jaRes.P, coefRes.P, posRes);
    //     return resMat;
    // }


    //-------------------------------------
    // DIVIDE BY SCALAR
    //-------------------------------------
    void divideByScalar(prec val, CSRMAT &resCSR)
    {
        std::shared_ptr<prec[]> res = resCSR.coef;
        VECTOR::divideByScalar(coef, val, nTerm, res);
    }
    //---
    CSRMAT operator / (prec val)
    {
        CSRMAT resCSR(nRow, nCol, nTerm);
        copyPattTo(resCSR);
        std::shared_ptr<prec[]> res = resCSR.coef;
        VECTOR::divideByScalar(coef, val, nTerm, res);
        return resCSR;
    }
    //---
    CSRMAT operator / (int val)
    {
        CSRMAT resCSR(nRow, nCol, nTerm);
        copyPattTo(resCSR);
        std::shared_ptr<prec[]> res = resCSR.coef;
        VECTOR::divideByScalar(coef, val, nTerm, res);
        return resCSR;
    }

    //-----------------------------------
    // GET DIAGONAL
    //-----------------------------------
    void diag(std::shared_ptr<prec[]> &diagonal)
    {
        if (diagonal == 0) diagonal = VECTOR::makePointer(nRow);
        for (int i = 0; i < nRow; i++)
        {
            for (int j = iat[i]; j < iat[i+1]; j++)
            {
                if (ja[j] == i)
                {
                    diagonal[i] = coef[j];
                    break;
                }
            }
        }
    }
    //---
    std::shared_ptr<prec[]> diag()
    {
        std::shared_ptr<prec[]> diagonal(new prec[nRow]);
        for (int i = 0; i < nRow; i++)
        {
            for (int j = iat[i]; j < iat[i+1]; j++)
            {
                if (ja[j] == i)
                {
                    diagonal[i] = coef[j];
                    break;
                }
            }
        }
        return diagonal;
    }
    //---
    //---
    void precondNSDiag(std::shared_ptr<prec[]> &diagonal)
    {
        if (diagonal == 0) diagonal = VECTOR::makePointer(nRow);
        VECTOR_INT vec;
        vec.setZeros(nRow);
        for (int i = 0; i < nRow; i++)
        {
            for (int j = iat[i]; j < iat[i+1]; j++)
            {
                if (ja[j] == i)
                {
                    vec[i] = 1;
                    diagonal[i] = coef[j];
                    break;
                }
            }
            if (vec[i] == 0) diagonal[i] = 1;
        }
    }
    //---
    std::shared_ptr<prec[]> precondNSDiag()
    {
        std::shared_ptr<prec[]> diagonal(new prec[nRow]);
        VECTOR vec;
        vec.setZeros(nRow);
        for (int i = 0; i < nRow; i++)
        {
            for (int j = iat[i]; j < iat[i+1]; j++)
            {
                if (ja[j] == i)
                {
                    vec[i] = 1;
                    diagonal[i] = coef[j];
                    break;
                }
            }
            if (abs(vec[i] - 1) > 1e-15) diagonal[i] = 1;
        }
        return diagonal;
    }
    //---
    CSRMAT RprodByDiag(VECTOR &vec)
    {
        if (nCol != vec.length) throw_line("ERROR: wrong dimensions in matrix by diag\n.");
        CSRMAT outMat(nRow, nCol, nTerm);;
        copyTo(outMat);
        std::shared_ptr<prec[]> coefOut = outMat.coef;
        for (int irow = 0; irow < nRow; irow++)
        {
            for (int pos = iat[irow]; pos < iat[irow+1]; pos++)
            {
                int tempCol = ja[pos];
                coefOut[pos] *= vec[tempCol];
            }
        }
        return outMat;
    }
    //---
    CSRMAT LprodByDiag(VECTOR &vec)
    {
        if (nRow != vec.length) throw_line("ERROR: wrong dimensions in matrix by diag\n.");
        CSRMAT outMat(nRow, nCol, nTerm);
        copyTo(outMat);
        std::shared_ptr<prec[]> coefOut = outMat.coef;
        for (int irow = 0; irow < nRow; irow++)
        {
            for (int pos = iat[irow]; pos < iat[irow+1]; pos++)
            {
                coefOut[pos] *= vec[irow];
            }
        }
        return outMat;
    }

    //-------------------------------------
    // JACOBI PRECONDITIONED CURRENT MATRIX   mat = mat * J^-1
    //-------------------------------------
    void jacobiPrecond()
    {
        std::shared_ptr<prec[]> diagonal = diag();
        for (int j = 0; j < nTerm; j++)
        {
            coef[j] = coef[j] / diagonal[ja[j]];
        }
    }
    //---
    static void jacobiPrecond(CSRMAT &mat)
    {
        std::shared_ptr<prec[]> coeff = mat.coef;
        std::shared_ptr<int[]> ja = mat.ja;
        std::shared_ptr<prec[]> diagonal = mat.diag();
        for (int j = 0; j < mat.nTerm; j++)
        {
            coeff[j] = coeff[j] / diagonal[ja[j]];
        }
    }
    //---
    void jacobiPrecond(std::shared_ptr<prec[]> &diagonal)
    {
        for (int j = 0; j < nTerm; j++)
        {
            coef[j] = coef[j] / diagonal[ja[j]];
        }
    }
    //---
    void jacobiPrecond(prec* &diagonal)
    {
        checkPtr(diagonal);
        for (int j = 0; j < nTerm; j++)
        {
            coef[j] = coef[j] / diagonal[ja[j]];
        }
    }
    //---
    void jacobiNSPrecond()
    {
        std::shared_ptr<prec[]> diagonal = precondNSDiag();
        for (int j = 0; j < nTerm; j++)
        {
            coef[j] = coef[j] / diagonal[ja[j]];
        }
    }
    //---
    static void jacobiNSPrecond(CSRMAT &mat)
    {
        std::shared_ptr<prec[]> coeff = mat.coef;
        std::shared_ptr<int[]> ja = mat.ja;
        std::shared_ptr<prec[]> diagonal = mat.precondNSDiag();
        for (int j = 0; j < mat.nTerm; j++)
        {
            coeff[j] = coeff[j] / diagonal[ja[j]];
        }
    }
    //---
    void jacobiNSPrecond(std::shared_ptr<prec[]> &diagonal)
    {
        for (int j = 0; j < nTerm; j++)
        {
            coef[j] = coef[j] / diagonal[ja[j]];
        }
    }
    //---------------------------------
    // ENLARGE ROWS
    //---------------------------------
    void enlargeRows(int N)
    {
        std::shared_ptr<int[]> newIat (new int[N+1]);

        for (int irow = 0; irow < nRow+1; irow++) newIat[irow] = iat[irow];

        for (int irow = nRow+1; irow < N+1; irow++) newIat[irow] = nTerm;

        iat.reset();

        iat = newIat;
        nRow = N;
    }

    //---------------------------------
    // COPY
    //---------------------------------
    void copyPattTo(CSRMAT& outMat)
    {
        //if (nRow != outMat.nRow || nCol != outMat.nCol) throw_line("ERROR: Trying to copy CSR matrices with different sizes\n");
        // outMat.iat.reset(); iat = 0; outMat.ja.reset(); ja = 0; outMat.coef.reset(); coef = 0;
        outMat.nRow = nRow; outMat.nCol = nCol; outMat.nTerm = nTerm;
        outMat.iat = VECTOR_INT::makePointer(nRow+1);
        outMat.ja  = VECTOR_INT::makePointer(nTerm);
        outMat.coef = VECTOR::makePointer(nTerm);
        VECTOR_INT::copy(iat, outMat.iat, nRow+1);
        VECTOR_INT::copy(ja, outMat.ja, nTerm);
    }
    //---
    void copyPatt(CSRMAT& inMat)
    {
        if (nRow != inMat.nRow || nCol != inMat.nCol) throw_line("ERROR: Trying to copy CSR matrices with different sizes\n")
        VECTOR::copy(inMat.iat, iat, nRow+1);
        VECTOR::copy(inMat.ja, ja, inMat.nTerm);
        nTerm = inMat.nTerm;
    }
    //---
    void copyTo(CSRMAT& outMat)
    {
        int nRowOut = outMat.nRow;
        int nColOut = outMat.nCol;
        int nTermOut = outMat.nTerm;
        if (nRow != nRowOut || nCol != nColOut || nTerm != nTermOut) throw_line("ERROR: Trying to copy CSR matrices with different sizes or terms \n");
        VECTOR::copy(iat, outMat.iat, nRow+1);
        VECTOR::copy(ja, outMat.ja, nTerm);
        VECTOR::copy(coef, outMat.coef, nTerm);
    }
    //---
    void copy(CSRMAT& inMat)
    {
        if (nRow != inMat.nRow || nCol != inMat.nCol || nTerm != inMat.nTerm) throw_line("ERROR: Trying to copy CSR matrices with different sizes or terms \n")
        VECTOR::copy(inMat.iat, iat, nRow+1);
        VECTOR::copy(inMat.ja, ja, nTerm);
        VECTOR::copy(inMat.coef, coef, nTerm);
    }
    //---
    void operator = (CSRMAT &inMat)
    {
        if (nRow != inMat.nRow || nCol != inMat.nCol || nTerm != inMat.nTerm)
        {
            if (iat != 0)
            {
                iat.reset(); iat = 0;
            } 
            if (ja != 0)
            {
                ja.reset(); ja = 0;
            } 
            if (coef != 0)
            {
                coef.reset(); coef = 0;
            } 
            initialize(inMat.nRow, inMat.nCol, inMat.nTerm);
        }
        VECTOR::copy(inMat.iat, iat, nRow+1);
        VECTOR::copy(inMat.ja, ja, nTerm);
        VECTOR::copy(inMat.coef, coef, nTerm);
    }
    //---
    bool operator == (CSRMAT &inMat)
    {
        bool isEqual = true;
        std::shared_ptr<int[]> inIat = inMat.iat;
        std::shared_ptr<int[]> inJa = inMat.ja;
        std::shared_ptr<prec[]> inCoef = inMat.coef;
        for (int i = 0; i < nTerm; i++)
        {
            if ((coef[i] != inCoef[i]) || (ja[i]!=inJa[i])) 
            {
                return false;
            }
        }
        for (int irow = 0; irow < nRow+1; irow++)
        {
            if (iat[irow] != inIat[irow]) return false;
        }
        return isEqual;
    }
    //---
    
    bool checkPattern(CSRMAT &inMat)
    {
        bool isEqual = true;
        std::shared_ptr<int[]> inIat = inMat.iat;
        std::shared_ptr<int[]> inJa = inMat.ja;
        for (int i = 0; i < nTerm; i++)
        {
            if (ja[i]!=inJa[i]) 
            {
                return false;
            }
        }
        for (int irow = 0; irow < nRow+1; irow++)
        {
            if (iat[irow] != inIat[irow]) return false;
        }
        return isEqual;
    }
    
    //--------------------------------------
    // GET IROW
    //--------------------------------------
    static std::shared_ptr<int[]> getIRow(CSRMAT &mat)
    {
        int nRows = mat.nRow;
        int nTerms = mat.nTerm;
        std::shared_ptr<int[]> iat = mat.iat;

        std::shared_ptr<int[]> iRow(new int[nTerms]);
        for (int i = 0; i < nRows; i++)
        {
            for (int j = iat[i]; j < iat[i+1]; j++)
            {
                iRow[j] = i;
            }
        }
        return iRow;
    }
    //-------------------------------------------------------
    // TRANSPOSE
    //-------------------------------------------------------
    void getTranspose(CSRMAT &TMat)
    {
        VECTOR_INT colCounter = VECTOR_INT::zeros(nCol);

        TMat.initialize(nCol, nRow, nTerm);
        std::shared_ptr<int[]> newIat = TMat.iat;
        std::shared_ptr<int[]> newJa  = TMat.ja;
        std::shared_ptr<prec[]> newCoef = TMat.coef;

        for (int pos = 0; pos < nTerm; pos++)
        {
            int tempCol = ja[pos];
            colCounter[tempCol]++;
        }
        //------------------------
        std::shared_ptr<int*[]> newJaPP    = VECTOR_INT::makeDoublePointer(nCol);
        std::shared_ptr<prec*[]> newCoefPP = VECTOR::makeDoublePointer(nCol);
        int count = 0;
        // build new iat
        newIat[0] = 0;
        for (int icol = 0; icol < nCol; icol++)
        {
            newJaPP[icol] = &newJa[count];
            newCoefPP[icol] = &newCoef[count];
            count += colCounter[icol];
            newIat[icol+1] = count;
        }
        //---
        colCounter.zeros();
        // build new ja and new coef
        for (int irow = 0; irow < nRow; irow++)
        {
            for (int pos = iat[irow]; pos < iat[irow+1]; pos++)
            {
                int tempCol = ja[pos];
                int tempCount = colCounter[tempCol];
                newJaPP[tempCol][tempCount] = irow;
                newCoefPP[tempCol][tempCount] = coef[pos];
                colCounter[tempCol]++;
            }
        }
    }
    //-------------------
    CSRMAT getTranspose()
    {
        VECTOR_INT colCounter = VECTOR_INT::zeros(nCol);

        CSRMAT TMat(nCol, nRow, nTerm);
        std::shared_ptr<int[]> newIat = TMat.iat;
        std::shared_ptr<int[]> newJa  = TMat.ja;
        std::shared_ptr<prec[]> newCoef = TMat.coef;

        for (int pos = 0; pos < nTerm; pos++)
        {
            int tempCol = ja[pos];
            colCounter[tempCol]++;
        }
        //------------------------
        std::shared_ptr<int*[]> newJaPP    = VECTOR_INT::makeDoublePointer(nCol);
        std::shared_ptr<prec*[]> newCoefPP = VECTOR::makeDoublePointer(nCol);
        int count = 0;
        // build new iat
        newIat[0] = 0;
        for (int icol = 0; icol < nCol; icol++)
        {
            newJaPP[icol] = &newJa[count];
            newCoefPP[icol] = &newCoef[count];
            count += colCounter[icol];
            newIat[icol+1] = count;
        }
        //---
        colCounter.zeros();
        // build new ja and new coef
        for (int irow = 0; irow < nRow; irow++)
        {
            for (int pos = iat[irow]; pos < iat[irow+1]; pos++)
            {
                int tempCol = ja[pos];
                int tempCount = colCounter[tempCol];
                newJaPP[tempCol][tempCount] = irow;
                newCoefPP[tempCol][tempCount] = coef[pos];
                colCounter[tempCol]++;
            }
        }
        return TMat;
    }
    //-------------------------------------------------------
    // CONCATENATE
    //-------------------------------------------------------
    static void concatDiagInRow(CSRMAT &mat, VECTOR &diag, CSRMAT &resMat, int flag = -1)
    {
        if (flag == -1)
        {
            if (mat.nRow != diag.length) throw_line("ERROR: Wrong dimensions in matrix concatenation\n");
        }
            
        std::shared_ptr<int[]> iatIn = mat.iat; std::shared_ptr<int[]> jaIn = mat.ja; std::shared_ptr<prec[]> coefIn = mat.coef; int nColIn = mat.nCol;
        int N = diag.length;
        int nRow = mat.nRow;

        int nTermRes = mat.nTerm + N; int nCol = nColIn + N;
        resMat.initialize(nRow, nCol, nTermRes);

        std::shared_ptr<prec[]> coefRes = resMat.coef;
        std::shared_ptr<int[]> iatRes = resMat.iat;
        std::shared_ptr<int[]> jaRes  = resMat.ja;

        int pos = 0;
        iatRes[0] = 0;
        for (int irow = 0; irow < nRow; irow++)
        {
            for (int tempPos = iatIn[irow]; tempPos < iatIn[irow+1]; tempPos++)
            {
                jaRes[pos] = jaIn[tempPos];
                coefRes[pos] = coefIn[tempPos];
                pos++;
            }

            jaRes[pos] = irow + nColIn;
            coefRes[pos] = diag[irow];
            pos++;
            iatRes[irow+1] = pos;
        }
    }
    //---
    static CSRMAT concatDiagInRow(CSRMAT &mat, VECTOR &diag, int flag = -1)
    {
        if (flag == -1)
        {
            if (mat.nRow != diag.length) throw_line("ERROR: Wrong dimensions in matrix concatenation\n");
        }
            
        std::shared_ptr<int[]> iatIn = mat.iat; std::shared_ptr<int[]> jaIn = mat.ja; std::shared_ptr<prec[]> coefIn = mat.coef; int nColIn = mat.nCol;
        int N = diag.length;
        int nRow = mat.nRow;

        int nTermRes = mat.nTerm + N; int nCol = nColIn + N;
        CSRMAT resMat(nRow, nCol, nTermRes);

        std::shared_ptr<prec[]> coefRes = resMat.coef;
        std::shared_ptr<int[]> iatRes = resMat.iat;
        std::shared_ptr<int[]> jaRes  = resMat.ja;

        int pos = 0;
        iatRes[0] = 0;
        for (int irow = 0; irow < nRow; irow++)
        {
            for (int tempPos = iatIn[irow]; tempPos < iatIn[irow+1]; tempPos++)
            {
                jaRes[pos] = jaIn[tempPos];
                coefRes[pos] = coefIn[tempPos];
                pos++;
            }

            jaRes[pos] = irow + nColIn;
            coefRes[pos] = diag[irow];
            pos++;
            iatRes[irow+1] = pos;
        }
        return resMat;
    }
     //---
    static void concatDiagInRow(VECTOR &diag, CSRMAT &mat, CSRMAT &resMat, int flag = -1)
    {
        if (flag == -1)
        {
            if (mat.nRow != diag.length) throw_line("ERROR: Wrong dimensions in matrix concatenation\n");
        }
        std::shared_ptr<int[]> iatIn = mat.iat; std::shared_ptr<int[]> jaIn = mat.ja; std::shared_ptr<prec[]> coefIn = mat.coef; int nColIn = mat.nCol;
        int N = diag.length;

        int nRow = mat.nRow;
        int nTermRes = mat.nTerm + N; int nCol = nColIn + N;
        resMat.initialize(nRow, nCol, nTermRes);

        std::shared_ptr<prec[]> coefRes = resMat.coef;
        std::shared_ptr<int[]> iatRes = resMat.iat;
        std::shared_ptr<int[]> jaRes  = resMat.ja;

        int pos = 0;
        iatRes[0] = 0;
        for (int irow = 0; irow < nRow; irow++)
        {
            jaRes[pos] = irow;
            coefRes[pos] = diag[irow];
            pos++;
            for (int tempPos = iatIn[irow]; tempPos < iatIn[irow+1]; tempPos++)
            {
                jaRes[pos] = jaIn[tempPos] + N;
                coefRes[pos] = coefIn[tempPos];
                pos++;
            }
            iatRes[irow+1] = pos;
        }
    }
    //---
    static CSRMAT concatDiagInRow(VECTOR &diag, CSRMAT &mat, int flag = -1)
    {
        if (flag == -1)
        {
            if (mat.nRow != diag.length) throw_line("ERROR: Wrong dimensions in matrix concatenation\n");
        }
        std::shared_ptr<int[]> iatIn = mat.iat; std::shared_ptr<int[]> jaIn = mat.ja; std::shared_ptr<prec[]> coefIn = mat.coef; int nColIn = mat.nCol;
        int N = diag.length;

        int nRow = mat.nRow;
        int nTermRes = mat.nTerm + N; int nCol = nColIn + N;
        CSRMAT resMat(nRow, nCol, nTermRes);

        std::shared_ptr<prec[]> coefRes = resMat.coef;
        std::shared_ptr<int[]> iatRes = resMat.iat;
        std::shared_ptr<int[]> jaRes  = resMat.ja;

        int pos = 0;
        iatRes[0] = 0;
        for (int irow = 0; irow < nRow; irow++)
        {
            jaRes[pos] = irow;
            coefRes[pos] = diag[irow];
            pos++;
            for (int tempPos = iatIn[irow]; tempPos < iatIn[irow+1]; tempPos++)
            {
                jaRes[pos] = jaIn[tempPos] + N;
                coefRes[pos] = coefIn[tempPos];
                pos++;
            }
            iatRes[irow+1] = pos;
        }
        return resMat;
    }
    //----------------
    static CSRMAT concatDiagInCol(CSRMAT &mat, VECTOR &diag, int flag = -1)
    {
        if (flag == 1)
        {
            if (mat.nCol != diag.length) throw_line("ERROR: Wrong dimensions in matrix concatenation\n");
        }
        else
        {
            if (flag <= 0) flag = 0;
        }
        std::shared_ptr<int[]> iatIn = mat.iat; std::shared_ptr<int[]> jaIn = mat.ja; std::shared_ptr<prec[]> coefIn = mat.coef; int nRowIn = mat.nRow;
        int N = diag.length;

        int nCol = mat.nCol;
        int nTermRes = mat.nTerm + N; int nRow = nRowIn + N;
        CSRMAT resMat(nRow, nCol, nTermRes);

        std::shared_ptr<prec[]>coefRes = resMat.coef;
        std::shared_ptr<int[]> iatRes = resMat.iat;
        std::shared_ptr<int[]> jaRes  = resMat.ja;

        int pos = 0;
        iatRes[0] = 0;
        for (int irow = 0; irow < nRowIn; irow++)
        {
            for (int tempPos = iatIn[irow]; tempPos < iatIn[irow+1]; tempPos++)
            {
                jaRes[pos] = jaIn[tempPos];
                coefRes[pos] = coefIn[tempPos];
                pos++;
            }
            iatRes[irow+1] = pos;
        }
        //
        for (int irow = 0; irow < N; irow++)
        {
            jaRes[pos] = irow + flag;
            coefRes[pos] = diag[irow];
            pos++;
            iatRes[irow+1+nRowIn] = pos;
        }
        return resMat;
    }
    //---
    static void concatDiagInCol(CSRMAT &mat, VECTOR &diag, CSRMAT &resMat, int flag = -1)
    {
        if (flag == 1)
        {
            if (mat.nCol != diag.length) throw_line("ERROR: Wrong dimensions in matrix concatenation\n");
        }
        else
        {
            if (flag <= 0) flag = 0;
        }
        std::shared_ptr<int[]> iatIn = mat.iat; std::shared_ptr<int[]> jaIn = mat.ja; std::shared_ptr<prec[]> coefIn = mat.coef; int nRowIn = mat.nRow;
        int N = diag.length;

        int nCol = mat.nCol;
        int nTermRes = mat.nTerm + N; int nRow = nRowIn + N;
        resMat.initialize(nRow, nCol, nTermRes);

        std::shared_ptr<prec[]>coefRes = resMat.coef;
        std::shared_ptr<int[]> iatRes = resMat.iat;
        std::shared_ptr<int[]> jaRes  = resMat.ja;

        int pos = 0;
        iatRes[0] = 0;
        for (int irow = 0; irow < nRowIn; irow++)
        {
            for (int tempPos = iatIn[irow]; tempPos < iatIn[irow+1]; tempPos++)
            {
                jaRes[pos] = jaIn[tempPos];
                coefRes[pos] = coefIn[tempPos];
                pos++;
            }
            iatRes[irow+1] = pos;
        }
        //
        for (int irow = 0; irow < N; irow++)
        {
            jaRes[pos] = irow + flag;
            coefRes[pos] = diag[irow];
            pos++;
            iatRes[irow+1+nRowIn] = pos;
        }
    }
    //---
    //---
    static void concatDiagInCol(VECTOR &diag, CSRMAT &mat, CSRMAT &resMat, int flag = -1)
    {
        if (flag == 1)
        {
            if (mat.nCol != diag.length) throw_line("ERROR: Wrong dimensions in matrix concatenation\n");
        }
        else 
        {
            if (flag <= 0) flag = 0;
        }
        std::shared_ptr<int[]> iatIn = mat.iat; std::shared_ptr<int[]> jaIn = mat.ja; std::shared_ptr<prec[]> coefIn = mat.coef; int nRowIn = mat.nRow;
        int N = diag.length;

        int nCol = mat.nCol;
        int nTermRes = mat.nTerm + N; int nRow = nRowIn + N;
        resMat.initialize(nRow, nCol, nTermRes);

        std::shared_ptr<prec[]> coefRes = resMat.coef;
        std::shared_ptr<int[]> iatRes = resMat.iat;
        std::shared_ptr<int[]> jaRes  = resMat.ja;

        int pos = 0;
        iatRes[0] = 0;
        for (int irow = 0; irow < N; irow++)
        {
            jaRes[pos] = irow + flag;
            coefRes[pos] = diag[irow];
            pos++;
            iatRes[irow+1] = pos;
        }
        //
        for (int irow = 0; irow < nRowIn; irow++)
        {
            for (int tempPos = iatIn[irow]; tempPos < iatIn[irow+1]; tempPos++)
            {
                jaRes[pos] = jaIn[tempPos];
                coefRes[pos] = coefIn[tempPos];
                pos++;
            }
            iatRes[N+irow+1] = pos;
        }
    }
    //---
    static CSRMAT concatDiagInCol(VECTOR &diag, CSRMAT &mat, int flag = -1)
    {
        if (flag == 1)
        {
            if (mat.nCol != diag.length) throw_line("ERROR: Wrong dimensions in matrix concatenation\n");
        }
        else 
        {
            if (flag <= 0) flag = 0;
        }
        std::shared_ptr<int[]> iatIn = mat.iat; std::shared_ptr<int[]> jaIn = mat.ja; std::shared_ptr<prec[]> coefIn = mat.coef; int nRowIn = mat.nRow;
        int N = diag.length;

        int nCol = mat.nCol;
        int nTermRes = mat.nTerm + N; int nRow = nRowIn + N;
        CSRMAT resMat(nRow, nCol, nTermRes);

        std::shared_ptr<prec[]> coefRes = resMat.coef;
        std::shared_ptr<int[]> iatRes = resMat.iat;
        std::shared_ptr<int[]> jaRes  = resMat.ja;

        int pos = 0;
        iatRes[0] = 0;
        for (int irow = 0; irow < N; irow++)
        {
            jaRes[pos] = irow + flag;
            coefRes[pos] = diag[irow];
            pos++;
            iatRes[irow+1] = pos;
        }
        //
        for (int irow = 0; irow < nRowIn; irow++)
        {
            for (int tempPos = iatIn[irow]; tempPos < iatIn[irow+1]; tempPos++)
            {
                jaRes[pos] = jaIn[tempPos];
                coefRes[pos] = coefIn[tempPos];
                pos++;
            }
            iatRes[N+irow+1] = pos;
        }
        return resMat;
    }
    //----------------
    static CSRMAT concatInRow(CSRMAT &mat1, CSRMAT &mat2)
    {
        if (mat1.nRow != mat2.nRow) throw_line("ERROR: Wrong dimensions in matrix concatenation\n");
        std::shared_ptr<int[]> iat1 = mat1.iat; std::shared_ptr<int[]> ja1 = mat1.ja; std::shared_ptr<prec[]> coef1 = mat1.coef; int nCol1 = mat1.nCol;
        std::shared_ptr<int[]> iat2 = mat2.iat; std::shared_ptr<int[]> ja2 = mat2.ja; std::shared_ptr<prec[]> coef2 = mat2.coef; int nCol2 = mat2.nCol;

        int nTermRes = mat1.nTerm+mat2.nTerm;
        int nRow = mat1.nRow; int nCol = nCol1 + nCol2;

        CSRMAT resMat(nRow, nCol, nTermRes);
        std::shared_ptr<prec[]> coefRes = resMat.coef;
        std::shared_ptr<int[]> iatRes = resMat.iat;
        std::shared_ptr<int[]> jaRes  = resMat.ja;

        int pos = 0;
        iatRes[0] = 0;
        for (int irow = 0; irow < nRow; irow++)
        {
            for (int tempPos = iat1[irow]; tempPos < iat1[irow+1]; tempPos++)
            {
                jaRes[pos] = ja1[tempPos];
                coefRes[pos] = coef1[tempPos];
                pos++;
            }
            for (int tempPos = iat2[irow]; tempPos < iat2[irow+1]; tempPos++)
            {
                jaRes[pos] = ja2[tempPos] + nCol1;
                coefRes[pos] = coef2[tempPos];
                pos++;
            }

            iatRes[irow+1] = pos;
        }
        return resMat;
    }
    //----------------------------
     static CSRMAT concatInCol(CSRMAT &mat1, CSRMAT &mat2)
    {
        if (mat1.nCol != mat2.nCol) throw_line("ERROR: Wrong dimensions in matrix concatenation\n");
        std::shared_ptr<int[]> iat1 = mat1.iat; std::shared_ptr<int[]> ja1 = mat1.ja; std::shared_ptr<prec[]> coef1 = mat1.coef; int nRow1 = mat1.nRow;
        std::shared_ptr<int[]> iat2 = mat2.iat; std::shared_ptr<int[]> ja2 = mat2.ja; std::shared_ptr<prec[]> coef2 = mat2.coef; int nRow2 = mat2.nRow;

        int nTermRes = mat1.nTerm+mat2.nTerm;
        int nRow = nRow1 + nRow2; int nCol = mat1.nCol;

        CSRMAT resMat(nRow, nCol, nTermRes);
        std::shared_ptr<prec[]> coefRes = resMat.coef;
        std::shared_ptr<int[]> iatRes = resMat.iat;
        std::shared_ptr<int[]> jaRes  = resMat.ja;

        int pos = 0;
        iatRes[0] = 0;
        for (int irow = 0; irow < nRow1; irow++)
        {
            for (int tempPos = iat1[irow]; tempPos < iat1[irow+1]; tempPos++)
            {
                jaRes[pos] = ja1[tempPos];
                coefRes[pos] = coef1[tempPos];
                pos++;
            }
            iatRes[irow+1] = pos;
        }
        //---
        for (int irow = 0; irow < nRow2; irow++)
        {
            for (int tempPos = iat2[irow]; tempPos < iat2[irow+1]; tempPos++)
            {
                jaRes[pos] = ja2[tempPos];
                coefRes[pos] = coef2[tempPos];
                pos++;
            }
            iatRes[nRow1+irow+1] = pos;
        }
        return resMat;
    }
    //-------------------------------------------------------
    // DELETE ROW
    //-------------------------------------------------------
    static CSRMAT deleteRow(CSRMAT &mat, int iRow)
    {
        int nRows = mat.nRow;
        int nCols = mat.nCol;
        int nTerms = mat.nTerm;
        std::shared_ptr<int[]> ja = mat.ja;
        std::shared_ptr<int[]> iat = mat.iat;
        std::shared_ptr<prec[]> coef = mat.coef; 

        int nElDel = iat[iRow+1]-iat[iRow];
        int nRedTerm = nTerms - nElDel;
        std::shared_ptr<int[]> redIat(new int[nRows]);
        std::shared_ptr<int[]> redJa(new int[nRedTerm]);
        std::shared_ptr<prec[]> redCoef(new prec[nRedTerm]);

        int countIat = 0;
        int countJa = 0;
        int nElInRow = 0;
        redIat[0] = 0;
        for (int i = 0; i < nRows ; i++)
        {
            if (i != iRow)
            {
                nElInRow = iat[i+1] - iat[i];
                for (int j = iat[i]; j < iat[i+1]; j++ )
                {
                    redJa[countJa] = ja[j];
                    redCoef[countJa] = coef[j];
                    countJa++;
                }
                redIat[countIat + 1] = redIat[countIat] + nElInRow;
                countIat ++;
            }                       
        }
        CSRMAT res(nRows-1, nCols, redIat, redJa, redCoef, nRedTerm);
        return res;
    }


    
    //-------------------------------------------------------
    // DELETE COL
    //-------------------------------------------------------
    static CSRMAT deleteCol(CSRMAT &mat, int iCol)
    {
        int nRows = mat.nRow;
        int nCols = mat.nCol;
        int nTerms = mat.nTerm;
        std::shared_ptr<int[]> ja = mat.ja;
        std::shared_ptr<int[]> iat = mat.iat;
        std::shared_ptr<prec[]> coef = mat.coef; 

        int* redIatP = (int*) calloc((nRows+1), sizeof(int));
        int* redTempJaP  = (int*) calloc(nTerms, sizeof(int));
        prec* redTempCoefP = (prec*) calloc(nTerms, sizeof(prec));

        std::shared_ptr<int[]> redIat = std::shared_ptr<int[]>(new(redIatP) int[nRows+1], [](int* redIatP){free(redIatP);});
        std::shared_ptr<int[]> redTempJa = std::shared_ptr<int[]>(new(redTempJaP) int[nTerms], [](int* redTempJaP){free(redTempJaP);});
        std::shared_ptr<prec[]> redTempCoef = std::shared_ptr<prec[]>(new(redTempCoefP) prec[nTerms], [](prec* redTempCoefP){free(redTempCoefP);});

        // int countIat = 0; ?????? la segna inutilizzata
        int countJa = 0;
        int nElInRow = 0;
        int erasedCol = 0;
        redIat[0] = 0;
        for (int i = 0; i < nRows; i++)
        {
            nElInRow = 0;
            erasedCol = 0;
            for (int j = iat[i]; j < iat[i+1]; j++)
            {
                if (ja[j] != iCol)
                {
                    if (ja[j] > iCol) erasedCol = 1;
                    redTempJa[countJa] = ja[j] - erasedCol;                         
                    redTempCoef[countJa] = coef[j];
                    countJa++;
                    nElInRow++;
                } 
            }
            redIat[i+1] = redIat[i] + nElInRow;
        } 
        
        int nRedTerm = countJa; 
        std::shared_ptr<int[]> redJa(new int[countJa]);
        std::shared_ptr<prec[]> redCoef(new prec[countJa]);
        
        VECTOR::copy(redTempJa, redJa, countJa);
        VECTOR::copy(redTempCoef, redCoef, countJa);
        //VECTOR::print2(redJa, redCoef, countJa, "ja", "coef");
        //VECTOR::print(redIat, nRows+1);
        CSRMAT res(nRows, nCols-1, redIat, redJa, redCoef, nRedTerm);
        return res;
    }


    //------------------------------------------------
    // DELETE ROW & COL
    //------------------------------------------------
        void deleteId(VECTOR_INT &ids)
    {
        int startRow = ids[0];
        int N = ids.length;

        int countRow = 1;
        int pos = 0;
        int countDelInRow;
        // rows up to first
        //std::cout << "\n up to first row \n";
        for (int oldrow = 0; oldrow < startRow; oldrow++)
        {
            for (int oldPos = iat[oldrow]; oldPos < iat[oldrow+1]; oldPos++)
            {
                int oldCol = ja[oldPos];
                bool isIn = ids.hasInOrd(oldCol, countDelInRow);
                if (isIn) continue;
                ja[pos] =  oldCol - countDelInRow;
                coef[pos] = coef[oldPos];
                //std::cout << "\n pos: " << pos << "\n";
                pos++;
            }
            //std::cout << "\n countRow: " << countRow << "\n";
            iat[countRow] = pos;
            countRow++;
        }
        // rows from first to last
        //std::cout << "\n from first to last ids row \n";
        for (int idsrow = 0; idsrow < N-1; idsrow++)
        {
            int currRowToDel = ids[idsrow]+1;
            int nextRowToDel = ids[idsrow+1];
            for (int oldrow = currRowToDel; oldrow < nextRowToDel; oldrow++)
            {
                for (int oldPos = iat[oldrow]; oldPos < iat[oldrow+1]; oldPos++)
                {
                    int oldCol = ja[oldPos];
                    bool isIn = ids.hasInOrd(oldCol, countDelInRow);
                    if (isIn) continue;
                    ja[pos] =  oldCol - countDelInRow;
                    coef[pos] = coef[oldPos];
                    //std::cout << "\n pos: " << pos << "\n";
                    pos++;
                }
                //std::cout << "\n countRow: " << countRow << "\n";
                iat[countRow] = pos;
                countRow++;
            }
            
        }
        // rows from last to nRow
        //std::cout << "\n from last ids to nRow\n";
        for (int oldRow = ids[N-1]+1; oldRow < nRow; oldRow++)
        {
            for (int oldPos = iat[oldRow]; oldPos < iat[oldRow+1]; oldPos++)
            {
                int oldCol = ja[oldPos];
                bool isIn = ids.hasInOrd(oldCol, countDelInRow);
                if (isIn) continue;
                ja[pos] =  oldCol - countDelInRow;
                coef[pos] = coef[oldPos];
                //std::cout << "\n pos: " << pos << "\n";
                pos++;
            }
            //std::cout << "\n countRow: " << countRow << "\n";
            iat[countRow] = pos;
            countRow++;
        }
        
        // VECTOR_INT::shrink(ja, pos); 
        // VECTOR::shrink(coef, pos);
        // VECTOR_INT::shrink(iat, countRow);
        nRow = countRow-1; nCol = countRow-1;
        nTerm = pos;
        // printf("pos: %d\n", pos);
        // VECTOR::print2(ja, coef, nTerm, "ja", "coef");
        // VECTOR::print(iat, nRow+1);
        //---
    }
    //---
    void deleteRow(VECTOR_INT &ids)
    {
        int startRow = ids[0];
        int N = ids.length;
        // delete rows
        // recalibrate ja 
        int pos = iat[startRow];
        for (int idsrow = 0; idsrow < N-1; idsrow++)
        {
            int currRowToDel = ids[idsrow]+1;
            int nextRowToDel = ids[idsrow+1];
            for (int oldPos = iat[currRowToDel]; oldPos < iat[nextRowToDel]; oldPos++)
            {
                ja[pos] = ja[oldPos]; 
                coef[pos] = coef[oldPos];
                pos++;
            }
        }
        for (int oldPos = iat[ids[N-1]+1]; oldPos < nTerm; oldPos++)
        {
            ja[pos] = ja[oldPos]; 
            coef[pos] = coef[oldPos];
            pos++;
        }
        VECTOR_INT::shrink(ja, pos); VECTOR::shrink(coef, pos);
        // recalibrate iat
        int row = startRow+1;
        int count = iat[startRow];
        for (int idsrow = 0; idsrow < N-1; idsrow++)
        {
            for (int irow = ids[idsrow]+1; irow < ids[idsrow+1]; irow++)
            {
                iat[row] = count + iat[irow+1] - iat[irow];
                row++;
            }
        }
        for (int irow = ids[N-1]+1; irow < nRow-1; irow++)
        {
            count +=  iat[irow+1] - iat[irow];
            iat[row] = count;
            row++;
        }
        iat[row] = pos;
        VECTOR_INT::shrink(iat, row+1);
        nRow = row;
        nTerm = pos;
        //---
    }
    //---


    static CSRMAT deleteID(CSRMAT &mat, int iRow, int iCol)
    {
        CSRMAT tempMat  = CSRMAT::deleteRow(mat, iRow);
        CSRMAT redMat   = CSRMAT::deleteCol(tempMat, iCol);
        return redMat;
    }
    
    //---------------------------------
    // GET FULL MATRIX
    //---------------------------------
    // static MATRIX getFullCSR(CSRMAT mat)
    // {
    //     int nRows = mat.nRow;
    //     int nCols = mat.nCol;
    //     int nTerms = mat.nTerm;
    //     std::shared_ptr<int[]> ja = mat.ja;
    //     std::shared_ptr<int[]> iRow = getIRow(mat);
    //     std::shared_ptr<prec[]> coef = mat.coef;
    //     MATRIX fullMat;
    //     fullMat.zeros(nRows, nCols);
    //     for(int i = 0; i < nTerms; i++)
    //     {
    //         fullMat[iRow[i]][ja[i]] = coef[i];
    //     }
    //     return fullMat;
    // }
    //---------------------------------
    // PRINT ON FILE
    //---------------------------------
    void writePatt(const char* outFileName)
    {
        // define local variables
        int beginRow, endRow;
        const char* frmt = "%5ld, %5ld, %5ld \n";
        FILE*  outFile = fopen(outFileName, "w");
        //-------------
        // set handles
        for (int i = 0; i < nRow; i++)
        {
            beginRow = iat[i];
            endRow   = iat[i+1];
            for (int j = beginRow; j < endRow; j++) fprintf(outFile, frmt, i, ja[j], 1);
        }
    }
    //---
    void write(const char* outFileName)
    {
        int beginRow, endRow;
        const char* frmt = "%5ld %5ld %20.5" format_e " \n";

        FILE* outFile = fopen(outFileName,"w");

        // Print matrix;
        for (int i = 0; i < nRow; i++)
        {
            beginRow = iat[i];
            endRow   = iat[i+1];
            for (int j = beginRow; j < endRow; j++)
            {
                fprintf(outFile, frmt, i, ja[j], coef[j]);
            }
        }
        fclose(outFile);
    }
    //---------------------------------
    // PRINT FULL MATRIX
    //---------------------------------
    void printFullCSR(const char* outFileName)
    {

        FILE* outFile = fopen(outFileName,"w");
        int pos = 0;
        fprintf(outFile, "%d %d\n", nRow, nCol);
        for (int irow = 0; irow < nRow; irow++)
        {
            for (int icol = 0; icol < nCol; icol++)
            {
                if (ja[pos] == icol && pos < iat[irow+1])
                {
                    fprintf(outFile, "%18.15" format " ", coef[pos]);
                    pos++;
                }
                else fprintf(outFile, "%1.0f ", 0.0);
            }
            fprintf(outFile, "\n");
        }
        fclose(outFile);
    }
    //---
    void printFullCSR()
    {
        int pos = 0;
        printf("\n---------------------\n");
        for (int irow = 0; irow < nRow; irow++)
        {
            for (int icol = 0; icol < nCol; icol++)
            {
                if (ja[pos] == icol && pos < iat[irow+1])
                {
                    printf("%8.4" format " ", coef[pos]);
                    pos++;
                }
                else printf("%8.0f ", 0.0);
            }
            printf("\n");
        }
        printf("---------------------\n");
    }
    static void printFullCSR(CSRMAT mat)
    {
        mat.printFullCSR();
    }
    //------------------------------------
    // DELETE
    //------------------------------------
    void dlt()
    {
        nRow = 0; nCol = 0; nTerm = 0;
        iat.reset(); ja.reset(); coef.reset();
    }
    //------------------------------------------------------------------
    //------------------------------------------------------------------------
    //SORT THE COLUMNS(TAKEN ONCE) AND SUMMING RESPECTIVE COEFFS FOR EACH ROW
    //------------------------------------------------------------------------
    static int sort(std::shared_ptr<int[]> &rowCol, std::shared_ptr<prec[]> &rowCoeff, int N)
    {
        // allocate space to sum the non zero coeff relative to the same column
        std::vector<prec> sumCoefByCol(N);
        for (int i = 0; i < N; i++) sumCoefByCol[i] = 0;
        // prec* sumCoefByCol = (prec*) calloc(N, sizeof(prec));

        std::shared_ptr<int[]> effectiveCol(new int[N]);

        // fulfil the pointers
        int nNZCol = 0;
        int index;
        for (int i = 0; i < N; i++)
        {
            int locCol = rowCol[i];
            index = nNZCol;
            bool check = VECTOR::isIn(effectiveCol, nNZCol, locCol, index);
            if (!check)
            {
                effectiveCol[nNZCol] = locCol;
                nNZCol ++;
            } 
            sumCoefByCol[index] += rowCoeff[i];
        }
        int repeatedColInRow = N - nNZCol;

        
        // sort the column in the "reduced" row and make the same movements for the relative coeff 
        std::shared_ptr<int[]> order = INSERTION_SORT::insertionSortWithCoeff(effectiveCol, nNZCol);
        
        for (int i = 0; i < nNZCol; i++) 
        {
            rowCoeff[i] = sumCoefByCol[order[i]]; 
            rowCol[i] = effectiveCol[i];
        }
        
        return repeatedColInRow;
    }
    //---
    static int sort(int* &rowCol, prec* &rowCoeff, int N)
    {
        // allocate space to sum the non zero coeff relative to the same column
        std::vector<prec> sumCoefByCol(N);
        for (int i = 0; i < N; i++) sumCoefByCol[i] = 0;
        // prec* sumCoefByCol = (prec*) calloc(N, sizeof(prec));

        std::shared_ptr<int[]> effectiveCol(new int[N]);

        // fulfil the pointers
        int nNZCol = 0;
        int index;
        for (int i = 0; i < N; i++)
        {
            int locCol = rowCol[i];
            index = nNZCol;
            bool check = VECTOR::isIn(effectiveCol, nNZCol, locCol, index);
            if (!check)
            {
                effectiveCol[nNZCol] = locCol;
                nNZCol ++;
            } 
            sumCoefByCol[index] += rowCoeff[i];
        }
        int repeatedColInRow = N - nNZCol;

        
        // sort the column in the "reduced" row and make the same movements for the relative coeff 
        std::shared_ptr<int[]> order = INSERTION_SORT::insertionSortWithCoeff(effectiveCol, nNZCol);
        
        for (int i = 0; i < nNZCol; i++) 
        {
            rowCoeff[i] = sumCoefByCol[order[i]]; 
            rowCol[i] = effectiveCol[i];
        }

        return repeatedColInRow;
    }
    //---

    
};





// PRODOTTO MATRICI
// if (nCol != mat.nRow) throw_line("ERROR: SPARSE MATRIX PRODUCT, inconstistend dimensions.\n")
//         // int nRowIn = mat.nRow;
//         int nColIn = mat.nCol;
//         std::shared_ptr<int[]> iatIn = mat.iat;
//         std::shared_ptr<int[]> jaIn = mat.ja; 
//         int nTermIn = mat.nTerm;
//         std::shared_ptr<prec[]> coefIn = mat.coef;
//         //----------------------
//         int rescale = nTerm;
//         int incumbMaxStor = 8*nTerm;
//         std::shared_ptr<int[]> iatRes(new int[nRow+1]);
//         VECTOR_INT jaRes(incumbMaxStor);
//         VECTOR coefRes(incumbMaxStor);

//         int pos = 0;
//         int posRes = 0;
//         iatRes[0] = 0;
//         for (int irow = 0; irow < nRow; irow++)
//         {
//             int nElInRow = iat[irow+1] - iat[irow];
//             VECTOR_INT tempColInRow(nElInRow);
//             VECTOR tempCoefInRow(nElInRow);
//             for (int el = 0; el < nElInRow; el++)
//             {
//                 tempColInRow[el]  = ja[pos];
//                 tempCoefInRow[el] = coef[pos];
//                 pos++;
//             }
//             VECTOR_INT colRes = VECTOR_INT::zeros(nColIn);
//             VECTOR tempCoefRes = VECTOR::zeros(nColIn);
//             VECTOR_INT nnzPos(nColIn); 
//             int nnz = 0;
//             //------------------------------------
//             for (int iter = 0; iter < nElInRow; iter++)
//             {
//                 int jrow = tempColInRow[iter];
//                 for (int posIn = iatIn[jrow]; posIn < iatIn[jrow+1]; posIn++)
//                 {
//                     int tempColIn = jaIn[posIn];
//                     if (colRes[tempColIn] == 0)
//                     {
//                         colRes[tempColIn] = 1;
//                         nnzPos[nnz] = tempColIn;
//                         nnz++;
//                     }
//                     tempCoefRes[tempColIn] += tempCoefInRow[iter] * coefIn[posIn];
//                 }
//             }
//             VECTOR_INT trueRow(nnz); VECTOR trueCoef(nnz);
//             for (int k = 0; k < nnz; k++) 
//             {
//                 int tempCol = nnzPos[k]; 
//                 trueRow[k] = tempCol;
//                 trueCoef[k] = tempCoefRes[tempCol];
//             }  
//             int repitedCol = sort(trueRow.P, trueCoef.P, nnz);
//             for (int tempPos = 0; tempPos < nnz; tempPos++)
//             {
//                 jaRes[posRes]   = trueRow[tempPos];
//                 coefRes[posRes] = trueCoef[tempPos];
//                 posRes++;
//                 if (posRes == incumbMaxStor)
//                 {
//                     incumbMaxStor += rescale;
//                     jaRes.enlarge(incumbMaxStor);
//                     coefRes.enlarge(incumbMaxStor);
//                 }
//             }
//             iatRes[irow+1] = posRes;
//             // free
//             //-----
//         }
//         jaRes.shrink(posRes); coefRes.shrink(posRes);
//         CSRMAT resMat(nRow, nColIn, iatRes, jaRes.P, coefRes.P, posRes);
//         return resMat;