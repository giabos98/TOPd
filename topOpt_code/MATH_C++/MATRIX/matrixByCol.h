#pragma once

#include "../../CODE_HEADERS/codeHeader.h"

class MATRIX_BY_COL
{
public:
    std::shared_ptr<prec*[]> PP = 0;
    std::shared_ptr<prec[]> P = 0;
    int nRow;
    int nCol;
    //------------------------------------------------
    // CONSTRUCTORS
    //------------------------------------------------
    MATRIX_BY_COL(int nRows, int nCols)
    {
        nRow = nRows;
        nCol = nCols;
        allocate(nRows, nCols, PP, P);  
    }

    MATRIX_BY_COL(int nRows, int nCols, std::shared_ptr<prec*[]> &matPP, std::shared_ptr<prec[]> &matP)
    {
        nRow = nRows;
        nCol = nCols;
        PP   = matPP;
        P    = matP;
        for (int i = 0; i < nCols; i++)
        {
           if (matPP[i] != &matP[i*nRows]) throw_line("ERROR: Allocating row matrix in column format.\n");
        }
    }
    
    //-----------------------------------------
    // ALLOCATE POINTERS
    //----------------------------------------- 
    static void allocate(int nRows, int nCols, std::shared_ptr<prec*[]> &Mat, std::shared_ptr<prec[]> &MatBuff)
    {
        int nEntries = nRows * nCols;

        MatBuff = VECTOR::makePointer(nEntries);
        Mat     = VECTOR::makeDoublePointer(nCols);

        int index = 0;
        for (int i = 0; i < nCols; i++)
        {
            Mat[i] = &MatBuff[index]; 
            index += nRows;
            if (Mat[i] == nullptr) throw_line("ERROR IN MATRIX BY COL ALLOCATION \n");
        } 
    }
    
    //-----------------------------------------
    // GET COLUMN
    //-----------------------------------------
    prec* operator [] (int colID)
    {
        prec* colPointer = PP[colID];
        return colPointer;
    }

    //-----------------------------------------
    // GET ROW COPY
    //-----------------------------------------
    std::shared_ptr<prec[]> operator () (int rowID)
    {
        std::shared_ptr<prec[]> Row(new prec[nCol]);

        for (int i = 0; i < nCol; i++) Row[i] = PP[i][rowID];

        return Row;
    }

    //---------------------------------------------------------
    // MODIFY COLUMN
    //---------------------------------------------------------
    void modifyCol(int colID, std::shared_ptr<prec[]> &newCol)
    {
        for (int i = 0; i < nCol; i++)
        {
            PP[colID][i] = newCol[i];
        }
    }
    //---
    void modifyCol(int colID, VECTOR &newCol)
    {
        std::shared_ptr<prec[]> P_in = newCol.P;
        for (int i = 0; i < nCol; i++)
        {
            PP[colID][i] = P_in[i];
        }
    }
    
    //---------------------------------
    // MATRIX SELF OPERATORS
    //---------------------------------
    void operator *= (prec coef)
    {
        for (int i = 0; i < nCol; i++)
        {
            for (int j = 0; j < nRow; j++)
            {
                PP[i][j] *= coef;
            }
        }
    }
    //---
    void operator /= (prec coef)
    {
        for (int i = 0; i < nCol; i++)
        {
            for (int j = 0; j < nRow; j++)
            {
                PP[i][j] /= coef;
            }
        }
    }

    //--------------------------------------------------------------------------
    // PRINT MATRIX
    //--------------------------------------------------------------------------
    void print()
    {    
        std::cout << " \n~~~~~~~~~~~ \n";
        for (int i = 0; i < nRow; i++)
        {
            for (int j = 0; j < nCol; j++)
            {
                printf("%7.3" format " ", PP[j][i]);
            }
            printf("\n");
        }
        std::cout << "~~~~~~~~~~~ \n";
    }
    //---
    static void print(MATRIX_BY_COL &mat)
    {       
        int nRows = mat.nRow;
        int nCols = mat.nCol;
        std::shared_ptr<prec*[]> PP_in = mat.PP;
        for (int i = 0; i < nRows; i++)
        {
            for (int j = 0; j < nCols; j++)
            {
                printf("%7.3" format " ", PP_in[j][i]);
            }
            printf("\n");
        }
    }
    //---
    void print(const char* outFileName)
    {

        FILE* outFile = fopen(outFileName,"w");
        
        //Print matrix;
        for (int i = 0; i < nRow; i++)
        {
            for (int j = 0; j < nCol; j++)
            {
                fprintf(outFile, "%7.3" format " ", PP[j][i]);
            }
            fprintf(outFile, "\n");
        }
        fclose(outFile);
    }
    //---
    static void print(MATRIX_BY_COL &mat, const char* outFileName)
    {

        FILE* outFile = fopen(outFileName,"w");
        int nRows = mat.nRow;
        int nCols = mat.nCol;
        std::shared_ptr<prec*[]> PP_in = mat.PP;
        //Print matrix;
        for (int i = 0; i < nRows; i++)
        {
            for (int j = 0; j < nCols; j++)
            {
                fprintf(outFile, "%7.3" format " ", PP_in[j][i]);
            }
            fprintf(outFile, "\n");
        }
        fclose(outFile);
    }
    //----------------------------------------------------
    // DEALLOC
    //----------------------------------------------------
    // void dlt()
    // {
    //     free(PP); free(P);
    // }
    // //---
    // static void dlt(MATRIX_BY_COL &mat)
    // {
    //     free(mat.PP); free(mat.P);
    // }
    // //---
    // void dltNull()
    // {
    //     free(PP); free(P);
    //     PP = nullptr; P = nullptr;
    // }
    //----------------------------------------------------
};