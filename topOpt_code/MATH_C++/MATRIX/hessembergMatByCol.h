#pragma once

#include "../../CODE_HEADERS/codeHeader.h"

class HESSEMBERG_BY_COL
{
    public:
    std::shared_ptr<prec*[]> PP = 0;
    std::shared_ptr<prec[]> P = 0;
    int m;
    int nEntries;

    //------------------------------------------------
    // CONSTRUCTORS
    //------------------------------------------------
    HESSEMBERG_BY_COL(int nCols)
    {
        m = nCols;
        nEntries = m*(m+3)/2;
        allocate(m, PP, P);  
    }
    //------------------------------------------------------------------
    static void allocate(int m, std::shared_ptr<prec*[]> &HessMat, std::shared_ptr<prec[]> &HessMatBuff)
    {
        int nEntries = m*(m+3)/2;
        HessMatBuff = VECTOR::makePointer(nEntries);
        HessMat     = VECTOR::makeDoublePointer(m);

        int index = 0;
        for (int i = 0; i < m; i++)
        {
            index += i;
            HessMat[i] = &HessMatBuff[index]; 
            index++;
            if (HessMat[i] == nullptr) throw_line("ERROR IN HESSEMBERG MATRIX ALLOCATION \n");
        } 
    }
    //-----------------------------------------
    // MODIFY ELEMENT
    //-----------------------------------------
    void el(int rowID, int colID, prec newValue)
    {
        int index = (colID )*(colID + 3)/2;
        P[index+rowID]= newValue;
    }
    //-----------------------------------------
    // GET COLUMN
    //-----------------------------------------
    prec* operator [] (int colID)
    {
        prec* colPointer = PP[colID];
        return colPointer;
    }

    //---------------------------------------------------------
    // MODIFY COLUMN
    //---------------------------------------------------------
    void modifyCol(int colID, std::shared_ptr<prec[]> &newCol)
    {
        for (int i = 0; i < colID+1; i++)
        {
            PP[colID][i] = newCol[i];
        }
    }
    //---
    void modifyCol(int colID, VECTOR &newCol)
    {
        std::shared_ptr<prec[]> P_in = newCol.P;
        for (int i = 0; i < colID+1; i++)
        {
            PP[colID][i] = P_in[i];
        }
    }
    
    //---------------------------------
    // MATRIX SELF OPERATORS
    //---------------------------------
    void operator *= (prec coef)
    {
        for (int i = 0; i < nEntries; i++)
        {
                P[i] *= coef;
        }
    }
    //---
    void operator /= (prec coef)
    {
        for (int i = 0; i < nEntries; i++)
        {
                P[i] /= coef;
        }
    }

    //--------------------------------------------------------------------------
    // PRINT MATRIX
    //--------------------------------------------------------------------------
    void print()
    {    
        std::cout << " \n~~~~~~~~~~~ \n";
        for (int j = 0; j < m; j++) printf("%7.3" format " ", PP[j][0]);
        printf("\n");

        for (int i = 1; i < m + 1; i++) //row
        {
            for (int j = 0; j < i - 1; j++) printf("%7.1d ", 0);
            for (int j = i - 1; j < m; j++) //col
            {
                printf("%7.3" format " ", PP[j][i]);
            }
            printf("\n");
        }
        std::cout << "~~~~~~~~~~~ \n";
    }
    //---
    static void print(HESSEMBERG_BY_COL &mat)
    {   
        std::cout << " \n~~~~~~~~~~~ \n";
        std::shared_ptr<prec*[]> PP = mat.PP;
        int m = mat.m;
        for (int j = 0; j < m; j++) printf("%7.3" format " ", PP[j][0]);
        printf("\n");

        for (int i = 1; i < m + 1; i++) //row
        {
            for (int j = 0; j < i - 1; j++) printf("%7.1d ", 0);
            for (int j = i - 1; j < m; j++) //col
            {
                printf("%7.3" format " ", PP[j][i]);
            }
            printf("\n");
        }
        std::cout << "~~~~~~~~~~~ \n";
    }
    //---
    void print(const char* outFileName)
    {

        FILE* outFile = fopen(outFileName,"w");
        
        //Print matrix;
        for (int j = 0; j < m; j++) fprintf(outFile, "%7.3" format " ", PP[j][0]);
        printf("\n");

        for (int i = 1; i < m + 1; i++) //row
        {
            for (int j = 0; j < i - 1; j++) fprintf(outFile, "%7.0d ", 0);
            for (int j = i - 1; j < m; j++) //col
            {
                printf("%7.3" format " ", PP[j][i]);
            }
            fprintf(outFile, "\n");
        }
        fclose(outFile);
    }
    //---
    static void print(HESSEMBERG_BY_COL &mat, const char* outFileName)
    {

        FILE* outFile = fopen(outFileName,"w");
        std::shared_ptr<prec*[]> PP = mat.PP;
        int m = mat.m;
        //Print matrix;
        for (int j = 0; j < m; j++) fprintf(outFile, "%7.3" format " ", PP[j][0]);
        printf("\n");

        for (int i = 1; i < m + 1; i++) //row
        {
            for (int j = 0; j < i - 1; j++) fprintf(outFile, "%7.0d ", 0);
            for (int j = i - 1; j < m; j++) //col
            {
                printf("%7.3" format " ", PP[j][i]);
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
    // static void dlt(HESSEMBERG_BY_COL &mat)
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