#pragma once

#include "../../CODE_HEADERS/codeHeader.h"

class UPPER_TRIANGULAR
{
    // to get the matrix element (i,j) must be written PP[i][N-j]
    public:
    std::shared_ptr<prec*[]> PP = 0;
    std::shared_ptr<prec[]> P = 0;
    int N;
    int nEntries;

    //------------------------------------------------
    // CONSTRUCTORS
    //------------------------------------------------
    UPPER_TRIANGULAR(int n)
    {
        N = n;
        nEntries = N*(N+1)/2;
        allocate(N, PP, P);  
    }
    //------------------------------------------------------------------
    // ALLOCATE
    //------------------------------------------------------------------
    static void allocate(int N, std::shared_ptr<prec*[]> &mat, std::shared_ptr<prec[]> &matBuff)
    {
        int nEntries = N*(N+1)/2;
        matBuff = VECTOR::makePointer(nEntries);
        mat     = VECTOR::makeDoublePointer(N);

        int index = 0;
        for (int i = 0; i < N; i++)
        {
            mat[i] = &matBuff[index]; 
            index += N - i;
            if (mat[i] == nullptr) throw_line("ERROR IN UPPER TRIANGULAR MATRIX ALLOCATION \n");
        } 
    }
    //-----------------------------------------
    // GET ROW
    //-----------------------------------------
    prec* operator [] (int rowID)
    {
        prec* rowPointer = PP[rowID];
        return rowPointer;
    }
    
    //-----------------------------------------
    // MODIFY ELEMENT
    //-----------------------------------------
    void el(int rowID, int colID, prec newValue)
    {
        PP[rowID][colID-rowID] = newValue;
    }
    prec el(int rowID, int colID)
    {
        return PP[rowID][colID-rowID];
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
        if (PP == nullptr) throw_line("ERROR: Matrix PP pointer is 'nullptr' \n");
        if (P == nullptr) throw_line("ERROR: Matrix P pointer is 'nullptr' \n")
        std::cout << " \n~~~~~~~~~~~ \n";

        for (int i = 0; i < N; i++) //row
        {
            for (int j = 0; j < i; j++) printf("%7.1d ", 0);
            for (int j = 0; j < N - i; j++) //col
            {
                printf("%7.3" format " ", PP[i][j]);
            }
            printf("\n");
        }
        std::cout << "~~~~~~~~~~~ \n";
    }
    //---
    static void print(UPPER_TRIANGULAR &mat)
    {   
        std::cout << " \n~~~~~~~~~~~ \n";
        std::shared_ptr<prec*[]> PP = mat.PP;
        int N = mat.N;
        for (int i = 0; i < N; i++) //row
        {
            for (int j = 0; j < i; j++) printf("%7.1d ", 0);
            for (int j = 0; j < N - i; j++) //col
            {
                printf("%7.3" format " ", PP[i][j]);
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
        for (int i = 0; i < N; i++) //row
        {
            for (int j = 0; j < i; j++) fprintf(outFile, "%7.0d ", 0);
            for (int j = 0; j < N - i; j++) //col
            {
                printf("%7.3" format " ", PP[i][j]);
            }
            fprintf(outFile, "\n");
        }
        fclose(outFile);
    }
    //---
    static void print(UPPER_TRIANGULAR &mat, const char* outFileName)
    {

        FILE* outFile = fopen(outFileName,"w");
        std::shared_ptr<prec*[]>  PP = mat.PP;
        int N = mat.N;
        //Print matrix;
        for (int i = 0; i < N; i++) //row
        {
            for (int j = 0; j < i; j++) fprintf(outFile, "%7.0d ", 0);
            for (int j = 0; j < N - i; j++) //col
            {
                printf("%7.3" format " ", PP[i][j]);
            }
            fprintf(outFile, "\n");
        }
        fclose(outFile);
    }

    //----------------------------------------------------
    // SOLVE LINEAR SYSTEM
    //----------------------------------------------------
    std::shared_ptr<prec[]> solveLS(std::shared_ptr<prec[]>  &b)
    {
        std::shared_ptr<prec[]> res(new prec[N]);
        for (int i = N-1; i > -1; i--)
        {
            prec tempRes = b[i];
            for (int j = i+1; j < N; j++)
            {
                tempRes -= PP[i][j] * res[j];
            }
            res[i] = tempRes/PP[i][0];
        }
        return res;
    }
    //---
    void solve(VECTOR &b, VECTOR &sol)
    {
        if (N != b.length) throw_line("ERROR: Linear System sizes are not compatible \n");
        std::shared_ptr<prec[]> bP = b.P;
        std::shared_ptr<prec[]> res(new prec[N]);
        for (int i = N-1; i > -1; i--)
        {
            prec tempRes = bP[i];
            for (int j = i+1; j < N; j++)
            {
                tempRes -= el(i,j) * res[j];
            }
            res[i] = tempRes/el(i,i);
        }
        sol.length = N;
        sol.P = res;
    }
    //--
    void solveLS(std::shared_ptr<prec[]> &b, std::shared_ptr<prec[]> &sol)
    {
        if (sol == 0) sol = VECTOR::makePointer(N);
        for (int i = N-1; i > -1; i--)
        {
            prec tempRes = b[i];
            for (int j = i+1; j < N; j++)
            {
                tempRes -= el(i,j) * sol[j];
            }
            sol[i] = tempRes/el(i,i);
        }
    }
    //--
    void solveLS(std::shared_ptr<prec[]> &b, std::shared_ptr<prec[]> &sol, int nMax)
    {
        if (sol == 0) sol = VECTOR::makePointer(nMax);
        for (int i = nMax-1; i > -1; i--)
        {
            prec tempRes = b[i];
            for (int j = i+1; j < nMax; j++)
            {
                tempRes -= el(i,j) * sol[j];
            }
            sol[i] = tempRes/el(i,i);
        }
    }
    //---
    void solveLS(VECTOR &b, VECTOR &sol)
    {
        if (sol.P == 0) sol.P = VECTOR::makePointer(N);
        if (sol.length == 0) sol.length = N;
        solve(b, sol);
    }
    //--
    static void solveLS(UPPER_TRIANGULAR &mat, VECTOR &b, VECTOR &sol)
    {
        mat.solve(b, sol);
    }
    //---
    static void solveLS(UPPER_TRIANGULAR &mat, std::shared_ptr<prec[]> &b, std::shared_ptr<prec[]> &sol)
    {
        mat.solveLS(b, sol);
    }

    //----------------------------------------------------
    // DEALLOC
    //----------------------------------------------------
    // void dlt()
    // {
    //     free(PP); free(P);
    // }
    // //---
    // static void dlt(UPPER_TRIANGULAR &mat)
    // {
    //     free(mat.PP); free(mat.P);
    // }
    // //---
    // void dltNull()
    // {
    //     free(PP); free(P);
    //     PP = nullptr; P = nullptr;
    // }
};