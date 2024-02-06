#pragma once

#include "../../CODE_HEADERS/codeHeader.h"

class LOWER_TRIANGULAR_INT
{
    // to get the matrix element (i,j) must be written PP[i][N-j]
    public:
    std::shared_ptr<int*[]> PP = 0;
    std::shared_ptr<int[]> P = 0;
    int N;
    int nEntries;

    //------------------------------------------------
    // CONSTRUCTORS
    //------------------------------------------------
    LOWER_TRIANGULAR_INT(int n)
    {
        N = n;
        nEntries = N*(N+1)/2;
        allocate(N, PP, P);  
    }
    //---
    LOWER_TRIANGULAR_INT(int n, std::shared_ptr<int*[]> &pp, std::shared_ptr<int[]> &p)
    {
        N = n; 
        PP = pp;
        P = p;
    }
    //------------------------------------------------------------------
    // ALLOCATE
    //------------------------------------------------------------------
    static LOWER_TRIANGULAR_INT zeros(int N)
    {
        int nEntries = N*(N+1)/2;

        std::shared_ptr<int[]> matBuff(new int[nEntries]);

        for (int i = 0; i < nEntries; i++) matBuff[i] = 0;
        std::shared_ptr<int*[]> mat(new int*[N]);

        int index = 0;
        for (int i = 0; i < N; i++)
        {
            index += i;
            mat[i] = &matBuff[index]; 
            if (mat[i] == nullptr) throw_line("ERROR IN LOWER TRIANGULAR MATRIX ALLOCATION \n");
        } 
        return LOWER_TRIANGULAR_INT(N, mat, matBuff);
    }
    //---
    static void allocate(int N, std::shared_ptr<int*[]>  &mat, std::shared_ptr<int[]> &matBuff)
    {
        int nEntries = N*(N+1)/2;
        matBuff = VECTOR_INT::makePointer(nEntries);
        mat     = VECTOR_INT::makeDoublePointer(N);

        int index = 0;
        for (int i = 0; i < N; i++)
        {
            index += i;
            mat[i] = &matBuff[index]; 
            if (mat[i] == nullptr) throw_line("ERROR IN UPPER TRIANGULAR MATRIX ALLOCATION \n");
        } 
    }
    //-----------------------------------------
    // GET ROW
    //-----------------------------------------
    int* operator [] (int rowID)
    {
        int* rowPointer = PP[rowID];
        return rowPointer;
    }

    //---------------------------------
    // MATRIX SELF OPERATORS
    //---------------------------------
    void operator *= (int coef)
    {
        for (int i = 0; i < nEntries; i++)
        {
                P[i] *= coef;
        }
    }
    //---
    void operator /= (int coef)
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
            for (int j = 0; j < i; j++) printf("%7.3d ", PP[i][j]);
            for (int j = i; j < N; j++) //col
            {
                printf("%7.1d ", 0);
            }
            printf("\n");
        }
        std::cout << "~~~~~~~~~~~ \n";
    }
    //---
    static void print(LOWER_TRIANGULAR_INT &mat)
    {   
        std::cout << " \n~~~~~~~~~~~ \n";
        std::shared_ptr<int*[]>  PP = mat.PP;
        int N = mat.N;
        for (int i = 0; i < N; i++) //row
        {
            for (int j = 0; j < i; j++) printf("%7.3d ", PP[i][j]);
            for (int j = i; j < N; j++) //col
            {
                printf("%7.1d ", 0);
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
            for (int j = 0; j < i; j++) printf("%7.3d ", PP[i][j]);
            for (int j = i; j < N; j++) //col
            {
                fprintf(outFile, "%7.0d ", 0);
            }
            fprintf(outFile, "\n");
        }
        fclose(outFile);
    }
    //---
    static void print(LOWER_TRIANGULAR_INT &mat, const char* outFileName)
    {

        FILE* outFile = fopen(outFileName,"w");
        std::shared_ptr<int*[]>  PP = mat.PP;
        int N = mat.N;
        //Print matrix;
        for (int i = 0; i < N; i++) //row
        {
            for (int j = 0; j < i; j++) printf("%7.3d ", PP[i][j]);
            for (int j = i; j < N; j++) //col
            {
                fprintf(outFile, "%7.0d ", 0);
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
    // static void dlt(LOWER_TRIANGULAR_INT &mat)
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


#pragma once

#include "../../CODE_HEADERS/codeHeader.h"

class LOWER_TRIANGULAR
{
    // to get the matrix element (i,j) must be written PP[i][N-j]
    public:
    std::shared_ptr<int_fast64_t*[]> PP = 0;
    std::shared_ptr<int_fast64_t[]> P = 0;
    int_fast64_t N;
    int_fast64_t nEntries;

    //------------------------------------------------
    // CONSTRUCTORS
    //------------------------------------------------
    LOWER_TRIANGULAR(int_fast64_t n)
    {
        N = n;
        nEntries = N*(N+1)/2;
        allocate(N, PP, P);  
    }
    //---
    LOWER_TRIANGULAR(int_fast64_t n, std::shared_ptr<int_fast64_t*[]> &pp, std::shared_ptr<int_fast64_t[]> &p)
    {
        N = n; 
        PP = pp;
        P = p;
    }
    //------------------------------------------------------------------
    // ALLOCATE
    //------------------------------------------------------------------
    static LOWER_TRIANGULAR zeros(int_fast64_t N)
    {
        int_fast64_t nEntries = N*(N+1)/2;

        std::shared_ptr<int_fast64_t[]> matBuff(new int_fast64_t[nEntries]);
        for (int_fast64_t i = 0; i < nEntries; i++) matBuff[i] = 0;
        std::shared_ptr<int_fast64_t*[]> mat(new int_fast64_t*[N]);


        int_fast64_t index = 0;
        for (int_fast64_t i = 0; i < N; i++)
        {
            index += i;
            mat[i] = &matBuff[index]; 
            if (mat[i] == nullptr) throw_line("ERROR IN LOWER TRIANGULAR MATRIX ALLOCATION \n");
        } 
        return LOWER_TRIANGULAR(N, mat, matBuff);
    }
    //---
    static void allocate(int_fast64_t N, std::shared_ptr<int_fast64_t*[]>  &mat, std::shared_ptr<int_fast64_t[]> &matBuff)
    {
        int_fast64_t nEntries = N*(N+1)/2;
        matBuff = std::shared_ptr<int_fast64_t[]>(new int_fast64_t[nEntries]);
        mat = std::shared_ptr<int_fast64_t*[]> (new int_fast64_t*[N]);

        int_fast64_t index = 0;
        for (int_fast64_t i = 0; i < N; i++)
        {
            index += i;
            mat[i] = &matBuff[index]; 
            if (mat[i] == nullptr) throw_line("ERROR IN UPPER TRIANGULAR MATRIX ALLOCATION \n");
        } 
    }
    //-----------------------------------------
    // GET ROW
    //-----------------------------------------
    int_fast64_t* operator [] (int_fast64_t rowID)
    {
        int_fast64_t* rowPointer = PP[rowID];
        return rowPointer;
    }

    //---------------------------------
    // MATRIX SELF OPERATORS
    //---------------------------------
    void operator *= (int_fast64_t coef)
    {
        for (int_fast64_t i = 0; i < nEntries; i++)
        {
                P[i] *= coef;
        }
    }
    //---
    void operator /= (int_fast64_t coef)
    {
        for (int_fast64_t i = 0; i < nEntries; i++)
        {
                P[i] /= coef;
        }
    }
};