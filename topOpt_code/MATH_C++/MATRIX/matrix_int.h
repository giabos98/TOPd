#pragma once

#include "../../CODE_HEADERS/codeHeader.h"

class MATRIX_INT
{
    public:
    std::shared_ptr<int*[]>  PP = 0;
    std::shared_ptr<int[]>  P = 0;
    int nRow;
    int nCol;

    //------------------------------------------------
    // CONSTRUCTORS
    //------------------------------------------------
public:
    // EMPTY CONSTRUCTOR
    MATRIX_INT() 
    {
        PP = 0;
        P = 0;
    }
    //---

    MATRIX_INT(int nRows, int nCols) // standard initialization !! TO DO !!!!!!!!!!!
    {
        nRow = nRows;
        nCol = nCols;
        allocate(nRows, nCols, PP, P);    
    }

    MATRIX_INT(int nRows, int nCols, std::shared_ptr<int*[]>  &matPP, std::shared_ptr<int[]>  &matP)
    {
        nRow = nRows;
        nCol = nCols;
        PP   = matPP;
        P    = matP;
        for (int i = 0; i < nRows; i++)
        {
           if (matPP[i] != &matP[i*nCols]) throw_line("ERROR: Allocating column MATRIX_INT in row format.\n");
        }
    }
    
    //-----------------------------------------
    // ALLOCATE SPACE
    //-----------------------------------------
    void initialize(int nRows, int nCols)
    {
        nRow = nRows;
        nCol = nCols;
        allocate(nRows, nCols, PP, P);
    }
    
    //-----------------------------------------
    // ALLOCATE POINTERS
    //-----------------------------------------
    static void allocate(int nRows, int nCols, std::shared_ptr<int*[]>  &Mat, std::shared_ptr<int[]>  &MatBuff)
    {
        int nEntries = nRows * nCols;

        MatBuff = VECTOR_INT::makePointer(nEntries);
        Mat     = VECTOR_INT::makeDoublePointer(nRows);

        int index = 0;
        for (int i = 0; i < nRows; i++)
        {
            Mat[i] = &MatBuff[index]; 
            index += nCols;
            if (Mat[i] == nullptr) throw_line("ERROR IN FULL MATRIX_INT ALLOCATION \n");
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

    //-----------------------------------------
    // GET COL COPY
    //-----------------------------------------
    std::shared_ptr<int[]> operator () (int colID)
    {
        std::shared_ptr<int[]> Col(new int[nRow]);

        for (int i = 0; i < nRow; i++) Col[i] = PP[i][colID];

        return Col;
    }

    //---------------------------------------------------------
    // MODIFY COLUMN
    //---------------------------------------------------------
    void modifyColumn(int colID, int* &newCol)
    {
        for (int i = 0; i < nRow; i++)
        {
            PP[i][colID] = newCol[i];
        }
    }
    //--
    void modifyColumn(int colID, VECTOR_INT &newCol)
    {
        std::shared_ptr<int[]> P_in = newCol.P;
        for (int i = 0; i < nRow; i++)
        {
            PP[i][colID] = P_in[i];
        }
    }
    //---------------------------------------------------------
    // MODIFY ROW
    //---------------------------------------------------------
    void modifyRow(int rowID, int* &newRow)
    {
        for (int i = 0; i < nCol; i++)
        {
            PP[rowID][i] = newRow[i];
        }
    }
    //---
    void modifyRow(int rowID, VECTOR_INT &newRow)
    {
        std::shared_ptr<int[]> P_in = newRow.P;
        for (int i = 0; i < nCol; i++)
        {
            PP[rowID][i] = P_in[i];
        }
    }
    
    //---------------------------------
    // MATRIX_INT SELF OPERATORS
    //---------------------------------
    void operator += (MATRIX_INT &mat)
    {
        if (nRow != mat.nRow  || nCol != mat.nCol) throw_line("ERROR: MATRIX_INT sizes are not the same");
        for (int i = 0; i < nRow; i++)
        {
            for (int j = 0; j < nCol; j++)
            {
                PP[i][j] += mat[i][j];
            }
        }
    }
    //---
    void operator -= (MATRIX_INT &mat)\
    {
        if (nRow != mat.nRow  || nCol != mat.nCol) throw_line("ERROR: MATRIX_INT sizes are not the same");
        for (int i = 0; i < nRow; i++)
        {
            for (int j = 0; j < nCol; j++)
            {
                PP[i][j] -= mat[i][j];
            }
        }
    }
    //---
    void operator *= (MATRIX_INT &mat)
    {
        int nCols = mat.nCol;
        int nRows = mat.nRow;
        if (nCol != nRows) throw_line("ERROR: MATRIX_INT sizes are not compatible");
        std::shared_ptr<int*[]> PP_in = mat.PP;
        MATRIX_INT res(nRow, nCols);       
        std::shared_ptr<int*[]> tempResPP = res.PP;
        for (int i = 0; i < nRow; i++)
        {            
            for (int k = 0; k < nCols; k++)
            {
                int tempRes = 0;
                for (int j = 0; j< nCol; j++)
                {
                    tempRes += PP[i][j]*PP_in[j][k];
                }
                tempResPP[i][k] = tempRes;
            }
        }
        P.reset(); PP.reset();
        nCol = nCols;
        P  = res.P;
        PP = tempResPP;
    }
    //---
    void operator *= (int coef)
    {
        for (int i = 0; i < nRow; i++)
        {
            for (int j = 0; j < nRow; j++)
            {
                PP[i][j] *= coef;
            }
        }
    }
    //---
    void operator /= (int coef)
    {
        for (int i = 0; i < nRow; i++)
        {
            for (int j = 0; j < nRow; j++)
            {
                PP[i][j] /= coef;
            }
        }
    }

    //---------------------------------
    // SUM MATRICES
    //---------------------------------
    MATRIX_INT operator + (MATRIX_INT &mat)
    {
        if (nRow != mat.nRow  || nCol != mat.nCol) throw_line("ERROR: MATRIX_INT sizes are not the same");
        MATRIX_INT res(nRow, nCol);
        std::shared_ptr<int*[]> PP_in = mat.PP;
        for (int i = 0; i < nRow; i++)
        {
            for (int j = 0; j < nRow; j++)
            {
                res.PP[i][j] = PP[i][j] + PP_in[i][j];
            }
        }
        return res;
    }

    //---------------------------------
    // SUBTRACT MATRICES
    //---------------------------------
    MATRIX_INT operator - (MATRIX_INT &mat)
    {
        if (nRow != mat.nRow  || nCol != mat.nCol) throw_line("ERROR: MATRIX_INT sizes are not the same");
        MATRIX_INT res(nRow, nCol);
        std::shared_ptr<int*[]> PP_in = mat.PP;
        for (int i = 0; i < nRow; i++)
        {
            for (int j = 0; j < nRow; j++)
            {
                res.PP[i][j] = PP[i][j] - PP_in[i][j];
            }
        }
        return res;
    }

    //---------------------------------
    // MULTIPLY MATRIX_INT BY VECTOR_INT
    //---------------------------------
    VECTOR_INT operator * (VECTOR_INT &vect)
    {
        if (nCol != vect.length) throw_line("ERROR: Sizes are not the same");
        std::shared_ptr<int[]> resP(new int[nRow]);
        for (int i = 0; i < nRow; i++)
        {
            int tempRes = 0;
            for (int j = 0; j < nCol; j++)
            {
                tempRes += PP[i][j] * vect[j];
            }
            resP[i] = tempRes;
        }
        VECTOR_INT resVect(resP, nRow);
        return resVect;
    }
    //---
    std::shared_ptr<int[]> operator * (std::shared_ptr<int[]> &vect)
    {
        std::shared_ptr<int[]> resP(new int[nRow]);
        for (int i = 0; i < nRow; i++)
        {
            int tempRes = 0;
            for (int j = 0; j < nCol; j++)
            {
                tempRes += PP[i][j] * vect[j];
            }
            resP[i] = tempRes;
        }
        return resP;
    }
    //---
    void prod (VECTOR_INT &vect, VECTOR_INT &resVect)
    {
        if (nCol != vect.length) throw_line("ERROR: Sizes are not the same");
        std::shared_ptr<int[]> resP  = resVect.P;
        std::shared_ptr<int[]> vectP = vect.P;
        for (int i = 0; i < nRow; i++)
        {
            int tempRes = 0;
            for (int j = 0; j < nCol; j++)
            {
                tempRes += PP[i][j] * vectP[j];
            }
            resP[i] = tempRes;
        }
    }
    //---
    void prod (VECTOR_INT &vect, int* &resVect)
    {
        if (nCol != vect.length) throw_line("ERROR: Sizes are not the same");
        std::shared_ptr<int[]> vectP = vect.P;
        for (int i = 0; i < nRow; i++)
        {
            int tempRes = 0;
            for (int j = 0; j < nCol; j++)
            {
                tempRes += PP[i][j] * vectP[j];
            }
            resVect[i] = tempRes;
        }
    }
    //---
    void prod (std::shared_ptr<int[]> &vect, std::shared_ptr<int[]> &resVect)
    {
        for (int i = 0; i < nRow; i++)
        {
            int tempRes = 0;
            for (int j = 0; j < nCol; j++)
            {
                tempRes += PP[i][j] * vect[j];
            }
            resVect[i] = tempRes;
        }
    }
    //---
    void prod (std::shared_ptr<int[]> &vect, VECTOR_INT &resVect)
    {
        std::shared_ptr<int[]> resP = resVect.P;
        for (int i = 0; i < nRow; i++)
        {
            int tempRes = 0;
            for (int j = 0; j < nCol; j++)
            {
                tempRes += PP[i][j] * vect[j];
            }
            resP[i] = tempRes;
        }
    }

    //-------------------------------------
    // MULTIPLY MATRICES
    //-------------------------------------
    MATRIX_INT operator * (MATRIX_INT &mat)
    {
        int nCols = mat.nCol;
        int nRows = mat.nRow;
        if (nCol != nRows) throw_line("ERROR: MATRIX_INT sizes are not compatible");
        MATRIX_INT res(nRow, nCols);
        std::shared_ptr<int*[]> PP_in = mat.PP;
        std::shared_ptr<int*[]> tempResPP = res.PP;
        for (int i = 0; i < nRow; i++)
        {
            
            for (int k = 0; k < nCols; k++)
            {
                int tempRes = 0;
                for (int j = 0; j< nCol; j++)
                {
                    tempRes += PP[i][j]*PP_in[j][k];
                }
                tempResPP[i][k] = tempRes;
            }
        }
        return res;
    }

    //-------------------------------------
    // MULTIPLY & DIVIDE BY SCALAR
    //-------------------------------------
    MATRIX_INT operator * (int coef)
    {
        MATRIX_INT res(nRow, nCol);
        std::shared_ptr<int*[]> resPP = res.PP;
        for (int i = 0; i < nRow; i++)
        {
            for (int j = 0; j < nCol; j++)
            {
                resPP[i][j] = PP[i][j] * coef;
            }
        }
        return res;
    }
    //---
    MATRIX_INT operator / (int coef)
    {
        MATRIX_INT res(nRow, nCol);
        std::shared_ptr<int*[]> resPP = res.PP;
        for (int i = 0; i < nRow; i++)
        {
            for (int j = 0; j < nCol; j++)
            {
                resPP[i][j] = PP[i][j] / coef;
            }
        }
        return res;
    }
    
    //------------------------------------------
    // COPY
    //------------------------------------------
    void operator = (MATRIX_INT &mat)
    {
        // dlt();
        P.reset(); PP.reset();
        int nRows = mat.nRow;
        int nCols = mat.nCol;
        initialize(nRows, nCols);
        //MATRIX_INT copyM(nRows, nCols);
        std::shared_ptr<int*[]> copyPP = mat.PP;
        for (int i = 0; i < nRows; i++)
        {
            for (int j = 0; j < nCols; j++)
            {
                PP[i][j] = copyPP[i][j];
            }
        }
    }
    //---
    void copyTo(MATRIX_INT &mat)
    {
        int nRows = mat.nRow;
        int nCols = mat.nCol;
        if (nRow != nRows || nCol != nCols) throw_line("ERROR: Matrices have different sizes");
        for (int i = 0; i < nRows; i++)
        {
            for (int j = 0; j < nCols; j++)
            {
                mat[i][j] = PP[i][j];
            }
        }
    }

    //------------------------------------------------------
    // GET TRANSPOSED MATRIX_INT
    //------------------------------------------------------
    MATRIX_INT transpose () 
    {
        MATRIX_INT tempM(nCol, nRow);
        for (int i = 0; i < nCol; i++)
        {
            for (int j = 0; j < nRow; j++)
            {
                tempM[i][j] = PP[j][i];
            }
        }
        return tempM;
    }
    //---
    static MATRIX_INT transpose(MATRIX_INT &mat)
    {
        int tempNCol = mat.nRow;
        int tempNRow = mat.nRow;
        MATRIX_INT tempM(tempNRow, tempNCol);
        for (int i = 0; i < tempNRow; i++)
        {
            for (int j = 0; j < tempNCol; j++)
            {
                tempM[i][j] = mat[j][i];
            }
        }
        return tempM;
    }

    void shrinkRows(int newNRow)
    {
        std::shared_ptr<int[]> newP; std::shared_ptr<int*[]> newPP;
        allocate(newNRow, nCol, newPP, newP);
        VECTOR_INT::copy(P, newP, newNRow*nCol);
        P.reset(); PP.reset();
        P = newP; 
        PP = newPP;
        nRow = newNRow;
    }

    //-------------------------------------------------------
    // DELETE ROW
    //-------------------------------------------------------
    static MATRIX_INT deleteRow(MATRIX_INT &mat, int iRow)
    {
        int nCols = mat.nCol;
        int nRows = mat.nRow;
        if (iRow > nRows-1 || iRow < 0) throw_line("ERROR: Not valid Row index \n");
        int redRow = 0;
        MATRIX_INT redMat((nRows-1), nCols);
        for (int i = 0; i < nRows; i++)
        {
            if (i != iRow) 
            {
                for (int j = 0; j < nCols; j++)
                {
                    redMat[redRow][j] = mat[i][j];
                }
                redRow ++;
            }
        }
        return redMat;
    }

    //-------------------------------------------------------
    // DELETE COL
    //-------------------------------------------------------
    static MATRIX_INT deleteCol(MATRIX_INT &mat, int iCol)
    {
        int nCols = mat.nCol;
        int nRows = mat.nRow;
        if (iCol > nCols-1 || iCol < 0) throw_line("ERROR: Not valid Column index \n");
        int redCol = 0;
        MATRIX_INT redMat(nRows, (nCols-1));
        for (int j = 0; j < nCols; j++)
        {
            if (j != iCol) 
            {
                for (int i = 0; i < nRows; i++)
                {
                    redMat[i][redCol] = mat[i][j];
                }
                redCol ++;
            }
        }
        return redMat;
    }

    //------------------------------------------------
    // DELETE ROW & COL
    //------------------------------------------------
    static MATRIX_INT deleteID(MATRIX_INT &mat, int iRow, int iCol)
    {
        MATRIX_INT tempMat  = MATRIX_INT::deleteRow(mat, iRow);
        MATRIX_INT redMat   = MATRIX_INT::deleteCol(tempMat, iCol);
        return redMat;
    }

    //-------------------------------------------------------
    // DETERMINANT
    //-------------------------------------------------------
    static int det(MATRIX_INT mat)
    {
        int res = 0;
        int nRows = mat.nRow;
        int nCols = mat.nCol;
        if (nRows != nCols) throw_line("ERROR: Determinant of a non-square MATRIX_INT \n");
        if (nRows == 2) 
        {
            res = mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0];
            return res;
        }
        int* row = mat[0];
        int sgn = -1;
        for (int i = 0; i < nCols; i++)
        {
            if (row[i] != 0)
            {
                sgn = -sgn;
                MATRIX_INT redMat = deleteID(mat, 0, i);
                res += sgn * row[i] * MATRIX_INT::det(redMat);
            }            
        }
        return res;
    }

    //--------------------------------------------------------------------------
    // PRINT MATRIX_INT
    //--------------------------------------------------------------------------
    void print()
    {    
        std::cout << " \n~~~~~~~~~~~ \n";
        for (int i = 0; i < nRow; i++)
        {
            for (int j = 0; j < nCol; j++)
            {
                printf("%6d", PP[i][j]);
            }
            printf("\n");
        }
        std::cout << "~~~~~~~~~~~ \n";
    }
    //---
    static void print(MATRIX_INT &mat)
    {       
        int nRows = mat.nRow;
        int nCols = mat.nCol;
        std::shared_ptr<int*[]> PP_in = mat.PP;
        std::cout << " \n~~~~~~~~~~~ \n";
        for (int i = 0; i < nRows; i++)
        {
            for (int j = 0; j < nCols; j++)
            {
                printf("%6d", PP_in[i][j]);
            }
            printf("\n");
        }
        std::cout << "~~~~~~~~~~~ \n";

    }
    //---
    static void print(int** &mat, int nRows, int nCols)
    {     
        std::cout << " \n~~~~~~~~~~~ \n";  
        for (int i = 0; i < nRows; i++)
        {
            for (int j = 0; j < nCols; j++)
            {
                printf("%6d", mat[i][j]);
            }
            printf("\n");
        }
        std::cout << "~~~~~~~~~~~ \n";
    }
    //---
    void print(const char* outFileName)
    {

        FILE* outFile = fopen(outFileName,"w");
        
        fprintf(outFile, "%6d %6d\n", nRow, nCol);
        //Print MATRIX_INT;
        for (int i = 0; i < nRow; i++)
        {
            for (int j = 0; j < nCol; j++)
            {
                fprintf(outFile, "%6d", PP[i][j]);
            }
            fprintf(outFile, "\n");
        }
        fclose(outFile);
    }
    //---
    void print3(const char* outFileName)
    {

        std::ofstream outFile;
        outFile.open(outFileName);
        //Print matrix;
        outFile << nRow << "\n";
        for (int i = 0; i < nRow; i++)
        {
            for (int j = 0; j < nCol; j++)
            {
                outFile << PP[i][j] << "\t";
            }
            outFile << "\n" ;
        }
        outFile.close();
    }
    //---
    static void print(MATRIX_INT &mat, const char* outFileName)
    {

        FILE* outFile = fopen(outFileName,"w");
        int nRows = mat.nRow;
        int nCols = mat.nCol;
        std::shared_ptr<int*[]> PP_in = mat.PP;
        //Print MATRIX_INT;
        for (int i = 0; i < nRows; i++)
        {
            for (int j = 0; j < nCols; j++)
            {
                fprintf(outFile, "%6d", PP_in[i][j]);
            }
            fprintf(outFile, "\n");
        }
        fclose(outFile);
    }
    //----------------------------------------------------
    // DEALLOC
    //----------------------------------------------------
    //-------------------------------
    // RESET TO EMPTY STATE
    //-------------------------------
    void complete_reset()
    {
        nRow = 0;
        nCol = 0;
        P.reset(); P = 0;
        PP.reset(); PP = 0;
    }
    // void dlt()
    // {
    //     free(PP); free(P);
    // }
    // //---
    // static void dlt(MATRIX_INT &mat)
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