#pragma once

#include "../../CODE_HEADERS/codeHeader.h"

class MATRIX
{
    public:
    std::shared_ptr<prec*[]> PP = 0;
    std::shared_ptr<prec[]> P = 0;
    int nRow = 0;
    int nCol = 0;

    //------------------------------------------------
    // CONSTRUCTORS
    //------------------------------------------------
public:
    // EMPTY CONSTRUCTOR
    MATRIX() 
    {
        PP = 0;
        P = 0;
    }
    //---
    MATRIX(int nRows, int nCols) // standard initialization !! TO DO !!!!!!!!!!!
    {
        nRow = nRows;
        nCol = nCols;
        allocate(nRows, nCols, PP, P);    
    }

    MATRIX(int nRows, int nCols, std::shared_ptr<prec*[]>  &matPP, std::shared_ptr<prec[]>  &matP)
    {
        nRow = nRows;
        nCol = nCols;
        PP   = matPP;
        P    = matP;
        for (int i = 0; i < nRows; i++)
        {
           if (matPP[i] != &matP[i*nCols]) throw_line("ERROR: Allocating column matrix in row format.\n");
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
    //---
    void define(int nRows, int nCols, std::shared_ptr<prec*[]>  &matPP, std::shared_ptr<prec[]>  &matP)
    {
        nRow = nRows;
        nCol = nCols;
        PP   = matPP;
        P    = matP;
    }
    //---
    void define(MATRIX &mat)
    {
        initialize(mat.nRow, mat.nCol);
        P = mat.P;
        PP = mat.PP;
    }
    //---
    void resetZeros()
    {
        for (int irow = 0; irow < nRow; irow++)
        {
            for (int icol = 0; icol < nCol; icol++)
            {
                PP[irow][icol] = 0;
            }
        }
    }
    //---
    void zeros(int nRows, int nCols)
    {
        nRow = nRows;
        nCol = nCols;
        allocateZero(nRows, nCols, PP, P);
    }
    //---
    void setZeros(int nRows, int nCols)
    {
        nRow = nRows;
        nCol = nCols;
        allocateZero(nRows, nCols, PP, P);
    }
    
    //-----------------------------------------
    // ALLOCATE POINTERS
    //-----------------------------------------
    static void allocate(int nRows, int nCols, std::shared_ptr<prec*[]>  &Mat, std::shared_ptr<prec[]>  &MatBuff)
    {
        int nEntries = nRows * nCols;

        MatBuff = VECTOR::makePointer(nEntries);
        Mat     = VECTOR::makeDoublePointer(nRows);

        int index = 0;
        for (int i = 0; i < nRows; i++)
        {
            Mat[i] = &MatBuff[index]; 
            index += nCols;
            if (Mat[i] == nullptr) throw_line("ERROR IN FULL MATRIX ALLOCATION \n");
        } 
    }
    //---
    static void allocateZero(int nRows, int nCols, std::shared_ptr<prec*[]> &Mat, std::shared_ptr<prec[]> &MatBuff)
    {
        int Nentries = nRows * nCols;

        prec* MatBuff_p = (prec*) calloc(Nentries, sizeof(prec));

        MatBuff = std::shared_ptr<prec[]>(new(MatBuff_p) prec[Nentries], [](prec* MatBuff_p){free(MatBuff_p);});
        Mat     = VECTOR::makeDoublePointer(nRows);

        int index = 0;
        for (int i = 0; i < nRows; i++)
        {
            Mat[i] = &MatBuff[index]; 
            index += nCols;
            if (Mat[i] == nullptr) throw_line("ERROR IN FULL MATRIX ALLOCATION \n");
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

    VECTOR get_row(int rowId)
    {
        VECTOR row(nCol);
        row = PP[rowId];
        return row;
    }

    void get_row(int rowId, VECTOR &row)
    {
        if (row.length != nCol)
        {
            throw_line("ERROR: trying to copy a row in a vector with length different from the number of columns\n");
        }
        row = PP[rowId];
    }

    static void get_row(int rowId, MATRIX &mat, VECTOR &row)
    {
        if (row.length != mat.nCol)
        {
            throw_line("ERROR: trying to copy a row in a vector with length different from the number of columns\n");
        }
        row = mat.PP[rowId];
    }

    //-----------------------------------------
    // GET COL COPY
    //-----------------------------------------
    std::shared_ptr<prec[]> operator () (int colID)
    {
        std::shared_ptr<prec[]> Col(new prec[nRow]);

        for (int i = 0; i < nRow; i++) Col[i] = PP[i][colID];

        return Col;
    }

    //---------------------------------------------------------
    // MODIFY COLUMN
    //---------------------------------------------------------
    void modifyColumn(int colID, prec* newCol)
    {
        for (int i = 0; i < nRow; i++)
        {
            PP[i][colID] = newCol[i];
        }
    }
    //--
    void modifyColumn(int colID, VECTOR &newCol)
    {
        std::shared_ptr<prec[]> P_in = newCol.P;
        for (int i = 0; i < nRow; i++)
        {
            PP[i][colID] = P_in[i];
        }
    }
    //---------------------------------------------------------
    // MODIFY ROW
    //---------------------------------------------------------
    void modifyRow(int rowID, std::shared_ptr<prec[]> &newRow)
    {
        for (int i = 0; i < nCol; i++)
        {
            PP[rowID][i] = newRow[i];
        }
    }
    //---
    void modifyRow(int rowID, VECTOR &newRow)
    {
        if (newRow.length != nCol) throw_line("ERROR: Modifying a row with a vector with wrong size w.r.t nCol\n");
        std::shared_ptr<prec[]> P_in = newRow.P;
        for (int i = 0; i < nCol; i++)
        {
            PP[rowID][i] = P_in[i];
        }
    }
    
    //---------------------------------
    // MATRIX SELF OPERATORS
    //---------------------------------
    void operator += (MATRIX &mat)
    {
        if (nRow != mat.nRow  || nCol != mat.nCol) throw_line("ERROR: Matrix sizes are not the same");
        for (int i = 0; i < nRow; i++)
        {
            for (int j = 0; j < nCol; j++)
            {
                PP[i][j] += mat[i][j];
            }
        }
    }
    //---
    void operator -= (MATRIX &mat)
    {
        if (nRow != mat.nRow  || nCol != mat.nCol) throw_line("ERROR: Matrix sizes are not the same");
        for (int i = 0; i < nRow; i++)
        {
            for (int j = 0; j < nCol; j++)
            {
                PP[i][j] -= mat[i][j];
            }
        }
    }
    //---
    void operator *= (MATRIX &mat)
    {
        int nCols = mat.nCol;
        int nRows = mat.nRow;
        if (nCol != nRows) throw_line("ERROR: Matrix sizes are not compatible");
        std::shared_ptr<prec*[]> PP_in = mat.PP;
        MATRIX res(nRow, nCols);       
        std::shared_ptr<prec*[]> tempResPP = res.PP;
        for (int i = 0; i < nRow; i++)
        {            
            for (int k = 0; k < nCols; k++)
            {
                prec tempRes = 0;
                for (int j = 0; j< nCol; j++)
                {
                    tempRes += PP[i][j]*PP_in[j][k];
                }
                tempResPP[i][k] = tempRes;
            }
        }
        // dlt();
        P.reset(); PP.reset();
        nCol = nCols;
        P  = res.P;
        PP = tempResPP;
    }
    //---
    void operator *= (prec coef)
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
    void operator *= (VECTOR &coef)
    {
        if (coef.length != nRow) throw_line("ERROR: incmpatible division of a matrix by a vector\n");
        for (int i = 0; i < nRow; i++)
        {
            prec temp_coef = coef[i];
            for (int j = 0; j < nRow; j++)
            {
                PP[i][j] *= temp_coef;
            }
        }
    }
    //---
    void operator /= (prec coef)
    {
        for (int i = 0; i < nRow; i++)
        {
            for (int j = 0; j < nRow; j++)
            {
                PP[i][j] /= coef;
            }
        }
    }
    //---
    void operator /= (VECTOR &coef)
    {
        if (coef.length != nRow) throw_line("ERROR: incmpatible division of a matrix by a vector\n");
        for (int i = 0; i < nRow; i++)
        {
            prec temp_coef = coef[i];
            if (abs(temp_coef) <= 1e-14) throw_line("ERROR: dividing by zero\n");
            for (int j = 0; j < nRow; j++)
            {
                PP[i][j] /= temp_coef;
            }
        }
    }

    //---------------------------------
    // SUM MATRICES
    //---------------------------------
    MATRIX operator + (MATRIX &mat)
    {
        if (nRow != mat.nRow  || nCol != mat.nCol) throw_line("ERROR: Matrix sizes are not the same");
        MATRIX res(nRow, nCol);
        std::shared_ptr<prec*[]> PP_in = mat.PP;
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
    MATRIX operator - (MATRIX &mat)
    {
        if (nRow != mat.nRow  || nCol != mat.nCol) throw_line("ERROR: Matrix sizes are not the same");
        MATRIX res(nRow, nCol);
        std::shared_ptr<prec*[]> PP_in = mat.PP;
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
    // MULTIPLY MATRIX BY VECTOR
    //---------------------------------
    VECTOR operator * (VECTOR &vect)
    {
        if (nCol != vect.length) throw_line("ERROR: Sizes are not the same");
        std::shared_ptr<prec[]> resP(new prec[nRow]);
        for (int i = 0; i < nRow; i++)
        {
            prec tempRes = 0;
            for (int j = 0; j < nCol; j++)
            {
                tempRes += PP[i][j] * vect[j];
            }
            resP[i] = tempRes;
        }
        VECTOR resVect(resP, nRow);
        return resVect;
    }
    //---
    std::shared_ptr<prec[]> operator * (std::shared_ptr<prec[]> &vect)
    {
        std::shared_ptr<prec[]> resP(new prec[nRow]);
        for (int i = 0; i < nRow; i++)
        {
            prec tempRes = 0;
            for (int j = 0; j < nCol; j++)
            {
                tempRes += PP[i][j] * vect[j];
            }
            resP[i] = tempRes;
        }
        return resP;
    }
    //---
    void prod (VECTOR &vect, VECTOR &resVect)
    {
        if (nCol != vect.length) throw_line("ERROR: Sizes are not the same");
        std::shared_ptr<prec[]> resP  = resVect.P;
        std::shared_ptr<prec[]> vectP = vect.P;
        for (int i = 0; i < nRow; i++)
        {
            prec tempRes = 0;
            for (int j = 0; j < nCol; j++)
            {
                tempRes += PP[i][j] * vectP[j];
            }
            resP[i] = tempRes;
        }
    }
    //---
    void prod (VECTOR &vect, std::shared_ptr<prec[]> &resVect)
    {
        if (nCol != vect.length) throw_line("ERROR: Sizes are not the same");
        std::shared_ptr<prec[]> vectP = vect.P;
        for (int i = 0; i < nRow; i++)
        {
            prec tempRes = 0;
            for (int j = 0; j < nCol; j++)
            {
                tempRes += PP[i][j] * vectP[j];
            }
            resVect[i] = tempRes;
        }
    }
    //---
    void prod (std::shared_ptr<prec[]> &vect, std::shared_ptr<prec[]> &resVect)
    {
        for (int i = 0; i < nRow; i++)
        {
            prec tempRes = 0;
            for (int j = 0; j < nCol; j++)
            {
                tempRes += PP[i][j] * vect[j];
            }
            resVect[i] = tempRes;
        }
    }
    //---
    void prod (std::shared_ptr<prec[]> &vect, VECTOR &resVect)
    {
        std::shared_ptr<prec[]> resP = resVect.P;
        for (int i = 0; i < nRow; i++)
        {
            prec tempRes = 0;
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
    MATRIX operator * (MATRIX &mat)
    {
        int nCols = mat.nCol;
        int nRows = mat.nRow;
        if (nCol != nRows) throw_line("ERROR: Matrix sizes are not compatible");
        MATRIX res(nRow, nCols);
        std::shared_ptr<prec*[]> PP_in = mat.PP;
        std::shared_ptr<prec*[]> tempResPP = res.PP;
        for (int i = 0; i < nRow; i++)
        {
            
            for (int k = 0; k < nCols; k++)
            {
                prec tempRes = 0;
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
    MATRIX operator * (prec coef)
    {
        MATRIX res(nRow, nCol);
        std::shared_ptr<prec*[]> resPP = res.PP;
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
    MATRIX operator / (prec coef)
    {
        MATRIX res(nRow, nCol);
        std::shared_ptr<prec*[]> resPP = res.PP;
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
    // BASOLUTE VALUE OF THE MATRIX
    //------------------------------------------
    MATRIX absMat()
    {
        MATRIX matOut(nRow, nCol);
        copyTo(matOut);
        std::shared_ptr<prec[]> matOutCoef = matOut.P;
        int nTerm = nRow*nCol;
        for (int i = 0; i < nTerm; i++)
        {
            matOutCoef[i] = abs(P[i]);
        }
        return matOut;
    }


    void enlargeRows(int N)
    {
        std::shared_ptr<prec*[]> Mat; std::shared_ptr<prec[]> MatBuff;
        allocate(N, nCol, Mat, MatBuff);

        for (int i = 0; i < nCol*nRow; i++) MatBuff[i] = P[i];

        P = MatBuff; PP = Mat; nRow = N;
    }
    
    //------------------------------------------
    // COPY
    //------------------------------------------
    void operator = (MATRIX &mat)
    {
        // dlt();
        P.reset();
        PP.reset();
        int nRows = mat.nRow;
        int nCols = mat.nCol;
        initialize(nRows, nCols);
        //MATRIX copyM(nRows, nCols);
        std::shared_ptr<prec*[]> copyPP = mat.PP;
        for (int i = 0; i < nRows; i++)
        {
            for (int j = 0; j < nCols; j++)
            {
                PP[i][j] = copyPP[i][j];
            }
        }
    }
    //---
    void copyTo(MATRIX &mat)
    {
        int nRows = mat.nRow;
        int nCols = mat.nCol;
        if (nRow != nRows || nCol != nCols) throw_line("ERROR in cpy: Matrices have different sizes");
        for (int i = 0; i < nRows; i++)
        {
            for (int j = 0; j < nCols; j++)
            {
                mat[i][j] = PP[i][j];
            }
        }
    }
    //------------------------------------------------------
    // NORM
    //------------------------------------------------------
    prec normFro()
    {
        prec norm = 0;
        for (int i = 0; i < nRow; i++)
        {
            for (int j = 0; j  < nCol; j++)
            {
                prec tempVal = PP[i][j];
                norm += tempVal*tempVal;
            }
        }
        norm = sqrt(norm);
        return norm;
    }
    //------------------------------------------------------
    // GET TRANSPOSED MATRIX
    //------------------------------------------------------
    //-----------------
    MATRIX transpose () 
    {
        MATRIX tempM(nCol, nRow);
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
    static MATRIX transpose(MATRIX &mat)
    {
        int tempNCol = mat.nRow;
        int tempNRow = mat.nRow;
        MATRIX tempM(tempNRow, tempNCol);
        for (int i = 0; i < tempNRow; i++)
        {
            for (int j = 0; j < tempNCol; j++)
            {
                tempM[i][j] = mat[j][i];
            }
        }
        return tempM;
    }

    //-------------------------------------------------------
    // DELETE ROW
    //-------------------------------------------------------
    static MATRIX deleteRow(MATRIX &mat, int iRow)
    {
        int nCols = mat.nCol;
        int nRows = mat.nRow;
        if (iRow > nRows-1 || iRow < 0) throw_line("ERROR: Not valid Row index \n");
        int redRow = 0;
        MATRIX redMat((nRows-1), nCols);
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
    static MATRIX deleteCol(MATRIX &mat, int iCol)
    {
        int nCols = mat.nCol;
        int nRows = mat.nRow;
        if (iCol > nCols-1 || iCol < 0) throw_line("ERROR: Not valid Column index \n");
        int redCol = 0;
        MATRIX redMat(nRows, (nCols-1));
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
    static MATRIX deleteID(MATRIX &mat, int iRow, int iCol)
    {
        MATRIX tempMat  = MATRIX::deleteRow(mat, iRow);
        MATRIX redMat   = MATRIX::deleteCol(tempMat, iCol);
        // tempMat.dlt();
        tempMat.P.reset(); tempMat.PP.reset();
        return redMat;
    }

    //----------------------------------
    // SHRINK
    //----------------------------------
    void shrinkRows(int nRows)
    {
        if (nRows > nRow || nRows < 0) throw_line("ERROR: Shrinking a matrix to a bigger one\n");
        nRow = nRows;
        std::shared_ptr<prec[]> newP =VECTOR::makePointer(nRow*nCol);
        std::shared_ptr<prec*[]> newPP(new prec*[nRow]);
        VECTOR::copy(P, newP, nRow*nCol);
        for (int i = 0; i < nRow; i++)
        {
            newPP[i] = &newP[i*nCol];
        }
        P.reset(); PP.reset();
        P = newP; PP = newPP;
    }

    //----------------------------------
    // SHRINK
    //----------------------------------
    void append_row(VECTOR new_row)
    {
        if (nCol != 0)
        {
            if (new_row.length != nCol) throw_line("ERROR: Appending a vector with wrong size to a matrix\n");
        }
        else
        {
            nCol = new_row.length;
        }
        
        std::shared_ptr<prec[]> newP =VECTOR::makePointer((nRow+1)*nCol);
        std::shared_ptr<prec*[]> newPP(new prec*[nRow+1]);
        VECTOR::copy(P, newP, nRow*nCol);
        VECTOR::copy(new_row.P, newP, nCol, nRow*nCol);
        nRow = nRow+1;
        for (int i = 0; i < nRow; i++)
        {
            newPP[i] = &newP[i*nCol];
        }
        P.reset(); PP.reset();
        P = newP; PP = newPP;
    }

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


    //-------------------------------------------------------
    // DETERMINANT
    //-------------------------------------------------------
    static prec det(MATRIX &mat)
    {
        prec res = 0;
        int nRows = mat.nRow;
        int nCols = mat.nCol;
        if (nRows != nCols) throw_line("ERROR: Determinant of a non-square matrix \n");
        if (nRows == 2) 
        {
            res = mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0];
            return res;
        }
        prec* row = mat[0];
        int sgn = -1;
        for (int i = 0; i < nCols; i++)
        {
            sgn = -sgn;
            if (row[i] != 0)
            {
                MATRIX redMat = deleteID(mat, 0, i);
                res += sgn * row[i] * MATRIX::det(redMat);
            }            
        }
        return res;
    }

    bool check_nan(bool stop = false)
    {
        bool there_are_nan = false;
        for (int i = 0; i < nRow; i++)
        {
            VECTOR temp_row = get_row(i);
            bool there_is_nan = temp_row.check_nan();
            if (there_is_nan)
            {
                std::cout << "\nrow: " << i << " contains a nan\n\n";
                there_are_nan = there_is_nan;
            }
        }
        if (stop && there_are_nan)
        {
            throw_line("ERROR: vector contains not wanted nan values\n");
        }
        return there_are_nan;
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
                printf("%7.3" format " ", PP[i][j]);
            }
            printf("\n");
        }
        std::cout << "~~~~~~~~~~~ \n";
    }
    //---
    static void print(MATRIX &mat)
    {       
        int nRows = mat.nRow;
        int nCols = mat.nCol;
        std::shared_ptr<prec*[]> PP_in = mat.PP;
        std::cout << " \n~~~~~~~~~~~ \n";
        for (int i = 0; i < nRows; i++)
        {
            for (int j = 0; j < nCols; j++)
            {
                printf("%7.3" format " ", PP_in[i][j]);
            }
            printf("\n");
        }
        std::cout << "~~~~~~~~~~~ \n";

    }
    //---
    static void print(std::shared_ptr<prec*[]> &mat, int nRows, int nCols)
    {     
        std::cout << " \n~~~~~~~~~~~ \n";  
        for (int i = 0; i < nRows; i++)
        {
            for (int j = 0; j < nCols; j++)
            {
                printf("%7.3" format " ", mat[i][j]);
            }
            printf("\n");
        }
        std::cout << "~~~~~~~~~~~ \n";
    }
    //---
    static void print(std::shared_ptr<int*[]> &mat, int nRows, int nCols)
    {       
        std::cout << " \n~~~~~~~~~~~ \n";
        for (int i = 0; i < nRows; i++)
        {
            for (int j = 0; j < nCols; j++)
            {
                printf("%d ", mat[i][j]);
            }
            printf("\n");
        }
        std::cout << "~~~~~~~~~~~ \n";
    }
    //---
    void print(const char* outFileName)
    {

        std::ofstream outFile;
        outFile.open(outFileName);
        //Print matrix;
        outFile << nRow << "\t" << nCol << "\n";
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
    void print3(const char* outFileName, int N)
    {

        std::ofstream outFile;
        outFile.open(outFileName);
        //Print matrix;
        outFile << N << "\n";
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < nCol; j++)
            {
                prec value = PP[i][j]; 
                if (abs(value) < 1e-16) value = 0;
                outFile << value << "\t";
            }
            outFile << "\n" ;
        }
        outFile.close();
    }
    //---
    static void print(MATRIX &mat, const char* outFileName)
    {

        FILE* outFile = fopen(outFileName,"w");
        int nRows = mat.nRow;
        int nCols = mat.nCol;
        std::shared_ptr<prec*[]> PP_in = mat.PP;

        fprintf(outFile, "%6d \t %6d\n", nRows, nCols);
        //Print matrix;
        for (int i = 0; i < nRows; i++)
        {
            for (int j = 0; j < nCols; j++)
            {
                fprintf(outFile, "%15.12" format " ", PP_in[i][j]);
            }
            fprintf(outFile, "\n");
        }
        fclose(outFile);
    }
    //---
    static void printForMatlab(MATRIX &mat, std::string name)
    {
        int nRows = mat.nRow;
        int nCols = mat.nCol;
        std::cout << "\n";
        //std::cout << "~~~~~~~~~~~ \n";
        std::cout << name << " = [";
        for (int i = 0; i < nRows; i++)
        {
            for (int j = 0; j < nCols; j++)
            {
                printf(" %" format, mat[i][j]);
            }
            printf(";");
            if (i == (nRows - 1)) printf("]");
        }
        std::cout << "\n";
        //std::cout << "~~~~~~~~~~~ \n";
    }
    //---
    void printExp()
    {    
        std::cout << " \n~~~~~~~~~~~ \n";
        for (int i = 0; i < nRow; i++)
        {
            for (int j = 0; j < nCol; j++)
            {
                printf("%7.3" format_e " ", PP[i][j]);
            }
            printf("\n");
        }
        std::cout << "~~~~~~~~~~~ \n";
    }
    //---
   
    //----------------------------------------------------
    // DEALLOC
    //----------------------------------------------------
    // void dlt()
    // {
    //     free(PP); free(P);
    //     P = nullptr; PP = nullptr;
    // }
    // //---
    // static void dlt(MATRIX &mat)
    // {
    //     free(mat.PP); free(mat.P);
    //     mat.PP = nullptr; mat.P = nullptr;
    // }
    // //---
    // void dltNull()
    // {
    //     free(PP); free(P);
    //     PP = nullptr; P = nullptr;
    // }
    //----------------------------------------------------
};