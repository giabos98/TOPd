#include "../CODE_HEADERS/codeHeader.h"

#pragma omp declare reduction(vsum:VECTOR \
                              : omp_out += omp_in) initializer(omp_priv = VECTOR::zeros(omp_orig.length))

void PARALLEL::set(VECTOR &vec, prec coef)
{
    int N = vec.length;

#pragma omp parallel for num_threads(nThread)
    for (int i = 0; i < N; i++)
        vec[i] = coef;
}

void PARALLEL::minus(VECTOR &vec)
{
    int N = vec.length;

#pragma omp parallel for num_threads(nThread)
    for (int i = 0; i < N; i++)
        vec[i] *= -1.0;
}

void PARALLEL::sum(VECTOR_INT &v1, VECTOR_INT &v2, VECTOR_INT &res)
{
    if (v1.length != v2.length)
        throw_line("ERROR: Summing vectors of different length\n");
    int vLength = v1.length;
    if (res.P == 0)
        res.initialize(vLength);
    omp_set_num_threads(nThread);
#pragma omp parallel for
    for (int i = 0; i < vLength; i++)
    {
        res[i] = v1[i] + v2[i];
    }
}

void PARALLEL::diff(VECTOR_INT &v1, VECTOR_INT &v2, VECTOR_INT &res)
{
    if (v1.length != v2.length)
        throw_line("ERROR: Summing vectors of different length\n");
    int vLength = v1.length;
    if (res.P == 0)
        res.initialize(vLength);
    omp_set_num_threads(nThread);
#pragma omp parallel for
    for (int i = 0; i < vLength; i++)
    {
        res[i] = v1[i] - v2[i];
    }
}

void PARALLEL::copy(VECTOR_INT &vecIn, VECTOR_INT &vecOut)
{
    int vLength = vecIn.length;
    if (vecOut.P == 0)
        vecOut.initialize(vLength);
    omp_set_num_threads(nThread);
#pragma omp parallel for schedule(static)
    for (int i = 0; i < vLength; i++)
    {
        vecOut[i] = vecIn[i];
    }
}

void PARALLEL::resetZeros(VECTOR_INT &vec)
{
    int length = vec.length;
    omp_set_num_threads(nThread);
#pragma omp parallel for schedule(auto)
    for (int i = 0; i < length; i++)
        vec[i] = 0;
}

void PARALLEL::sum(VECTOR &v1, VECTOR &v2, VECTOR &res)
{
    if (v1.length != v2.length)
        throw_line("ERROR: Summing vectors of different length\n");
    int vLength = v1.length;
    if (res.P == 0)
        res.initialize(vLength);
    omp_set_num_threads(nThread);
#pragma omp parallel for
    for (int i = 0; i < vLength; i++)
    {
        res[i] = v1[i] + v2[i];
    }
}
//---
void PARALLEL::sum(prec *&v1, prec *&v2, int &N, prec *&res)
{
    omp_set_num_threads(nThread);
#pragma omp parallel for
    for (int i = 0; i < N; i++)
    {
        res[i] = v1[i] + v2[i];
    }
}

void PARALLEL::diff(VECTOR &v1, VECTOR &v2, VECTOR &res)
{
    if (v1.length != v2.length)
        throw_line("ERROR: Summing vectors of different length\n");
    int vLength = v1.length;
    if (res.P == 0)
        res.initialize(vLength);
    omp_set_num_threads(nThread);
#pragma omp parallel for
    for (int i = 0; i < vLength; i++)
    {
        res[i] = v1[i] - v2[i];
    }
}
//---
void PARALLEL::diff(prec *&v1, prec *&v2, int &N, prec *&res)
{
    omp_set_num_threads(nThread);
#pragma omp parallel for
    for (int i = 0; i < N; i++)
    {
        res[i] = v1[i] - v2[i];
    }
}

void PARALLEL::minus(CSRMAT &mat)
{
    std::shared_ptr<prec[]> coef = mat.coef;
#pragma omp parallel for num_threads(nThread)
    for (int i = 0; i < mat.nTerm; i++)
    {
        coef[i] = -coef[i];
    }
}

void PARALLEL::pointProd(VECTOR &v1, VECTOR &v2, VECTOR &res)
{
    if (v1.length != v2.length)
        throw_line("ERROR: Summing vectors of different length\n");
    int vLength = v1.length;
    if (res.P == 0)
        res.initialize(vLength);
    omp_set_num_threads(nThread);
#pragma omp parallel for
    for (int i = 0; i < vLength; i++)
    {
        res[i] = v1[i] * v2[i];
    }
}
//---
void PARALLEL::pointProd(prec *&v1, prec *&v2, int &N, prec *&res)
{
    omp_set_num_threads(nThread);
#pragma omp parallel for
    for (int i = 0; i < N; i++)
    {
        res[i] = v1[i] * v2[i];
    }
}

void PARALLEL::pointDiv(VECTOR &v1, VECTOR &v2, VECTOR &res)
{
    if (v1.length != v2.length)
        throw_line("ERROR: Summing vectors of different length\n");
    int vLength = v1.length;
    if (res.P == 0)
        res.initialize(vLength);
    omp_set_num_threads(nThread);
#pragma omp parallel for
    for (int i = 0; i < vLength; i++)
    {
        res[i] = v1[i] / v2[i];
    }
}
//---
void PARALLEL::pointDiv(prec *&v1, prec *&v2, int &N, prec *&res)
{
    omp_set_num_threads(nThread);
#pragma omp parallel for
    for (int i = 0; i < N; i++)
    {
        res[i] = v1[i] / v2[i];
    }
}

void PARALLEL::multiplyByScalar(VECTOR &v1, prec &coef, VECTOR &res)
{
    int vLength = v1.length;
    if (res.P == 0)
        res.initialize(vLength);
    omp_set_num_threads(nThread);
#pragma omp parallel for
    for (int i = 0; i < vLength; i++)
    {
        res[i] = v1[i] * coef;
    }
}
//---
void PARALLEL::multiplyByScalar(prec *&v1, prec &coef, int &N, prec *&res)
{
    omp_set_num_threads(nThread);
#pragma omp parallel for
    for (int i = 0; i < N; i++)
    {
        res[i] = v1[i] * coef;
    }
}

void PARALLEL::divideByScalar(VECTOR &v1, prec &coef, VECTOR &res)
{
    if (coef == 0)
        throw_line("ERROR: Dividing by zero\n");
    int vLength = v1.length;
    if (res.P == 0)
        res.initialize(vLength);
    omp_set_num_threads(nThread);
#pragma omp parallel for
    for (int i = 0; i < vLength; i++)
    {
        res[i] = v1[i] / coef;
    }
}
//---
void PARALLEL::divideByScalar(prec *&v1, prec &coef, int &N, prec *&res)
{
    if (coef == 0)
        throw_line("ERROR: Dividing by zero\n");
    omp_set_num_threads(nThread);
#pragma omp parallel for
    for (int i = 0; i < N; i++)
    {
        res[i] = v1[i] / coef;
    }
}

prec PARALLEL::dot(VECTOR &v1, VECTOR &v2)
{
    if (v1.length != v2.length)
        throw_line("ERROR: Scalar product with vectors of different length\n");
    int vLength = v1.length;
    prec res = 0;
    omp_set_num_threads(nThread);
#pragma omp parallel for reduction(+ \
                                   : res)
    for (int i = 0; i < vLength; i++)
    {
        res += v1[i] * v2[i];
    }
    return res;
}
//---
prec PARALLEL::dot(prec *&v1, prec *&v2, int &N)
{
    prec res = 0;
    omp_set_num_threads(nThread);
#pragma omp parallel for reduction(+ \
                                   : res)
    for (int i = 0; i < N; i++)
    {
        res += v1[i] * v2[i];
    }
    return res;
}

prec PARALLEL::norm(VECTOR &vec)
{
    prec res = dot(vec, vec);
    res = sqrt(res);
    return res;
}
//---
prec PARALLEL::norm(prec *&vec, int &N)
{
    prec res;
    res = dot(vec, vec, N);
    res = sqrt(res);
    return res;
}

void PARALLEL::resetZeros(VECTOR &vec)
{
    int length = vec.length;
    omp_set_num_threads(nThread);
#pragma omp parallel for schedule(auto)
    for (int i = 0; i < length; i++)
        vec[i] = 0;
}
//---
void PARALLEL::resetZeros(prec *&vec, int &N)
{
    omp_set_num_threads(nThread);
#pragma omp parallel for schedule(auto)
    for (int i = 0; i < N; i++)
        vec[i] = 0;
}

void PARALLEL::copy(VECTOR &vecIn, VECTOR &vecOut)
{
    int vLength = vecIn.length;
    if (vecOut.P == 0)
        vecOut.initialize(vLength);
    omp_set_num_threads(nThread);
#pragma omp parallel for
    for (int i = 0; i < vLength; i++)
    {
        vecOut[i] = vecIn[i];
    }
}
//---
void PARALLEL::copy(prec *&vecIn, int &N, prec *&vecOut)
{
    omp_set_num_threads(nThread);
#pragma omp parallel for
    for (int i = 0; i < N; i++)
    {
        vecOut[i] = vecIn[i];
    }
}
//---
void PARALLEL::copy(int *&vecIn, int &N, int *&vecOut)
{
    omp_set_num_threads(nThread);
    #pragma omp parallel for
    for (int i = 0; i < N; i++)
    {
        vecOut[i] = vecIn[i];
    }
}

//-----
// ALLOCATE MATRIX BEFORE
void PARALLEL::prod(CSRMAT &mat, prec &coef, CSRMAT &res)
{
    int nTerm = mat.nTerm;
    PARALLEL::copyPatt(mat, res);
    std::shared_ptr<prec[]> coefRes = res.coef;
    std::shared_ptr<prec[]> coefM = mat.coef;
    omp_set_num_threads(nThread);
#pragma omp parallel for
    for (int i = 0; i < nTerm; i++)
    {
        coefRes[i] = coefM[i] * coef;
    }
}

void PARALLEL::prod(CSRMAT &mat, VECTOR &vec, VECTOR &res)
{
    int nCol = mat.nCol;
    int nRow = mat.nRow;
    std::shared_ptr<int[]> iat = mat.iat;
    std::shared_ptr<int[]> ja = mat.ja;
    std::shared_ptr<prec[]> coef = mat.coef;
    if (vec.length != nCol)
    {
        throw_line("\nERROR, wrong matrix by product dimensions\n");
    }
    std::shared_ptr<prec[]> vecP = vec.P;
    // local variables
    if (res.P == 0) res.initialize(nRow);

    //--------------------
    #pragma omp parallel for num_threads(nThread) schedule(auto)
    for (int k = 0; k < nRow; k++)
    {
        int beginRow = iat[k];
        int beginNextRow = iat[k + 1];
        prec tempRes = 0.0;
        for (int i = beginRow; i < beginNextRow; i++)
        {
            int tempCol = ja[i];
            tempRes += coef[i] * vecP[tempCol];
        }
        res[k] = tempRes;
    }
}
//---
void PARALLEL::prod(CSRMAT &mat, prec *&vec, int &N, prec *&res)
{
    int nCol = mat.nCol;
    int nRow = mat.nRow;
    std::shared_ptr<int[]> iat = mat.iat;
    std::shared_ptr<int[]> ja = mat.ja;
    std::shared_ptr<prec[]> coef = mat.coef;
    if (N != nCol)
    {
        throw_line("\nERROR, wrong matrix by vector product dimensions\n");
    }

    // local variables
    int beginRow, beginNextRow;

    //--------------------
    omp_set_num_threads(nThread);
    #pragma omp parallel for private(beginRow, beginNextRow)
    for (int k = 0; k < nRow; k++)
    {
        beginRow = iat[k];
        beginNextRow = iat[k + 1];
        res[k] = 0.0;
        for (int i = beginRow; i < beginNextRow; i++)
        {
            int tempCol = ja[i];
            res[k] = res[k] + coef[i] * vec[tempCol];
        }
    }
}

void PARALLEL::prod(CSRMAT &mat1, CSRMAT &mat2, CSRMAT &res)
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
    int incumbMaxStor = nTerm1 / nThread + nThread;

    incumbMaxStor *= 3;
    
    std::shared_ptr<int[]> iatRes(new int[nRow1 + 1]);
    VECTOR_INT jaRes(incumbMaxStor);
    VECTOR coefRes(incumbMaxStor);

    // int pos = 0;
    // int posRes = 0;
    iatRes[0] = 0;

    std::vector<VECTOR_INT> jaResThread(nThread);
    std::vector<VECTOR> coefResThread(nThread);
    std::vector<int> threadFromRow(nRow1);
    std::vector<int> posRowInThread(nRow1);
    int rangenRow = int(nRow1 / nThread) + nThread;
    std::vector<VECTOR_INT> falseIat(nThread);
    VECTOR nnzInRow(nRow1);

    VECTOR_INT myCountPos(nThread);
    PARALLEL::resetZeros(myCountPos);

    int myCountRow = 0;
    // int myCountPos = 0;
    #pragma omp parallel for num_threads(nThread)
    for (int thread = 0; thread < nThread; thread++)
    {
        jaResThread[thread].initialize(incumbMaxStor);
        coefResThread[thread].initialize(incumbMaxStor);
        falseIat[thread].initialize(rangenRow);

        falseIat[thread][0] = 0;
    }

    #pragma omp parallel for num_threads(nThread) firstprivate(myCountRow, incumbMaxStor) 
    for (int irow = 0; irow < nRow1; irow++)
    {
        int myId = omp_get_thread_num();
        int myPos = myCountPos[myId];
        int oldCountPos = myPos;
        threadFromRow[irow] = myId;
        posRowInThread[irow] = myCountRow;
        myCountRow++;

        int nElInRow = iat1[irow + 1] - iat1[irow];
        VECTOR_INT tempColInRow(nElInRow);
        VECTOR tempCoefInRow(nElInRow);
        int pos = iat1[irow];
        for (int el = 0; el < nElInRow; el++)
        {
            tempColInRow[el] = ja1[pos];
            tempCoefInRow[el] = coef1[pos];
            pos++;
        }
        VECTOR_INT colInResRow(nElInRow * nCol2);
        VECTOR coefInResRow(nElInRow * nCol2);
        int tempCount = 0;
        //------------------------------------
        for (int iter = 0; iter < nElInRow; iter++)
        {
            int jrow = tempColInRow[iter];
            for (int posIn = iat2[jrow]; posIn < iat2[jrow + 1]; posIn++)
            {
                int tempColIn = ja2[posIn];
                colInResRow[tempCount] = tempColIn;
                coefInResRow[tempCount] = tempCoefInRow[iter] * coef2[posIn];
                tempCount++;
            }
        }

        int repitedCol = CSRMAT::sort(colInResRow.P, coefInResRow.P, tempCount);


        int tempNnz = tempCount - repitedCol;
        for (int tempPos = 0; tempPos < tempNnz; tempPos++)
        {
            jaResThread[myId][myPos] = colInResRow[tempPos];
            coefResThread[myId][myPos] = coefInResRow[tempPos];
            myPos++;
            if (myPos == incumbMaxStor)
            {
                int maxRow = (int) nRow1/nThread; maxRow += nThread;
                prec perc = (prec) myCountRow/maxRow;
                incumbMaxStor += 3*incumbMaxStor*(1-perc)/perc;

                jaResThread[myId].enlarge(incumbMaxStor);
                coefResThread[myId].enlarge(incumbMaxStor);
            }
        }

        int *tempP = &(falseIat[myId][0]);
        tempP[myCountRow] = myPos;

        nnzInRow[irow] = myPos - oldCountPos;
        myCountPos[myId] = myPos;
    }

    // REASSEMBLE COMPLETE MATRIX
    int count = 0;
    int myPos = myCountPos.sum();
    jaRes.initialize(myPos);
    coefRes.initialize(myPos);
    for (int i = 0; i < nRow1; i++)
    {
        count += nnzInRow[i];
        iatRes[i + 1] = count;
    }

    #pragma omp parallel for num_threads(nThread)
    for (int irow = 0; irow < nRow1; irow++)
    {
        int posRes = iatRes[irow];
        int threadId = threadFromRow[irow];
        int posRow = posRowInThread[irow];
        int initPos = falseIat[threadId][posRow];
        int endPos = falseIat[threadId][posRow + 1];
        for (int jpos = initPos; jpos < endPos; jpos++)
        {
            jaRes[posRes] = jaResThread[threadId][jpos];
            coefRes[posRes] = coefResThread[threadId][jpos];
            posRes++;
        }
    }
    res.initialize(nRow1, nCol2, iatRes, jaRes.P, coefRes.P, myPos);
}

//--------------------------------

void PARALLEL::diag(CSRMAT &mat, VECTOR &res)
{
    int nRow = mat.nRow;
    std::shared_ptr<int[]> iat = mat.iat;
    std::shared_ptr<int[]> ja = mat.ja;
    std::shared_ptr<prec[]> coef = mat.coef;
    if (res.P == 0)
        res.initialize(nRow);
    PARALLEL::resetZeros(res);

    omp_set_num_threads(nThread);
#pragma omp parallel for schedule(auto)
    for (int i = 0; i < nRow; i++)
    {
        for (int j = iat[i]; j < iat[i + 1]; j++)
        {
            if (ja[j] == i)
            {
                res[i] = coef[j];
                break;
            }
        }
    }
}

void PARALLEL::RprodByDiag(CSRMAT &mat, VECTOR &diag, CSRMAT &res)
{
    int nRow = mat.nRow;
    int nCol = mat.nCol;
    // int nTerm = mat.nTerm;
    std::shared_ptr<int[]> iat = mat.iat;
    std::shared_ptr<int[]> ja = mat.ja;
    std::shared_ptr<prec[]> coef = mat.coef;
    if (nCol != diag.length)
        throw_line("ERROR: wrong dimensions in matrix by diag\n.");
    PARALLEL::copyPatt(mat, res);

    std::shared_ptr<prec[]> coefOut = res.coef;
    //---
    omp_set_num_threads(nThread);
#pragma omp parallel for schedule(auto)
    for (int irow = 0; irow < nRow; irow++)
    {
        for (int pos = iat[irow]; pos < iat[irow + 1]; pos++)
        {
            int tempCol = ja[pos];
            coefOut[pos] = coef[pos] * diag[tempCol];
        }
    }
}

void PARALLEL::LprodByDiag(CSRMAT &mat, VECTOR &diag, CSRMAT &res)
{
    int nRow = mat.nRow;
    // int nCol = mat.nCol;
    // int nTerm = mat.nTerm;
    std::shared_ptr<int[]> iat = mat.iat;
    std::shared_ptr<int[]> ja = mat.ja;
    std::shared_ptr<prec[]> coef = mat.coef;
    //
    if (nRow != diag.length)
        throw_line("ERROR: wrong dimensions in matrix by diag\n.");
    PARALLEL::copyPatt(mat, res);

    std::shared_ptr<prec[]> coefOut = res.coef;
    //---
    omp_set_num_threads(nThread);
#pragma omp parallel for schedule(auto)
    for (int irow = 0; irow < nRow; irow++)
    {
        for (int pos = iat[irow]; pos < iat[irow + 1]; pos++)
        {
            coefOut[pos] = coef[pos] * diag[irow];
        }
    }
}

void PARALLEL::getTranspose(CSRMAT &mat, CSRMAT &res)
{
    int nRow = mat.nRow;
    int nCol = mat.nCol;
    int nTerm = mat.nTerm;

    res.initialize(nCol, nRow, nTerm);
    std::shared_ptr<int[]> iat = mat.iat;
    std::shared_ptr<int[]> ja = mat.ja;
    std::shared_ptr<prec[]> coef = mat.coef;

    VECTOR_INT colCounter(nCol);
    PARALLEL::resetZeros(colCounter);

    std::shared_ptr<int[]> newIat = res.iat;
    std::shared_ptr<int[]> newJa = res.ja;
    std::shared_ptr<prec[]> newCoef = res.coef;

    for (int pos = 0; pos < nTerm; pos++)
    {
        int tempCol = ja[pos];
        colCounter[tempCol]++;
    }

    //------------------------
    std::shared_ptr<int *[]> newJaPP = VECTOR_INT::makeDoublePointer(nCol);
    std::shared_ptr<prec *[]> newCoefPP = VECTOR::makeDoublePointer(nCol);
    // build new iat
    newIat[0] = 0;
    int count = 0;
    for (int icol = 0; icol < nCol; icol++)
    {
        newJaPP[icol] = &newJa[count];
        newCoefPP[icol] = &newCoef[count];
        count += colCounter[icol];
        newIat[icol + 1] = count;
    }
    //---
    PARALLEL::resetZeros(colCounter);
    // build new ja and new coef2
    for (int irow = 0; irow < nRow; irow++)
    {
#pragma omp parallel for num_threads(nThread)
        for (int pos = iat[irow]; pos < iat[irow + 1]; pos++)
        {
            int tempCol = ja[pos];
            int tempCount = colCounter[tempCol];
            newJaPP[tempCol][tempCount] = irow;
            newCoefPP[tempCol][tempCount] = coef[pos];
            colCounter[tempCol]++;
        }
    }
}

void PARALLEL::copyPatt(CSRMAT &matIn, CSRMAT &matOut)
{
    int nRow = matIn.nRow;
    int nCol = matIn.nCol;
    int nTerm = matIn.nTerm;
    VECTOR_INT iat;
    iat.length = nRow + 1;
    iat.P = matIn.iat;
    VECTOR_INT ja;
    ja.length = nTerm;
    ja.P = matIn.ja;

    VECTOR_INT iatRes(nRow + 1);
    VECTOR_INT jaRes(nTerm);
    PARALLEL::copy(iat, iatRes);
    PARALLEL::copy(ja, jaRes);

    std::shared_ptr<int[]> iatResP = iatRes.P;
    std::shared_ptr<int[]> jaResP = jaRes.P;
    matOut.initialize(nRow, nCol, nTerm, iatResP, jaResP);
}

void PARALLEL::concatInRow(CSRMAT &mat1, CSRMAT &mat2, CSRMAT &resMat)
{
    int nRow1 = mat1.nRow;
    int nCol1 = mat1.nCol;
    std::shared_ptr<int []> iat1  = mat1.iat;
    std::shared_ptr<int []> ja1   = mat1.ja;
    std::shared_ptr<prec[]> coef1 = mat1.coef;

    int nRow2 = mat2.nRow;
    int nCol2 = mat2.nCol;
    std::shared_ptr<int []> iat2  = mat2.iat;
    std::shared_ptr<int []> ja2   = mat2.ja;
    std::shared_ptr<prec[]> coef2 = mat2.coef;
    if (nRow1 != nRow2) throw_line("ERROR: Different nRow in row concatenation of CSR matrices.");

    int nRowRes = nRow1;
    int nColRes = nCol1+nCol2;
    int nTermRes = mat1.nTerm + mat2.nTerm;
    std::shared_ptr<int []> iatRes(new int[nRowRes+1]);
    std::shared_ptr<int []> jaRes(new int[nTermRes]);
    std::shared_ptr<prec[]> coefRes(new prec[nTermRes]);
    iatRes[0] = 0;
    int posRes = 0;
    for (int irow = 0; irow < nRow1; irow++)
    {
        for (int pos1 = iat1[irow]; pos1 < iat1[irow+1]; pos1++)
        {
            jaRes[posRes] = ja1[pos1]; 
            coefRes[posRes] = coef1[pos1];
            posRes++;
        }
        //
        for (int pos2 = iat2[irow]; pos2 < iat2[irow+1]; pos2++)
        {
            jaRes[posRes] = ja2[pos2] + nCol1; 
            coefRes[posRes] = coef2[pos2];
            posRes++;
        }
        iatRes[irow+1] = posRes;
    }

    resMat.initialize(nRowRes, nColRes, iatRes, jaRes, coefRes, nTermRes);
}


void PARALLEL::concatInCol(CSRMAT &mat1, CSRMAT &mat2, CSRMAT &resMat, int delay)
{
    int nRow1 = mat1.nRow;
    int nCol1 = mat1.nCol;
    std::shared_ptr<int []> iat1  = mat1.iat;
    std::shared_ptr<int []> ja1   = mat1.ja;
    std::shared_ptr<prec[]> coef1 = mat1.coef;

    int nRow2 = mat2.nRow;
    int nCol2 = mat2.nCol;
    std::shared_ptr<int []> iat2  = mat2.iat;
    std::shared_ptr<int []> ja2   = mat2.ja;
    std::shared_ptr<prec[]> coef2 = mat2.coef;
    if (nCol1 != nCol2) throw_line("ERROR: Different nCol in col concatenation of CSR matrices.");

    int nRowRes = nRow1 + nRow2;
    int nColRes = nCol1;
    int nTermRes = mat1.nTerm + mat2.nTerm;
    std::shared_ptr<int []> iatRes(new int[nRowRes+1]);
    std::shared_ptr<int []> jaRes(new int[nTermRes]);
    std::shared_ptr<prec[]> coefRes(new prec[nTermRes]);
    iatRes[0] = 0;
    int posRes = 0;
    for (int irow = 0; irow < nRow1; irow++)
    {
        for (int pos1 = iat1[irow]; pos1 < iat1[irow+1]; pos1++)
        {
            jaRes[posRes] = ja1[pos1]; 
            coefRes[posRes] = coef1[pos1];
            posRes++;
        }
        iatRes[irow+1] = posRes;
    }
    //
    for (int irow = 0; irow < nRow2; irow++)
    {
        for (int pos2 = iat2[irow]; pos2 < iat2[irow+1]; pos2++)
        {
            jaRes[posRes] = ja2[pos2]+delay; 
            coefRes[posRes] = coef1[pos2];
            posRes++;
        }
        iatRes[nRow1+irow+1] = posRes;
    }

    resMat.initialize(nRowRes, nColRes, iatRes, jaRes, coefRes, nTermRes);
}
//---------------------------
// BLOCK MATRIX
//---------------------------
void PARALLEL::createBlockMatrix(CSRMAT &refMat, std::vector<std::vector<VECTOR*>> &coef, CSRMAT &matRes)
{
    std::shared_ptr<prec[]> coefH = refMat.coef;
    int dim = coef.size();
    int nRow = refMat.nRow;
    int nCol = refMat.nCol;
    int nTerm = refMat.nTerm;
    std::shared_ptr<int[]> iat = refMat.iat;
    std::shared_ptr<int[]> ja = refMat.ja;

    int nRowRes = nRow*dim; 
    int nColRes = nCol*dim;
    int nTermRes = nTerm*dim*dim;
    
    std::shared_ptr<int[]> iatRes(new int[nRowRes+1]);
    std::shared_ptr<int[]> jaRes(new int[nTermRes]);
    std::shared_ptr<prec[]> coefRes(new prec[nTermRes]);

    iatRes[0] = 0;


    #pragma omp parallel for num_threads(nThread)
    for (int icomp = 0; icomp < dim; icomp++)
    {
        int posGlob = nTerm*(icomp*dim);
        int globRow = icomp*nRow;
        for (int irow = 0; irow < nRow; irow++)
        {
            int startPos = iat[irow]; int endPos = iat[irow+1];
            //----------------------------------------
            for (int jcomp = 0; jcomp < dim; jcomp++)
            {
                VECTOR tempCoef = *(coef[icomp][jcomp]);
                if (tempCoef.length != nTerm) throw_line("ERROR: wrong dimensions in BLOCK MATRIX creation.");
                //-----------
                for (int pos = startPos; pos < endPos; pos++)
                {
                    jaRes[posGlob] = ja[pos] + nCol*jcomp;          
                    coefRes[posGlob] = tempCoef[pos];
                    posGlob++;
                }
                //------------
            }
            iatRes[globRow+1] = posGlob;
            globRow++;
        }
    }
    matRes.initialize(nRowRes, nColRes, iatRes, jaRes, coefRes, nTermRes);
}

//--------------------------------------------------------
// SOLVER GMRES
//--------------------------------------------------------
void PARALLEL::gmres_p(CSRMAT &mat, VECTOR &b, VECTOR &sol, VECTOR &x_0, prec &bNorm, prec toll, int p, int maxIt)
{
    int N = mat.nRow;

    if (bNorm > 1e12) bNorm = 1;

    // check
    if (maxIt > N)
        maxIt = N;
    if (p > N)
        p = N;
    if (N != mat.nCol)
        throw_line("ERROR, in matrix dimension, non square matrix\n");

    PARALLEL::copy(x_0, sol);
    //---------------------------------------------------
    // alloc of everything useful
    //---------------------------------------------------
    VECTOR w(N);       // new vector linear independent
    VECTOR d(p + 1);   // Q matrix (in QR decomposition) first row
    VECTOR hv_j(N);    // store temporary product h_{i,j}*V[j];
    VECTOR tempVec(N); // store temporary vector result
    VECTOR z(p);       // store results of Rz = d

    prec *hv_jP = &hv_j[0];
    prec *tempVecP = &tempVec[0];
    prec *dP = &d[0];
    // prec* zP = &z[0];

    HESSEMBERG_BY_COL H(p + 1);
    MATRIX_BY_COL V(N, p + 1);

    // alloc for QR factorization
    UPPER_TRIANGULAR R(p + 1);
    HESSEMBERG_BY_COL Q(p + 1);

    //----------------------------------------
    // INITIALIZE PARAMETERS
    //----------------------------------------
    VECTOR lastV; // pointer to last vector of V

    VECTOR Ax_old(N);
    PARALLEL::prod(mat, x_0, Ax_old);
    VECTOR r_old(N);
    PARALLEL::diff(b, Ax_old, r_old);

    prec *r_oldP = &(r_old[0]);

    prec beta = PARALLEL::norm(r_old);

    int it = 0;
    int k = 0;
    prec res;
    //------------------------------
    // OUTER LOOP
    //------------------------------
    while (beta / bNorm > toll && it < maxIt)
    {
        prec *colAddress = V[0];
        PARALLEL::divideByScalar(r_oldP, beta, N, colAddress);
        maxIt += k;
        k = 0;
        res = toll + 1;
        //---------------------------------------
        // INNER LOOP
        //--------------------------------------
        while (res > toll && k < p && it < maxIt)
        {
            prec *wP = &(w.P[0]);

            prec *VkP = V[k];
            PARALLEL::prod(mat, VkP, N, wP);
            for (int j = 0; j < k + 1; j++)
            {
                prec *tempP = V[j];
                prec h = PARALLEL::dot(wP, tempP, N);
                PARALLEL::multiplyByScalar(tempP, h, N, hv_jP);
                PARALLEL::diff(wP, hv_jP, N, wP);
                H[k][j] = h;
            }
            prec hh = PARALLEL::norm(wP, N); // h_{k+1, k}
            prec *tempP = V[k + 1];
            PARALLEL::divideByScalar(wP, hh, N, tempP); // store V last vector
            H[k][k + 1] = hh;

            res = PARALLEL::QRforGMRES(H, R, Q, dP, k + 1);

            //-----
            if (res < 0)
            {
                printf("it: %d res: %" format_e "  WARNING: NEGATIVE RESIDUAL, PROBABLY DUE TO MALCONDITIONATE MATRIX.\nRESTARTED GMRES APPLIED.\n", it, res);
                res = toll + 1;
                break;
            }
            //-----
            res = beta * res / bNorm;
            k++;
            it++;
        }
        //-----------------------------------------------
        PARALLEL::multiplyByScalar(dP, beta, k, tempVecP);

        // R.solveLS(tempVec.P, z.P, k); // find z sol
        PARALLEL::solveLS(R, tempVec, z, k);

        // compute solution
        int opXThread = k / nThread;
        #pragma omp parallel num_threads(nThread) reduction(vsum : sol)
        {
            VECTOR tempVec2(N);
            int myId = omp_get_thread_num();
            for (int i = myId * opXThread; i < (myId + 1) * opXThread; i++)
            {
                prec *tempP = V[i];
                VECTOR::multiplyByScalar(tempP, z[i], N, tempVec2.P);
                sol += tempVec2;
            }
        }

        for (int i = nThread * opXThread; i < k; i++)
        {
            prec *tempP = V[i];
            PARALLEL::multiplyByScalar(tempP, z[i], N, tempVecP);
            PARALLEL::sum(sol, tempVec, sol);
        }

        //
        VECTOR Asol;
        PARALLEL::prod(mat, sol, Asol);

        PARALLEL::diff(b, Asol, r_old);
        beta = PARALLEL::norm(r_old); // actual residuals
        if (beta / bNorm > toll)
        PARALLEL::copy(sol, x_0);
    }
    // //-------------------------------------------
    // printf("final it: %d relative res: %" format_e "\n", it, beta/bNorm);
}

prec PARALLEL::QRforGMRES(HESSEMBERG_BY_COL &mat, UPPER_TRIANGULAR &R, HESSEMBERG_BY_COL &Q, prec *&d, int nCols)
{
    //--------------------------------------------------
    int k = nCols - 1;
    int maxRowID = 2 + k;

    //-----------------------------
    // alloc something
    //-----------------------------
    VECTOR hv_j(maxRowID);
    prec *hv_jP = &hv_j[0];
    VECTOR w(maxRowID);
    prec *wP = &w[0];

    prec *tempRawP = mat[k];
    VECTOR tempP(maxRowID);
    tempP = tempRawP;
    prec *tempPP = &(tempP[0]);
    PARALLEL::copy(tempPP, maxRowID, wP);
    for (int j = 0; j < k; j++)
    {
        int locMax = 2 + j;
        prec *tempPointer = Q[j];
        VECTOR pointColj(locMax);
        pointColj = tempPointer;
        prec *pointColjP = &(pointColj[0]);
        prec h = PARALLEL::dot(wP, pointColjP, locMax);
        PARALLEL::multiplyByScalar(pointColjP, h, locMax, hv_jP);
        PARALLEL::diff(wP, hv_jP, locMax, wP);
        R.el(j, k, h);
    }
    prec hh = PARALLEL::norm(wP, maxRowID);
    R.el(k, k, hh);
    PARALLEL::divideByScalar(wP, hh, maxRowID, wP);
    prec *pointColk = Q[k];

    PARALLEL::copy(wP, maxRowID, pointColk);
    // VECTOR::copy(w.P, pointColk, maxRowID);
    d[k] = pointColk[0];

    prec res = 0.0;
#pragma omp parallel for num_threads(nThread) reduction(+ : res)
    for (int i = 0; i < nCols; i++)
        res += d[i] * d[i];

    res = (1 - res);

    // res = abs(res);
    // if (res < 0 && res > -1e-16) res = std::abs(res);
    if (res < 0)
    {
        return res;
    }
    res = sqrt(res);
    return res;
}
//--------------------------
void PARALLEL::solveLS(UPPER_TRIANGULAR &mat, VECTOR &b, VECTOR &sol, int &N)
{

    for (int i = N - 1; i > 0; i--)
    {
        prec tempRes = b[i];
        tempRes /= mat.el(i, i);

        b[i - 1] -= mat.el(i - 1, i) * tempRes;

#pragma omp parallel for num_threads(nThread)
        for (int j = i - 2; j > -1; j--)
        {
            b[j] -= mat.el(j, i) * tempRes;
        }
        sol[i] = tempRes;
    }

    sol[0] = b[0] / mat.el(0, 0);
}

void PARALLEL::SIMPLEgmres_p(CSRMAT &mat, CSRMAT &A, CSRMAT &B, CSRMAT &BT, CSRMAT &S, CSRMAT &M1BT,
                                VECTOR &b, VECTOR &sol, VECTOR &x_0, prec &bNorm, prec &finalRes, prec toll, int p, int maxIt)
{
    int N = x_0.length;

    if (bNorm > 1e12) bNorm = 1;

    // check
    if (maxIt > N)
        maxIt = N;
    if (p > N)
        p = N;

    PARALLEL::copy(x_0, sol);
    //---------------------------------------------------
    // alloc of everything useful
    //---------------------------------------------------
    VECTOR w(N);       // new vector linear independent
    VECTOR d(p + 1);   // Q matrix (in QR decomposition) first row
    VECTOR hv_j(N);    // store temporary product h_{i,j}*V[j];
    VECTOR tempVec(N); // store temporary vector result
    VECTOR y(p);       // store results of Rz = d

    prec *hv_jP = &hv_j[0];
    prec *tempVecP = &tempVec[0];
    prec *dP = &d[0];
    // prec* zP = &z[0];

    HESSEMBERG_BY_COL H(p + 1);
    MATRIX_BY_COL V(N, p + 1);

    MATRIX_BY_COL Z(N, p + 1);

    // alloc for QR factorization
    UPPER_TRIANGULAR R(p + 1);
    HESSEMBERG_BY_COL Q(p + 1);

    //----------------------------------------
    // NEW MATRICES FOR SIMPLE
    //----------------------------------------
    int nDof_v = A.nRow;
    int nDof_p = B.nRow;
    VECTOR M1(nDof_v);
    PARALLEL::diag(A, M1);
    M1.invert();
    CSRMAT A_star;
    PARALLEL::RprodByDiag(A, M1, A_star);

    VECTOR vu_k(nDof_v);
    VECTOR vp_k(nDof_p);
    prec *vu_kP = &(vu_k[0]);
    prec *vp_kP = &(vp_k[0]);

    VECTOR u0(nDof_v);
    PARALLEL::resetZeros(u0);
    VECTOR p0(nDof_p);
    PARALLEL::resetZeros(p0);

    VECTOR tempz_u(nDof_v);
    VECTOR z_u(nDof_v);
    prec *z_uP = &(z_u[0]);
    VECTOR newRhsP(nDof_p);
    VECTOR z_p(nDof_p);
    prec *z_pP = &(z_p[0]);

    CSRMAT S_star;
    VECTOR M2(nDof_p);
    PARALLEL::diag(S, M2);
    M2.invert();
    PARALLEL::RprodByDiag(S, M2, S_star);

    VECTOR M1BTz_p(nDof_v);
    //----------------------------------------
    // INITIALIZE PARAMETERS
    //----------------------------------------
    VECTOR lastV; // pointer to last vector of V

    VECTOR Ax_old(N);
    PARALLEL::prod(mat, x_0, Ax_old);
    VECTOR r_old(N);
    PARALLEL::diff(b, Ax_old, r_old);

    prec *r_oldP = &(r_old[0]);

    prec beta = PARALLEL::norm(r_old);

    int it = 0;
    int k = 0;
    prec res;
    //------------------------------
    // OUTER LOOP
    //------------------------------
    while (beta / bNorm > toll && it < maxIt)
    {
        prec *colAddress = V[0];
        PARALLEL::divideByScalar(r_oldP, beta, N, colAddress);
        maxIt += k;
        k = 0;
        res = toll + 1;
        //---------------------------------------
        // INNER LOOP
        //--------------------------------------
        while (res > toll && k < p && it < maxIt)
        {
            prec *wP = &(w.P[0]);
            prec *VkP = V[k];
            prec *ZkP = Z[k];
            //*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
            // COMPUTE Z[K] (SIMPLE ALGOTIHM)
            //*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
            // STEP 1: solve A*tempz_u = vu_k
            prec *startV = &VkP[0];
            prec *startP = startV + nDof_v;
            PARALLEL::copy(startV, nDof_v, vu_kP);
            PARALLEL::copy(startP, nDof_p, vp_kP);
            prec tempbNorm = PARALLEL::norm(vu_k);
            PARALLEL::gmres_p(A_star, vu_k, tempz_u, u0, tempbNorm, 1e-3, 500, 4);
            PARALLEL::pointProd(tempz_u, M1, tempz_u);
             //*-*-*-*-*
            // STEP 2: solve S*z_p = vp_k - B*tempz_u
            //*-*-*-*-*
            PARALLEL::prod(B, tempz_u, newRhsP);
            PARALLEL::diff(vp_k, newRhsP, newRhsP);

            tempbNorm = PARALLEL::norm(newRhsP);

            PARALLEL::resetZeros(p0);
            PARALLEL::gmres_p(S_star, newRhsP, z_p, p0, tempbNorm, 1e-3, 500, 8);
            PARALLEL::pointProd(z_p, M2, z_p);
            //*-*-*-*
            // STEP3: z_u = tempz_u - BT*z_p
            //*-*-*-*
            PARALLEL::prod(M1BT, z_p, M1BTz_p);
            PARALLEL::diff(tempz_u, M1BTz_p, z_u);
            // UPDATE LASTSOL
            startV = ZkP;
            startP = startV + nDof_v;
            PARALLEL::copy(z_uP, nDof_v, startV);
            PARALLEL::copy(z_pP, nDof_p, startP);

            //*-*-*-*-*-*-*-*-*-*-*-
            PARALLEL::prod(mat, ZkP, N, wP);
            for (int j = 0; j < k + 1; j++)
            {
                prec *tempP = V[j];
                prec h = PARALLEL::dot(wP, tempP, N);
                PARALLEL::multiplyByScalar(tempP, h, N, hv_jP);
                PARALLEL::diff(wP, hv_jP, N, wP);
                H[k][j] = h;
            }
            prec hh = PARALLEL::norm(wP, N); // h_{k+1, k}
            prec *tempP = V[k + 1];
            PARALLEL::divideByScalar(wP, hh, N, tempP); // store V last vector
            H[k][k + 1] = hh;

            res = PARALLEL::QRforGMRES(H, R, Q, dP, k + 1);

            //-----
            if (res < 0)
            {
                printf("it: %d res: %" format_e "  WARNING: NEGATIVE RESIDUAL, PROBABLY DUE TO MALCONDITIONATE MATRIX.\nRESTARTED GMRES APPLIED.\n", it, res);
                res = toll + 1;
                break;
            }
            //-----
            res = beta * res / bNorm;
            k++;
            it++;
        }
        //-----------------------------------------------
        PARALLEL::multiplyByScalar(dP, beta, k, tempVecP);

        // R.solveLS(tempVec.P, y.P, k); // find z sol

        PARALLEL::solveLS(R, tempVec, y, k);

        // compute solution
        int opXThread = k / nThread;
        #pragma omp parallel num_threads(nThread) reduction(vsum : sol)
        {
            VECTOR tempVec2(N);
            int myId = omp_get_thread_num();
            for (int i = myId * opXThread; i < (myId + 1) * opXThread; i++)
            {
                prec *tempP = Z[i];
                VECTOR::multiplyByScalar(tempP, y[i], N, tempVec2.P);
                sol += tempVec2;
            }
        }

        for (int i = nThread * opXThread; i < k; i++)
        {
            prec *tempP = Z[i];
            PARALLEL::multiplyByScalar(tempP, y[i], N, tempVecP);
            PARALLEL::sum(sol, tempVec, sol);
        }

        //
        VECTOR Asol;
        PARALLEL::prod(mat, sol, Asol);

        PARALLEL::diff(b, Asol, r_old);
        beta = PARALLEL::norm(r_old); // actual residuals
        if (beta / bNorm > toll)
            PARALLEL::copy(sol, x_0);
    }
    // //-------------------------------------------
    // printf("final it: %d relative res: %" format_e "\n", it, beta / bNorm);

    finalRes = beta/bNorm;
}

//---------------------------------------------------------------------------------------------
void PARALLEL::SCHURgmres_p(CSRMAT &mat, CSRMAT &M, CSRMAT &MR, VECTOR &b, VECTOR &sol, VECTOR &x_0, prec &bNorm, prec &finalRes, prec toll, int p, int maxIt)
{
    int N = x_0.length;

    if (bNorm > 1e12) bNorm = 1;

    // check
    if (maxIt > N)
        maxIt = N;
    if (p > N)
        p = N;

    PARALLEL::copy(x_0, sol);
    //---------------------------------------------------
    // alloc of everything useful
    //---------------------------------------------------
    VECTOR w(N);       // new vector linear independent
    VECTOR d(p + 1);   // Q matrix (in QR decomposition) first row
    VECTOR hv_j(N);    // store temporary product h_{i,j}*V[j];
    VECTOR tempVec(N); // store temporary vector result
    VECTOR y(p);       // store results of Rz = d

    prec *hv_jP = &hv_j[0];
    prec *tempVecP = &tempVec[0];
    prec *dP = &d[0];
    // prec* zP = &z[0];

    HESSEMBERG_BY_COL H(p + 1);
    MATRIX_BY_COL V(N, p + 1);

    MATRIX_BY_COL Z(N, p + 1);

    // alloc for QR factorization
    UPPER_TRIANGULAR R(p + 1);
    HESSEMBERG_BY_COL Q(p + 1);

    //----------------------------------------
    // INITIALIZE PARAMETERS
    //----------------------------------------
    VECTOR lastV(N); // pointer to last vector of V
    VECTOR lastZ(N);

    VECTOR u0(N); u0.resetZeros();

    VECTOR Ax_old(N);
    PARALLEL::prod(mat, x_0, Ax_old);
    VECTOR r_old(N);
    PARALLEL::diff(b, Ax_old, r_old);

    prec *r_oldP = &(r_old[0]);

    prec beta = PARALLEL::norm(r_old);

    int it = 0;
    int k = 0;
    prec res;
    //------------------------------
    // OUTER LOOP
    //------------------------------
    while (beta / bNorm > toll && it < maxIt)
    {
        prec *colAddress = V[0];
        PARALLEL::divideByScalar(r_oldP, beta, N, colAddress);
        maxIt += k;
        k = 0;
        res = toll + 1;
        //---------------------------------------
        // INNER LOOP
        //--------------------------------------
        while (res > toll && k < p && it < maxIt)
        {
            prec *wP = &(w.P[0]);
            prec *VkP = V[k];
            prec *ZkP = Z[k];

            lastV = VkP;
            //*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
            // COMPUTE Z[K] (SCHUR ALGOTIHM)
            //*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
            // STEP 1: solve A*tempz_u = vu_k
            prec tempbNorm = PARALLEL::norm(VkP, N);

            PARALLEL::gmres_p(M, lastV, lastZ, u0, tempbNorm, 1e-3, 500, 3);

            // PARALLEL::prod(MR, lastZ, lastV);
            
            // UPDATE LASTSOL
            prec* newZP = &(lastZ[0]);
            PARALLEL::copy(newZP, N, ZkP);

            //*-*-*-*-*-*-*-*-*-*-*-
            PARALLEL::prod(mat, ZkP, N, wP);
            for (int j = 0; j < k + 1; j++)
            {
                prec *tempP = V[j];
                prec h = PARALLEL::dot(wP, tempP, N);
                PARALLEL::multiplyByScalar(tempP, h, N, hv_jP);
                PARALLEL::diff(wP, hv_jP, N, wP);
                H[k][j] = h;
            }
            prec hh = PARALLEL::norm(wP, N); // h_{k+1, k}
            prec *tempP = V[k + 1];
            PARALLEL::divideByScalar(wP, hh, N, tempP); // store V last vector
            H[k][k + 1] = hh;

            res = PARALLEL::QRforGMRES(H, R, Q, dP, k + 1);

            //-----
            if (res < 0)
            {
                printf("it: %d res: %" format_e "  WARNING: NEGATIVE RESIDUAL, PROBABLY DUE TO MALCONDITIONATE MATRIX.\nRESTARTED GMRES APPLIED.\n", it, res);
                res = toll + 1;
                break;
            }
            //-----
            res = beta * res / bNorm;
            k++;
            it++;
        }
        //-----------------------------------------------
        PARALLEL::multiplyByScalar(dP, beta, k, tempVecP);

        // R.solveLS(tempVec.P, y.P, k); // find z sol

        PARALLEL::solveLS(R, tempVec, y, k);

        // compute solution
        int opXThread = k / nThread;
        #pragma omp parallel num_threads(nThread) reduction(vsum : sol)
        {
            VECTOR tempVec2(N);
            int myId = omp_get_thread_num();
            for (int i = myId * opXThread; i < (myId + 1) * opXThread; i++)
            {
                prec *tempP = Z[i];
                VECTOR::multiplyByScalar(tempP, y[i], N, tempVec2.P);
                sol += tempVec2;
            }
        }

        for (int i = nThread * opXThread; i < k; i++)
        {
            prec *tempP = Z[i];
            PARALLEL::multiplyByScalar(tempP, y[i], N, tempVecP);
            PARALLEL::sum(sol, tempVec, sol);
        }

        //
        VECTOR Asol;
        PARALLEL::prod(mat, sol, Asol);

        PARALLEL::diff(b, Asol, r_old);
        beta = PARALLEL::norm(r_old); // actual residuals
        if (beta / bNorm > toll)
            PARALLEL::copy(sol, x_0);
    }
    // //-------------------------------------------
    printf("final it: %d relative res: %" format_e "\n", it, beta / bNorm);

    finalRes = beta/bNorm;
}

//---------------------------
// BLOCK MATRIX GENERAL CASE
//---------------------------
void PARALLEL::createBlockMatrix(std::vector<std::vector<CSRMAT*>> &refMat, CSRMAT &matRes)
{
    
    int nRowBlock = refMat.size(); int nColBlock = refMat[0].size();
    // std::cout << "\nrowB: " << nRowBlock << "\tcolB: " << nColBlock << "\n";
    int nRowRes = 0;
    int nColRes = 0;
    int nTermRes = 0;
    for (int icomp = 0; icomp < nRowBlock; icomp++)
    {
        nRowRes += (*refMat[icomp][0]).nRow;
        for (int jcomp = 0; jcomp < nColBlock; jcomp++)
        {
                if ((*refMat[icomp][jcomp]).nRow != (*refMat[icomp][0]).nRow) throw_line("\nERROR: concatenating block matrices with different nRow\n");
                nTermRes += (*refMat[icomp][jcomp]).nTerm;
        }
    }
    for (int jcomp = 0; jcomp < nColBlock; jcomp++)
    {
        nColRes += (*refMat[0][jcomp]).nCol;
        for (int icomp = 0; icomp < nRowBlock; icomp++)
        {
                if ((*refMat[icomp][jcomp]).nCol != (*refMat[0][jcomp]).nCol) throw_line("\nERROR: concatenating block matrices with different nCol\n");
        }
    }

    // std::cout << "\nrowRes: " << nRowRes << "\tcolRes: " << nColRes<< "\n";
    // std::cout << "\ntermRes: " << nTermRes << "\n";

    // std::vector<std::vector<int>> compStartPos(dim); for (int i = 0; i < dim; i++) compStartPos[i].resize(dim);
    VECTOR_INT startBlockRows(nRowBlock+1);
    int countiRows = 0;
    startBlockRows[0] = 0;
    for (int i = 0; i < nRowBlock; i++)
    {
        countiRows += (*refMat[i][0]).nRow;
        startBlockRows[i+1] = countiRows;
    } 
    // startBlockRows.printRow("startBlockRows");

    VECTOR_INT blockColsBefore(nColBlock+1);
    int countiCols = 0;
    blockColsBefore[0] = 0;
    for (int j = 0; j < nColBlock; j++)
    {
        countiCols += (*refMat[0][j]).nCol;
        blockColsBefore[j+1] = countiCols;
    } 
    // blockColsBefore.printRow("startBlockCols");
    
    std::shared_ptr<int[]> iatRes(new int[nRowRes+1]);
    std::shared_ptr<int[]> jaRes(new int[nTermRes]);
    std::shared_ptr<prec[]> coefRes(new prec[nTermRes]);

    iatRes[0] = 0;

    int posGlob = 0;
    int rowGlob = 0;
    // int colGlob = 0;
    for (int icomp = 0; icomp < nRowBlock; icomp++)
    {   
        int nRowTemp = startBlockRows[icomp+1] - startBlockRows[icomp];
        // std::cout << "\ntempRow: " << nRowTemp << "\n"; 
        for (int irow = 0; irow < nRowTemp; irow++)
        {
            for (int jcomp = 0; jcomp < nColBlock; jcomp++)
            {
                int nColBefore = blockColsBefore[jcomp];
                CSRMAT* tempResMat = refMat[icomp][jcomp];
                // (*tempResMat).printFullCSR();
                std::shared_ptr<int[]> tempIat   = (*tempResMat).iat;
                std::shared_ptr<int[]> tempJa    = (*tempResMat).ja;
                std::shared_ptr<prec[]> tempCoef = (*tempResMat).coef;
                // int nColTemp = (*tempResMat).nCol;
                // int nTermTemp = (*tempResMat).nTerm;
                int startPos = tempIat[irow]; int endPos = tempIat[irow+1];
                //----------------------------------------
                for (int pos = startPos; pos < endPos; pos++)
                {
                    jaRes[posGlob] = tempJa[pos] + nColBefore;         
                    coefRes[posGlob] = tempCoef[pos];
                    // std::cout << "rowGlob: " << rowGlob << "\tposGlob: " << posGlob << "\n";
                    posGlob++;
                }
            }
            rowGlob++;
            iatRes[rowGlob] = posGlob;
        }
    }  
    // VECTOR::print(iatRes, nRowRes+1);
    // VECTOR::print2(jaRes, coefRes, nTermRes, "jaRes", "coefRes");
    //VECTOR::print(iatRes, nRowRes+1);
    matRes.initialize(nRowRes, nColRes, iatRes, jaRes, coefRes, nTermRes);
}

// int PARALLEL::nThread = 4;
// int main()
// {
//     int nRow = 3;
//     int nCol = 3;
//     int i_sparse[6] = {0, 0, 1, 1, 2, 2};
//     int j_sparse[6] = {0, 1, 1, 2, 0, 2};
//     prec coef_sparse[6] = {1, 3, 2, 5, 7, 4};

//     std::shared_ptr<int[]> iS(i_sparse);
//     std::shared_ptr<int[]> jS(j_sparse);
//     std::shared_ptr<prec[]> coefS(coef_sparse);

//     CSRMAT mat(nRow, nCol, 6, iS, jS, coefS);

//     mat.printFullCSR();

//     int dim = 2;
//     std::vector<std::vector<VECTOR*>> prova(dim);
//     VECTOR vec1(6);
//     vec1[0] = 1; vec1[1] = 3; vec1[2] = 2; vec1[3] = 5; vec1[4] = 7; vec1[5] = 4;

//     VECTOR vec2(6);
//     vec2[0] = -4; vec2[1] = 4; vec2[2] = -4; vec2[3] = 4; vec2[4] = -4; vec2[5] = 4;

//     int nTerm = 5;
//     for (int icomp = 0; icomp < dim; icomp++)
//     {
//         prova[icomp].resize(dim);
//         for (int jcomp = 0; jcomp < dim; jcomp++)
//         {
//             if (jcomp == 0) (prova[icomp][jcomp]) = &vec1;
//             else  (prova[icomp][jcomp]) = &vec2;
//         }
//     }

//     CSRMAT res;
//     printf("non\n");
//     PARALLEL::createBlockMatrix(mat, prova, res);

//     res.printFullCSR();
//     // int nRow = 5e0;
//     // int nCol = 5e0;
//     // int num = 10e0;
//     // PARALLEL::nThread = 3;
//     // VECTOR vec1(nCol);
//     // VECTOR vec2(nCol);

//     // for (int i = 0; i < nCol; i++)
//     // {
//     //     vec1[i] = 0;
//     //     vec2[i] = 2;
//     // }

//     // #pragma omp parallel num_threads(3) reduction(vsum: vec1)
//     // {
//     //     vec1.reset(1.0);
//     //     #pragma omp critical
//     //     {
//     //     vec1.print();
//     //     }
//     // }
//     // vec1.print();

//     // CSRMAT mat(nRow, nCol, num);
//     // std::shared_ptr<int[]> iat = mat.iat;
//     // std::shared_ptr<int[]> ja = mat.ja;
//     // std::shared_ptr<prec[]> coef = mat.coef;

//     // iat[0] = 0;
//     // int pos = 0;
//     // for (int irow = 0; irow < nRow; irow++)
//     // {
//     //     for (int j = 0; j < 2; j++)
//     //     {
//     //         ja[pos] = irow+j;
//     //         if (irow == nRow-1) ja[pos] = 3+j;
//     //         coef[pos] = 1.0+j;
//     //         pos++;
//     //     }
//     //     iat[irow+1] = pos;
//     // }

//     // int nTerm = 7;
//     // VECTOR_INT irow(nTerm);
//     // irow[0] = 0; irow[1] = 0; irow[2] = 1;
//     // irow[3] = 2; irow[4] = 2; irow[5] = 3;
//     // irow[6] = 3;

//     // VECTOR_INT icol(nTerm);
//     // icol[0] = 0; icol[1] = 1; icol[2] = 1;
//     // icol[3] = 1; icol[4] = 2; icol[5] = 0; icol[6] = 3;

//     // VECTOR coef(nTerm);
//     // for (int i = 0; i < nTerm; i++) coef[i] = i+1;

//     // CSRMAT prova(4,4,nTerm, irow.P, icol.P, coef.P);

//     // CSRMAT prova2 = prova;

//     // mat.printFullCSR();

//     // for (int irow = nRow-100; irow < nRow; irow++)
//     // {
//     //     for (int j = 0; j < 100; j++)
//     //     {
//     //         ja[pos] = j+irow-100;
//     //         coef[pos] = 1.0;
//     //         pos++;
//     //     }
//     //     iat[irow+1] = pos;
//     // }

//     // if (pos != num) throw_line("SOMETHING WENT WRONG\n");

//     // mat.writePatt("path.txt");
//     // Record start time
//     // auto startTime = std::chrono::high_resolution_clock::now();

//     // VECTOR vec3 = PARALLEL::copy(vec1);
//     // CSRMAT res = PARALLEL::copyPatt(mat);

//     // prova.printFullCSR();

//     // CSRMAT res = PARALLEL::prod(mat, mat);

//     // // CSRMAT res = prova*prova;

//     // // CSRMAT res = mat.getTranspose();
//     // // Record end time
//     // auto finishTime = std::chrono::high_resolution_clock::now();

//     // // res.printFullCSR();
//     // std::chrono::duration<double> deltaT = finishTime - startTime;
//     // std::cout << "TIME: " << deltaT.count() << "\n";

//     // res.printFullCSR();
//     // res.writePatt("pathRes.txt");

//     // prec test = mat.nTerm;//(8,9);

//     // printf("valore: %f\n", test);
//     // printf("res: %" format "\n", res);

//     // res.print();
//     // for (int i = 0; i < nRow; i++)
//     // {
//     //     // if (res[i] != 1.0) throw_line("SOMETHING WENT WRONG\n");
//     // }

//     return 1;
// }