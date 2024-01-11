#include "../../CODE_HEADERS/codeHeader.h"
// #include "CSRMat.h"

// CSRMAT getBFromSYSMAT(CSRMAT SYS)
//     {
//         int nTerm = 8;
//         CSRMAT B(3, 4, nTerm);
//         std::shared_ptr<int[]> iatB = B.iat; std::shared_ptr<int[]> jaB = B.ja; std::shared_ptr<prec[]> coefB = B.coef;

//         std::shared_ptr<int[]> iat = SYS.iat; std::shared_ptr<int[]> ja = SYS.ja; std::shared_ptr<prec[]> coef = SYS.coef;
//         int lastRow = 7;

//         int pos = 0;
//         iatB[0] = 0;
//         int irowB = 0;
//         for (int irow = 4; irow < lastRow; irow++)
//         {
//             for (int tempPos = iat[irow]; tempPos < iat[irow+1]; tempPos++)
//             {
//                 jaB[pos] = ja[tempPos];
//                 coefB[pos]  = coef[tempPos];
//                 pos++;
//             }
//             iatB[irowB+1] = pos;
//             irowB++;
//         }
//         return B;
//     } 

int PARALLEL::nThread = 4;

int main()
{
    srand((unsigned) time(NULL));

    int nRow = 4;
    int nCol = 4;
    int nTerm = 3*3;

    int nMat = 2;
    std::vector<CSRMAT> SYS(nMat);

    for (int imat = 0; imat < nMat; imat++)
    {
        VECTOR_INT iSparse(nTerm);
        VECTOR_INT jSparse(nTerm);
        VECTOR coef(nTerm);
        for (int i = 0; i < nTerm; i++)
        {
            iSparse[i] = (rand() % nRow);
        }
        for (int i = 0; i < nTerm; i++)
        {
            jSparse[i] = (rand() % nCol);
        }
        for (int i = 0; i < nTerm; i++)
        {
            int randCoef = (rand() % 100);
            coef[i] = randCoef / (double)10.0;
        }

        // VECTOR checkRow(nRow); checkRow.resetZeros();
        // for (int i = 0; i < nTerm; i++)
        // {
        //     checkRow[iSparse[i]] = 1;
        // }
        // for (int i = 0; i < nRow; i++)
        // {
        //     if (checkRow[i] != 1) std::cout << "\nrow to zero: " << i << "\n";
        // }

        
        // // VECTOR_INT::print2(iSparse.P, jSparse.P, nTerm, "iSparse", "jSparse");
        // // coef.printRow("coef");


        // std::shared_ptr<int[]> isP(new int[nTerm]); 
        // std::shared_ptr<int[]> jsP(new int[nTerm]);
        // std::shared_ptr<prec[]> coefsP(new prec[nTerm]);

        // for (int i = 0; i < nTerm; i++)
        // {
        //     isP[i] = iSparse[i];
        //     jsP[i] = jSparse[i];
        //     coefsP[i] = coef[i];
        // }

        SYS[imat].initialize(nRow, nCol, nTerm, iSparse.P, jSparse.P, coef.P);
        // std::string tempName = "MATLAB_TESTING/mat" + std::to_string(imat+1) + ".txt"; const char* tempNameC = tempName.c_str();
        // SYS[imat].write(tempNameC);

        VECTOR::print2(iSparse.P, jSparse.P, nTerm, "IGLOB", "JGLOB");
        coef.print();


        VECTOR_INT::print(SYS[imat].iat, nRow+1);
        VECTOR::print2(SYS[imat].ja, SYS[imat].coef, SYS[imat].nTerm, "ja", "coef");
        // VECTOR::print(SYS[0].coef, SYS[0].nTerm);
        SYS[imat].printFullCSR(); 
    }

    CSRMAT SUM = SYS[0] + SYS[1];

    SYS[0].printFullCSR(); 
    SYS[1].printFullCSR();
    SUM.printFullCSR();
}