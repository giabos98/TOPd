#pragma once

// #include "../../CODE_HEADERS/codeHeader.h"

class INSERTION_SORT
{
    public:

    //---------------------------------------------------------------------------------
    // INSERTION SORT ALGORITHM: sort columns and coeffs w.r.t. columns for each row
    //---------------------------------------------------------------------------------
    static std::shared_ptr<int[]> insertionSortWithCoeff(std::shared_ptr<int[]> &colInRow, int N)
    {
        std::shared_ptr<int[]> reorder(new int[N]);
        for (int i = 0; i < N; i++)
        {
            reorder[i] = i;
        }
        int j;
        for (int i = 1; i < N; i++)
        {
            j = i;
            while (j > 0 && colInRow[j-1] > colInRow[j])
            {
                swap(colInRow, j-1, j);
                j = j - 1;
            }
            moveTo(reorder, i, j);
        } 
        return reorder;
    }
    
    //---------------------------------------------------------  
    // INSERTION SORT ALGORITHM
    //---------------------------------------------------------  
    static void insertionSort(std::shared_ptr<int[]> &A, int N)
    {
        int j;
        for (int i = 1; i < N; i++)
        {
            j = i;
            while (j > 0 && A[j-1] > A[j])
            {
                swap(A, j-1, j);
                j = j - 1;
            }
        } 
    }
    //-------------------------------------------
    // SWAP THE VALUES OF 2 INDICIES IN A POINTER
    //-------------------------------------------
    static void swap(std::shared_ptr<int[]> &P, int id1, int id2)
    {
        if (id1 == id2) throw_line("ERROR: Trying to swap or insert using the same index");
        //if (id1 > id2) swapVal(id1, id2);  
        int temp = P[id2];
        P[id2] = P[id1];
        P[id1] = temp;
    }

    //--------------------------------------------------------------
    // MOVE A VALUE IN A POINTER INTO ANOTHER ONE SCALING THE OTHERS
    //--------------------------------------------------------------
    static void moveTo(std::shared_ptr<int[]> &reorder, int oldID, int newID)
    {
        if (oldID != newID)
        {
            int temp = reorder[oldID];
            for (int i = oldID; i > newID; i--)
            {
                reorder[i] = reorder[i - 1];
            }   
            reorder[newID] =  temp;
        }
    }

    //---------------------------------
    // SWAP 2 VALUES
    //---------------------------------
    static void swapVal(int &a, int &b)
    {
        int temp = b;
        b = a;
        a = temp;
    }
    //---
    static void swapVal(prec &a, prec &b)
    {
        prec temp = b;
        b = a;
        a = temp;
    }

};