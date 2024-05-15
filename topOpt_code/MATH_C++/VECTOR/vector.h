#pragma once

#include "vector_int.h"
#include "../../CODE_HEADERS/codeHeader.h"

class VECTOR
{
public:
    std::shared_ptr<prec[]> P = 0;
    int length = 0;
    static int nThread;
    //-------------------------------------------------------
    // CONSTRUCTOR
    //-------------------------------------------------------
    // EMPTY CONSTRUCTOR
    VECTOR() 
    {
        P = 0;
    }	
    //---
    VECTOR(int N)
    {   
        length = N;
        P = std::shared_ptr<prec[]>(new prec[N]);
    }
    //---
    VECTOR(prec* &vecPointer, int N)
    {  
        P = std::shared_ptr<prec[]>(new(vecPointer) prec[N], [](prec* vecPointer){free(vecPointer);});
        length = N;
    }
    //---
    VECTOR(std::shared_ptr<prec[]> &vecPointer, int N)
    {  
        P = vecPointer;
        length = N;
    }
    // //---
    // VECTOR(VECTOR &inVec)
    // {
    //     length = inVec.length;
    //     prec* P_in = inVec.P;
    //     P = (prec*) malloc(length * sizeof(prec));
    //     for (int i = 0; i < length; i++)
    //     {
    //         P[i] = P_in[i];
    //     }
    // }
    ~VECTOR()
    {}

    //-----------------------
    // INITIALIZE
    //-----------------------
    void initialize(int N)
    {
        length = N;
        P = std::shared_ptr<prec[]>(new prec[N]);
    }
    void initialize(prec* &vecPointer, int N)
    {  
        initialize(N);
        for (int i = 0; i < N; i++)
        {
            P[i] = vecPointer[i];
        }
        // length = N;
        // P = std::shared_ptr<prec[]>(new(vecPointer) prec[N], [](prec* vecPointer){free(vecPointer);});
    }
    //---
    static std::shared_ptr<prec[]> makePointer(int N)
    {
        std::shared_ptr<prec[]> Pout = 0;
        Pout = std::shared_ptr<prec[]>(new prec[N]);
        if (Pout == 0) throw_line("ERROR in ptr allocation, probably insufficient memory\n");

        return Pout;
    }
    //---
    static std::shared_ptr<prec*[]> makeDoublePointer(int N)
    {
        std::shared_ptr<prec*[]> Pout = 0;
        Pout = std::shared_ptr<prec*[]>(new prec*[N]);
        if (Pout == 0) throw_line("ERROR in ptr allocation, probably insufficient memory\n");
        return Pout;
    }
    //---
    static VECTOR zeros(int N)
    {
        VECTOR vec(N);
        for (int i = 0; i < N; i++) vec[i] = 0;
        return vec;
    }
    //---
    void setZeros(int N)
    {
        P.reset(); P = 0;
        length = N;
        P = std::shared_ptr<prec[]>(new prec[N]);
        for (int i = 0; i < length; i++) P[i] = 0;
    }
    //---
    
    void define(int N, std::shared_ptr<prec[]> &inVec)
    {
        initialize(N);
        for (int i = 0; i < length; i++)
        {
            P[i] = inVec[i];
        }
    }
    //---
    void set2Values(prec val0, prec val1)
    {
        P[0] = val0; P[1] = val1;
    }

    void set3Values(prec val0, prec val1, prec val2)
    {
        P[0] = val0; P[1] = val1; P[2] = val2;
    }

    void set4Values(prec val0, prec val1, int val2, prec val3)
    {
        P[0] = val0; P[1] = val1; P[2] = val2; P[3] = val3;
    }
    //---
    static void set2Values(prec* &P, prec val0, prec val1)
    {
        P[0] = val0; P[1] = val1;
    }

    static void set3Values(prec* &P, prec val0, prec val1, prec val2)
    {
        P[0] = val0; P[1] = val1; P[2] = val2;
    }

    static void set4Values(prec* &P, prec val0, prec val1, int val2, prec val3)
    {
        P[0] = val0; P[1] = val1; P[2] = val2; P[3] = val3;
    }
    //---

    //-----------------------------------
    // RESET THE COMPS TO ZERO
    //-----------------------------------
    void resetZeros()
    {
        for (int i = 0; i < length; i++)
        {
            P[i] = 0;
        }
    }
    //---
    void reset(prec factor)
    {
        for (int i = 0; i < length; i++)
        {
            P[i] = factor;
        }
    }
    //---

    //-------------------------------------------------------
    // GET SINGLE ENTRY
    //-------------------------------------------------------
    prec& operator [] (int index)
    {
        if (index >= length || index < 0) throw_line("ERROR: index out of bounds\n");
        return P[index];
    }

    //-------------------------------------------------------
    // GET ENTRY (!!!works also with negative entries: -1=length-1, -length=0)
    //-------------------------------------------------------
    prec get(int id)
    {
        if (id >= length || id < (-1)*length)
        {
            throw_line("ERROR: index out of bounds\n");
        } 
        else if (id >= 0)
        {
            return P[id];
        }
        else 
        {
            int real_id = length + id;
            return P[real_id];
        }
    }

    //-------------------------------------------------------
    // GET LAST ENTRY
    //-------------------------------------------------------
    prec get_last()
    {
        return P[length-1];
    }

    //-------------------------------------------------------
    // DOT PRODUCT
    //-------------------------------------------------------
    prec dot(VECTOR inVec)
    {
        if (length != inVec.length) printf("ERROR: the vectors have different lengths \n;");
        prec res = 0;
        std::shared_ptr<prec[]> pIn = inVec.P;
        for (int i = 0; i < length; i++)
        {
            res += P[i] * pIn[i];
        }
        return res;
    } 
    //---
    prec dot(std::shared_ptr<prec[]> &inVec)
    {
        prec res = 0;
        for (int i = 0; i < length; i++)
        {
            res += P[i] * inVec[i];
        }
        return res;
    } 

    //-----------------------------------------------------
    // DOT PRODUCT OF POINTERS UNTIL COMPONENT N
    //-----------------------------------------------------  
    static prec dot(std::shared_ptr<prec[]> &vec1, std::shared_ptr<prec[]> &vec2, int N)
    {
        prec res = 0;
        for (int i = 0; i < N; i++)
        {
            res += vec1[i] * vec2[i];
        }
        return res;
    }
    //---
    static prec dot(std::shared_ptr<prec[]> &vec1, prec* &vec2, int N)
    {
        prec res = 0;
        for (int i = 0; i < N; i++)
        {
            res += vec1[i] * vec2[i];
        }
        return res;
    }

    //-----------------------------------------------------
    // DOT PRODUCT WITH BEGIN & END
    //-----------------------------------------------------
    prec dot(VECTOR &inVec, int begin, int end)
    {
        if (length != inVec.length) printf("ERROR: the vectors have different lengths \n;");
        if (begin < 0 || end > length) printf("ERROR: The indices are not valid \n");
        prec res = 0;
        std::shared_ptr<prec[]>  P_in = inVec.P;
        for (int i = begin; i < end; i++)
        {
            res += P[i] * P_in[i];
        }
        return res;
    } 
    
    //-----------------------------------------------------
    // POINTWISE VECTOR PRODUCT
    //-----------------------------------------------------
    VECTOR operator * (VECTOR inVec)
    {
        if (length != inVec.length) printf("ERROR: the vectors have different lengths \n;");
        VECTOR res(length);    
        std::shared_ptr<prec[]> resP = res.P;        
        std::shared_ptr<prec[]>  P_in = inVec.P;
        for (int i = 0; i < length; i++)
        {
            resP[i] = P[i] * P_in[i];
        }
        return res;
    }
    //---
    void pointdot(prec* &inVec, VECTOR &res)
    {   
        if (length != res.length) printf("ERROR: the vectors have different lengths \n;");
        for (int i = 0; i < length; i++)
        {
            res.P[i] = P[i] * inVec[i];
        }
    }
    //---
    void pointdot(std::shared_ptr<prec[]>  &inVec, std::shared_ptr<prec[]>  &res)
    {   
        if(res == 0)
        {
           res = std::shared_ptr<prec[]>(new prec[length]);
        }

        for (int i = 0; i < length; i++)
        {
            res[i] = P[i] * inVec[i];
        }
    }
    //---
    static void pointdot(std::shared_ptr<prec[]> &vec1, std::shared_ptr<prec[]> &vec2, int N, std::shared_ptr<prec[]> &res)
    {
        if(res == 0)
        {
            res = std::shared_ptr<prec[]>(new prec[N]);
        } 
        for (int i = 0; i < N; i++)
        {
            res[i] = vec1[i] * vec2[i];
        }
    }
    //---
    // POINT DIVIDE
    //---
    VECTOR pointdiv(VECTOR inVec)
    {   
        VECTOR res(length);
        if (length != inVec.length) printf("ERROR: the vectors have different lengths \n;");
        for (int i = 0; i < length; i++)
        {
            res.P[i] = P[i] / inVec[i];
        }
        return res;
    }
    //---
    void pointdiv(VECTOR &inVec, VECTOR &res)
    {   
        if (length != res.length) printf("ERROR: the vectors have different lengths \n;");
        std::shared_ptr<prec[]> Pin = inVec.P;
        for (int i = 0; i < length; i++)
        {
            res.P[i] = P[i] / Pin[i];
        }
    }
    //---
    void pointdiv(std::shared_ptr<prec[]> &inVec, VECTOR &res)
    {   
        if (length != res.length) throw_line("ERROR: the vectors have different lengths \n;");
        for (int i = 0; i < length; i++)
        {
            res.P[i] = P[i] / inVec[i];
        }
    }
    //---
    void pointdiv(std::shared_ptr<prec[]> &inVec, std::shared_ptr<prec[]> &res)
    {   
        if(res == 0)
        {
            res = std::shared_ptr<prec[]>(new prec[length]);
        } 
        for (int i = 0; i < length; i++)
        {
            res[i] = P[i] / inVec[i];
        }
    }
    //---
    static void pointdiv(std::shared_ptr<prec[]> &vec1, std::shared_ptr<prec[]> &vec2, int N, std::shared_ptr<prec[]> &res)
    {
        if(res == 0)
        {
            res = std::shared_ptr<prec[]>(new prec[N]);
        } 
        for (int i = 0; i < N; i++)
        {
            res[i] = vec1[i] / vec2[i];
        }
    }
    //----------------------------------------------------------
    // SELF OPERATIONS
    //----------------------------------------------------------
    //---
    void operator += (VECTOR inVec)
    {
        if (length != inVec.length) throw_line("ERROR: the vectors have different lengths \n;")
        std::shared_ptr<prec[]> P_in = inVec.P;
        for (int i = 0; i < length; i++)
        {
            P[i] += P_in[i];
        }
    }
    //------------------------
    void operator += (prec inScalar)
    {
        for (int i = 0; i < length; i++)
        {
            P[i] += inScalar;
        }
    }

    void operator -= (VECTOR &inVec)
    {
        if (length != inVec.length) throw_line("ERROR: the vectors have different lengths \n;")
        std::shared_ptr<prec[]> P_in = inVec.P;
        for (int i = 0; i < length; i++)
        {
            P[i] -= P_in[i];
        }
    }

    void operator /= (prec coef)
    {
        for (int i = 0; i < length; i++)
        {
            P[i] /= coef;
        }
    }

    void operator *= (prec coef)
    {
        for (int i = 0; i < length; i++)
        {
            P[i] *= coef;
        }
    }

    void operator *= (VECTOR &vec)
    {
        if (vec.length != length) throw_line("ERROR: Pointwise product of vectors of different lenght\n");
        for (int i = 0; i < length; i++)
        {
            P[i] *= vec[i];
        }
    }

    void operator /= (VECTOR vec)
    {
        if (vec.length != length) throw_line("ERROR: Pointwise division of vectors of different lenght\n");
        for (int i = 0; i < length; i++)
        {
            if (vec[i] == 0) throw_line("ERROR: Dividing by zero in pointwise divsion of vectors\n");
            P[i] /= vec[i];
        }
    }
    
    //----------------------------------------------------------
    // VECTOR SUM
    //----------------------------------------------------------
    VECTOR operator + (VECTOR inVec)
    {
        if (length != inVec.length) throw_line("ERROR: the vectors have different lengths \n;")
        VECTOR res(length);
        std::shared_ptr<prec[]> resP = res.P;
        std::shared_ptr<prec[]> P_in = inVec.P;
        for (int i = 0; i < length; i++)
        {
            resP[i] = P[i] + P_in[i];
        }
        return res;
    }
    //----------------------------------------------------------
    VECTOR operator + (prec coef)
    {
        VECTOR res(length);
        std::shared_ptr<prec[]> resP = res.P;
        for (int i = 0; i < length; i++)
        {
            resP[i] = P[i] + coef;
        }
        return res;
    }
    //---
    prec sum()
    {
        prec res = 0;
        for (int i = 0; i < length; i++) res += P[i];
        return res;
    }
    //---
    prec sum_comps()
    {
        prec res = 0;
        for (int i = 0; i < length; i++) res += P[i];
        return res;
    }
    //---
    void sum(std::shared_ptr<prec[]> &inVec, VECTOR &res)
    {
        for (int i = 0; i < length; i++)
        {
            res.P[i] = P[i] + inVec[i];
        }
    }
    //---
    void sum(std::shared_ptr<prec[]> &inVec, std::shared_ptr<prec[]> &res)
    {
        if(res == 0) 
        {
            res = std::shared_ptr<prec[]>(new prec[length]);
        }
        for (int i = 0; i < length; i++)
        {
            res[i] = P[i] + inVec[i];
        }
    }
    //---
    static void sum(std::shared_ptr<prec[]> &vec1, std::shared_ptr<prec[]> &vec2, int N, std::shared_ptr<prec[]> &res)
    {
        if(res == 0) 
        {
            res = std::shared_ptr<prec[]>(new prec[N]);
        }
        for (int i = 0; i < N; i++)
        {
            res[i] = vec1[i] + vec2[i];
        }
    }
    //---
    static void sum(VECTOR &vec1, VECTOR &vec2, VECTOR &res)
    {
        int N = vec1.length;
        if(res.length == 0) 
        {
            res.initialize(N);
        }
        for (int i = 0; i < N; i++)
        {
            res[i] = vec1[i] + vec2[i];
        }
    }
    //---
    static void sum(prec* &vec1, std::shared_ptr<prec[]> &vec2, int N, prec* &res)
    {
        if(res == 0) 
        {
            throw_line("ERROR somewhere\n");
        }
        for (int i = 0; i < N; i++)
        {
            res[i] = vec1[i] + vec2[i];
        }
    }

    //----------------------------------------------------------
    // VECTOR MEAN
    //----------------------------------------------------------
    static prec mean(VECTOR &vec)
    {
        if (vec.length == 0) throw_line("ERROR: trying to evaluate the mean of an empty vector.\n");
        if (vec.P == 0) throw_line("ERROR: trying to evaluate the mean of an empty vector.\n");
        prec tempMean = vec.sum();
        tempMean = tempMean / vec.length;
        return tempMean;
    }

    static prec weightedMean(VECTOR &vec, VECTOR &weights)
    {
        if (vec.length == 0) throw_line("ERROR: trying to evaluate the artihmetic weighted mean of an empty vector.\n");
        if (vec.P == 0) throw_line("ERROR: trying to evaluate the mean of an empty vector.\n");
        if (vec.length != weights.length) throw_line("ERROR vectors and comp. weights have different lengths\n");
        prec tempMean = 0.0;
        for (int i = 0; i < vec.length; i++)
        {
            tempMean += vec[i]*weights[i];
        }
        return tempMean;
    }

    static prec geometricMean(VECTOR &vec, VECTOR &weights)
    {
        if (vec.length == 0) throw_line("ERROR: trying to evaluate the geometric mean of an empty vector.\n");
        if (vec.P == 0) throw_line("ERROR: trying to evaluate the mean of an empty vector.\n");
        if (vec.length != weights.length) throw_line("ERROR vectors and comp. weights have different lengths\n");
        prec tempMean = 1;
        for (int i = 0; i < vec.length; i++)
        {
            tempMean *= pow(vec[i], weights[i]);
        }
        return tempMean;
    }

    static prec meanByIdSet(VECTOR &vec, VECTOR_INT &idSet)
    {
        if (idSet.length > vec.length) throw_line("ERROR in meanByIdSet: the set size is bigger than the vector one\n")
        int maxId = idSet.max();
        if (maxId > vec.length) throw_line("ERROR in meanByIdSet: the set size maxId is bigger than the vector one\n")
        VECTOR tempVec(idSet.length);
        for (int i = 0; i < tempVec.length; i++)
        {
            tempVec[i] = vec[idSet[i]];
        }
        prec tempMean = VECTOR::mean(tempVec);
        return tempMean;
    }

    static prec weigthedMeanByIdSet(VECTOR &vec, VECTOR_INT &idSet, VECTOR &weights)
    {
        if (idSet.length > vec.length) throw_line("ERROR in meanByIdSet: the set size is bigger than the vector one\n")
        if (idSet.length != weights.length) throw_line("ERROR in meanByIdSet: the idSet and weights sizes are different\n")
        int maxId = idSet.max();
        if (maxId > vec.length) throw_line("ERROR in meanByIdSet: the set size maxId is bigger than the vector one\n")
        VECTOR tempVec(idSet.length);
        for (int i = 0; i < tempVec.length; i++)
        {
            tempVec[i] = vec[idSet[i]];
        }
        prec tempMean = VECTOR::weightedMean(tempVec, weights);
        return tempMean;
    }

    static prec geometricMeanByIdSet(VECTOR &vec, VECTOR_INT &idSet, VECTOR &weights)
    {
        if (idSet.length > vec.length) throw_line("ERROR in meanByIdSet: the set size is bigger than the vector one\n")
        if (idSet.length != weights.length) throw_line("ERROR in meanByIdSet: the idSet and weights sizes are different\n")
        int maxId = idSet.max();
        if (maxId > vec.length) throw_line("ERROR in meanByIdSet: the set size maxId is bigger than the vector one\n")
        VECTOR tempVec(idSet.length);
        for (int i = 0; i < tempVec.length; i++)
        {
            tempVec[i] = vec[idSet[i]];
        }
        prec tempMean = VECTOR::geometricMean(tempVec, weights);
        return tempMean;
    }

    //----------------------------------------------------------
    // VECTOR DIFFERENCE
    //----------------------------------------------------------
    VECTOR operator - (VECTOR inVec)
    {
        if (length != inVec.length) throw_line("ERROR in VECTOR DIFFERENCE: the vectors have different lengths \n")
        VECTOR res(length);
        std::shared_ptr<prec[]> resP = res.P;
        std::shared_ptr<prec[]> P_in = inVec.P;
        for (int i = 0; i < length; i++)
        {
            resP[i] = P[i] - P_in[i];
        }
        return res;
    }
    //---
    VECTOR operator - (prec coef)
    {
        VECTOR res(length);
        std::shared_ptr<prec[]> resP = res.P;
        for (int i = 0; i < length; i++)
        {
            resP[i] = P[i] - coef;
        }
        return res;
    }
    //---
    void diff(std::shared_ptr<prec[]> &inVec, VECTOR &res)
    {
        for (int i = 0; i < length; i++)
        {
            res.P[i] = P[i] - inVec[i];
        }
    }
    //---
    void diff(std::shared_ptr<prec[]> &inVec, std::shared_ptr<prec[]> &res)
    {
        if(res == 0) 
        {
            res = std::shared_ptr<prec[]>(new prec[length]);
        }
        for (int i = 0; i < length; i++)
        {
            res[i] = P[i] - inVec[i];
        }
    }
    //---
    static void diff(std::shared_ptr<prec[]> &vec1, std::shared_ptr<prec[]> &vec2, int N, std::shared_ptr<prec[]> &res)
    {
        if(res == 0) 
        {
            res = std::shared_ptr<prec[]>(new prec[N]);
        }
        for (int i = 0; i < N; i++)
        {
            res[i] = vec1[i] - vec2[i];
        }
    }

    //---------------------------------------------------------------
    // DIVIDE BY SCALAR
    //---------------------------------------------------------------
    VECTOR operator / (prec scalar)
    {
        VECTOR res(length);
        std::shared_ptr<prec[]> resP = res.P;
        for (int i = 0; i < length; i++)
        {
            resP[i] = P[i] / scalar;
        }
        return res;
    }
    //---
    void divideByScalar(prec scalar, VECTOR &res)
    {
        for (int i = 0; i < length; i++)
        {
            res.P[i] = P[i]/scalar;
        }
    }
    //---
    void divideByScalar(prec scalar, std::shared_ptr<prec[]> &res)
    {
        if(res == 0) 
        {
            res = std::shared_ptr<prec[]>(new prec[length]);
        }      
        for (int i = 0; i < length; i++)
        {
            res[i] = P[i]/scalar;
        }
    }
    //---
    static void divideByScalar( std::shared_ptr<prec[]> &vec, prec scalar, int N,  std::shared_ptr<prec[]> &res)
    {
        if(res == 0) 
        {
            res = std::shared_ptr<prec[]>(new prec[N]);
        }
        for (int i = 0; i < N; i++)
        {
            res[i] = vec[i]/scalar;
        }
    }
    //---
    static void divideByScalar( std::shared_ptr<prec[]> &vec, prec scalar, int N,  prec* & res)
    {
        if(res == 0 || res == nullptr) 
        {
            std::shared_ptr<prec[]> resP(new prec[N]);
            res = &(resP[0]);
        }
        for (int i = 0; i < N; i++)
        {
            res[i] = vec[i]/scalar;
        }
    }

    void squared_root()
    {
        for (int i = 0; i < length; i++)
        {
            P[i] = sqrt(P[i]);
        }
    }

    void squared_root(VECTOR &res)
    {
        if (res.length != length) 
        {
            throw_line("ERROR: vectors have different legnths\n");
        }
        else
        {
            for (int i = 0; i < length; i++)
            {
                res[i] = sqrt(P[i]);
            }
        }
    }

    void power(prec power)
    {
        for (int i = 0; i < length; i++)
        {
            prec temp = P[i];
            P[i] = std::pow(temp, power);
        }
    }

    void power(int power, VECTOR &res)
    {
        if (res.length != length) 
        {
            throw_line("ERROR: vectors have different legnths\n");
        }
        else
        {
            for (int i = 0; i < length; i++)
            {
                res[i] = std::pow(P[i], power);
            }
        }
    }

    // //---------------------------------------------------
    // // MULTIPLY TRANSPOSE VECTOR BY MATRIX
    // //---------------------------------------------------
    // VECTOR operator & (MATRIX mat)
    // {
    //     int nRows = mat.nRow;
    //     if (length != nRows) throw_line("ERROR: vector and matrix sizes are not compatible \n")        
    //     int nCols = mat.nCol;
    //     VECTOR res(nCols);
    //     prec* resP = res.P;
    //     for (int i = 0; i < nCols; i++)
    //     {
    //         for (int j = 0; j < nRows; j++)
    //         {
    //             resP[i] = P[j] * mat[j][i];
    //         }            
    //     }
    //     return res;
    // }
    // //---
    
    //---------------------------------------------------
    // MULTIPLY BY SCALAR
    //---------------------------------------------------
    VECTOR operator * (prec scalar)
    {
        VECTOR res(length);
        std::shared_ptr<prec[]> resP = res.P;
        for (int i = 0; i < length; i++)
        {
            resP[i] = P[i] * scalar;
        }
        return res;
    }
    //---
    void multiplyByScalar(prec scalar, VECTOR &res)
    {
        for (int i = 0; i < length; i++)
        {
            res.P[i] = P[i]*scalar;
        }
    }
    //---
    void multiplyByScalar(prec scalar,  std::shared_ptr<prec[]> &res)
    {
        if(res == 0) 
        {
            res = std::shared_ptr<prec[]>(new prec[length]);
        }
        for (int i = 0; i < length; i++)
        {
            res[i] = P[i]*scalar;
        }
    }
    //---
    static void multiplyByScalar(std::shared_ptr<prec[]> &vec, prec scalar, int N, std::shared_ptr<prec[]> &res)
    {
        if(res == 0) 
        {
            res = std::shared_ptr<prec[]>(new prec[N]);
        }
        for (int i = 0; i < N; i++)
        {
            res[i] = vec[i]*scalar;
        }
    }
    //---
    static void multiplyByScalar(prec* &vec, prec scalar, int N, std::shared_ptr<prec[]> &res)
    {
        if(res == 0) 
        {
            res = std::shared_ptr<prec[]>(new prec[N]);
        }
        for (int i = 0; i < N; i++)
        {
            res[i] = vec[i]*scalar;
        }
    }
    //---
    static void multiplyByScalar(std::shared_ptr<prec[]> &vec, prec scalar, int N, VECTOR &res)
    {
        std::shared_ptr<prec[]> P_in = res.P;
        for (int i = 0; i < N; i++)
        {
            P_in[i] = vec[i]*scalar;
        }
    }
    //---
    static void multiplyByScalar(prec* &vec, prec scalar, int N, VECTOR &res)
    {
        std::shared_ptr<prec[]> P_in = res.P;
        for (int i = 0; i < N; i++)
        {
            P_in[i] = vec[i]*scalar;
        }
    }

    //------------------------------------------------------
    // CROSS PRODUCT IN 3D
    //------------------------------------------------------
    static VECTOR cross(VECTOR &v1, VECTOR &v2)
    {
        if (v1.length != 3 || v2.length != 3) throw_line("ERROR: wrong dimension in cross product\n");
        VECTOR res(3);
        res[0] = (v1[1] * v2[2]) - (v1[2] * v2[1]); if (std::abs(res[0]) < 1e-16) res[0] = 0.0;
        res[1] = (v1[2] * v2[0]) - (v1[0] * v2[2]); if (std::abs(res[1]) < 1e-16) res[1] = 0.0;
        res[2] = (v1[0] * v2[1]) - (v1[1] * v2[0]); if (std::abs(res[2]) < 1e-16) res[2] = 0.0;
        return res;
    }
    //---
    static void cross(VECTOR &v1, VECTOR &v2, VECTOR &res)
    {
        if (v1.length != 3 || v2.length != 3) throw_line("ERROR: wrong dimension in cross product\n");
        res.P.reset();
        res.P = 0;
        res.initialize(3);
        res[0] = (v1[1] * v2[2]) - (v1[2] * v2[1]); if (std::abs(res[0]) < 1e-16) res[0] = 0.0;
        res[1] = (v1[2] * v2[0]) - (v1[0] * v2[2]); if (std::abs(res[1]) < 1e-16) res[1] = 0.0;
        res[2] = (v1[0] * v2[1]) - (v1[1] * v2[0]); if (std::abs(res[2]) < 1e-16) res[2] = 0.0;
    }
    //---
    VECTOR cross(VECTOR &v2)
    {
        if (length != 3 || v2.length != 3) throw_line("ERROR: wrong dimension in cross product\n");
        VECTOR res(3);
        res[0] = (P[1] * v2[2]) - (P[2] * v2[1]); if (std::abs(res[0]) < 1e-16) res[0] = 0.0;
        res[1] = (P[2] * v2[0]) - (P[0] * v2[2]); if (std::abs(res[1]) < 1e-16) res[1] = 0.0;
        res[2] = (P[0] * v2[1]) - (P[1] * v2[0]); if (std::abs(res[2]) < 1e-16) res[2] = 0.0;
        return res;
    }
    //---
    static VECTOR cross(prec* &v1, prec* &v2)
    {
        VECTOR res(3);
        res[0] = (v1[1] * v2[2]) - (v1[2] * v2[1]); if (std::abs(res[0]) < 1e-16) res[0] = 0.0;
        res[1] = (v1[2] * v2[0]) - (v1[0] * v2[2]); if (std::abs(res[1]) < 1e-16) res[1] = 0.0;
        res[2] = (v1[0] * v2[1]) - (v1[1] * v2[0]); if (std::abs(res[2]) < 1e-16) res[2] = 0.0;
        return res;
    }

    bool check_nan(bool stop = false)
    {
        bool there_are_nan = false;
        for (int i = 0; i < length; i++)
        {
            if (std::isnan(P[i]))
            {
                std::cout << "\ni: " << i << " is nan";
                there_are_nan = true;
            }
        }
        if (stop && there_are_nan)
        {
            throw_line("ERROR: vector contains not wanted nan values\n");
        }
        return there_are_nan;
    }

    //------------------------------------------------------
    // CHECK IF THERE IS AN ELEMENT IN A VECTOR
    //------------------------------------------------------
    bool hasIn(prec value)
    {
        bool answer = false;
        int i;
        for (i = 0; i < length; i++)
        {
            if (P[i] == value)
            {
                answer = true;
                return answer;
            } 
        }
        return answer;
    }
    //---
    bool isIn(prec value, int* idAddress)
    {
        bool answer = false;
        int i;
        for (i = 0; i < length; i++)
        {
            if (P[i] == value)
            {
                answer = true;
                break;
            } 
        }
        *idAddress = i;
        return answer;
    }
    //---
    static bool isIn(std::shared_ptr<prec[]> &V_in, int N, prec value, int idAddress)
    {
        bool answer = false;
        int i;
        for (i = 0; i < N; i++)
        {
            if (V_in[i] == value)
            {
                answer = true;
                idAddress = i;
                break;
            } 
        }
        
        return answer;
    }
    //---
    static bool isIn(std::shared_ptr<int[]> &V_in, int N, int value, int &idAddress)
    {
        bool answer = false;
        int i;
        for (i = 0; i < N; i++)
        {
            if (V_in[i] == value)
            {
                answer = true;
                idAddress = i;
                break;
            } 
        }       
        return answer;
    }
    //------------------------------------------------------
    // INVERT
    //------------------------------------------------------
    void invert()
    {
        for (int i = 0; i < length; i++)
        {
            if (P[i] == 0) throw_line("ERROR, inverting a vector with zero entries\n");
            P[i] = 1/P[i];
        } 
    } 

    //------------------------------------------------------
    // REVERSE THE VECTOR COMPONENTS
    //------------------------------------------------------
    void reverse()
    {
        std::shared_ptr<prec[]> temp_P = 0;
        temp_P = std::shared_ptr<prec[]>(new prec[length]);
        for (int i = 0; i < length; i++)
        {
            temp_P[i] = P[i];
        }
        for (int i = 0; i < length; i++)
        {
            P[i] = temp_P[length-1-i];
        }
    }
    

    //------------------------------------------------------
    // VECTOR NORM
    //------------------------------------------------------
    prec norm()
    {
        prec res = 0;
        for (int i = 0; i < length; i++)
        {
                res += P[i] * P[i];
        }
        res = sqrt(res);
        return res;
    }
    //---
    static prec norm(VECTOR vec)
    {
        prec res = 0;
        std::shared_ptr<prec[]> P_in = vec.P;
        for (int i = 0; i < vec.length; i++)
        {
                res += P_in[i] * P_in[i];
        }
        res = sqrt(res);
        return res;
    }
    //
    static prec norm(std::shared_ptr<prec[]> &vec, int N)
    {
        prec res = 0;
        for (int i = 0; i < N; i++)
        {
                res += vec[i] * vec[i];
        }
        res = sqrt(res);
        return res;
    }

    //-----------------------------------------------------
    // NORM WITH BEGIN & END
    //-----------------------------------------------------
    prec norm(int begin, int end)
    {
        if (begin < 0 || end > length) throw_line("ERROR: The indices are not valid \n");
        prec res = 0;
        for (int i = begin; i < end; i++)
        {
                res += P[i] * P[i];
        }
        res = sqrt(res);
        return res;
    }

    //------------------------------------------------------
    // COPY VECTOR INTO OTHER VECTOR
    //------------------------------------------------------
    void exact_copy(VECTOR inVec)
    {
        length = inVec.length;
        P.reset(); P = 0;
        P = std::shared_ptr<prec[]>(new prec[length]);
        for (int i = 0; i < length; i++)
        {
            P[i] = inVec[i];
        }
    }

    void operator = (VECTOR inVec) //(works also to copy vectors of different sizes)
    {
        //if (length != inVec.length) throw_line("ERROR: the vectors have different lengths \n;")
        if (length == 0)
        {
            length = inVec.length;
            P = std::shared_ptr<prec[]>(new prec[length]);
        }

        std::shared_ptr<prec[]> P_in = inVec.P;
        for (int i = 0; i < length; i++)
        {
            P[i] = P_in[i];
        }
    }
    //---
    void operator = (std::shared_ptr<prec[]> &inVec)
    {
        for (int i = 0; i < length; i++)
        {
            P[i] = inVec[i];
        }
    }
    //---
    void operator = (prec* &inVec)
    {
        for (int i = 0; i < length; i++)
        {
            P[i] = inVec[i];
        }
    }
    //---
    void copyTo(VECTOR &res)
    {
        std::shared_ptr<prec[]> P_in = res.P;
        for (int i = 0; i < length; i++) P_in[i] = P[i];
    }
    //---
    void copyTo(std::shared_ptr<prec[]> &res)
    {
        if(res == 0)
        {
            res = std::shared_ptr<prec[]>(new prec[length]);
        } 
        for (int i = 0; i < length; i++) res[i] = P[i];
    }
    //---
    void copyTo(prec* &res)
    {
        if(res == 0 || res == nullptr)
        {
            throw_line("ERROR of some type\n");
        } 
        for (int i = 0; i < length; i++) res[i] = P[i];
    }
    //---
    void copy(std::shared_ptr<prec[]> &vec)
    {
        if(vec == nullptr) throw_line("ERROR: COPYING NULL PTR\n");
        for (int i = 0; i < length; i++) P[i] = vec[i];
    }
    //---
    void copy(prec* &vec)
    {
        if(vec == nullptr) throw_line("ERROR: COPYING NULL PTR\n");
        for (int i = 0; i < length; i++) P[i] = vec[i];
    }
    //---
    static void copy(std::shared_ptr<prec[]> &inVec, std::shared_ptr<prec[]> &outVec, int N, int startPos = 0)
    {
        for (int i = 0; i < N; i++) outVec[i+startPos] = inVec[i];
    }
    //---
    static void copy(std::shared_ptr<prec[]> &inVec, prec* &outVec, int N, int startPos = 0)
    {
        for (int i = 0; i < N; i++) outVec[i+startPos] = inVec[i];
    }
    //---
    static void copy(prec* &inVec, prec* &outVec, int N, int startPos = 0)
    {
        for (int i = 0; i < N; i++) outVec[i+startPos] = inVec[i];
    }
    //---
    static void copy(std::shared_ptr<int[]> &inVec, std::shared_ptr<int[]> &outVec, int N, int startPos = 0)
    {
        for (int i = 0; i < N; i++) outVec[i+startPos] = inVec[i];
    }
    //----------------------------------------

    //----------------------------------------
    // CHECK IF TWO VECTORS ARE EQUAL
    //----------------------------------------
    static bool equal(VECTOR vec1, VECTOR vec2)
    {
        bool state = true;
        if (vec1.length != vec2.length)
        {
            state = false;
        }
        else
        {
            int vecL = vec1.length;
            for (int i = 0; i < vecL; i++)
            {
                if (vec1[i] != vec2[i])
                {
                    state = false;
                    break;
                }
            }
        }
        return state;
    }
    //----------------------------------------

    //----------------------------------------
    // MAX AND MIN
    //----------------------------------------
    VECTOR vecAbs()
    {
        VECTOR outVec(length);
        std::shared_ptr<prec[]> outP = outVec.P;
        for (int i = 0; i < length; i++)
        {
            prec value = P[i];
            outP[i] = abs(value);
        }
        return outVec;
    }
    //---
    prec max()
    {
        prec incumbentMax = -1e50;
        for (int i = 0; i < length; i++)
        {
            if (P[i] > incumbentMax) 
            {
                incumbentMax = P[i];
            }
        }
        return incumbentMax;
    }
    //---
    prec min()
    {
        prec incumbentMin = 1e50;
        for (int i = 0; i < length; i++)
        {
            if (P[i] < incumbentMin) 
            {
                incumbentMin = P[i];
            }
        }
        return incumbentMin;
    }

    //----------------------------------------
    // APPEND VECTOR
    //----------------------------------------
    void append1(prec coef)
    {
        int newL = length + 1;
        std::shared_ptr<prec[]> newP(new prec[newL]);

        for (int i = 0; i < length; i++) newP[i] = P[i];
        newP[length] = coef;

        P.reset(); P = 0;
        P = newP;
        length = newL;
    }
    //---
    void append(VECTOR &vecIn)
    {
        int newL = length + vecIn.length;

        std::shared_ptr<prec[]> newP(new prec[newL]);

        for (int i = 0; i < length; i++) newP[i] = P[i];

        copy(vecIn.P, newP, vecIn.length, length);

        P.reset(); P = 0;
        P = newP;
        length = newL;
    }
    //---
    void append (std::shared_ptr<prec[]> &vecIn, int vecInL)
    {
        int newL = length + vecInL;

        std::shared_ptr<prec[]> newP(new prec[newL]);

        for (int i = 0; i < length; i++) newP[i] = P[i];
        copy(vecIn, newP, vecInL, length);

        P.reset(); P = 0;
        P = newP;
        length = newL;
    }
    //---
    static void append (VECTOR &vec, VECTOR vecIn)
    {
        vec.append(vecIn);
    }
    //---
    static void append (std::shared_ptr<prec[]> &vec, int vecL, std::shared_ptr<prec[]> &vecIn, int vecInL)
    {
        int newL = vecL + vecInL;
        std::shared_ptr<prec[]> newP(new prec[newL]);
        for (int i = 0; i < vecL; i++) newP[i] = vec[i];
        copy(vecIn, newP, vecInL, vecL);

        vec.reset(); vec = 0;
        vec = newP;
    }
    //---
    //-------------------------------
    // APPEND ONE VALUE
    //-------------------------------
    void append (prec val)
    {
        int newL = length + 1;

        std::shared_ptr<prec[]> newP(new prec[newL]);

        for (int i = 0; i < length; i++) newP[i] = P[i];
        newP[length]= val;

        P.reset(); P = 0;
        P = newP;
        length = newL;
    }
    //-----------------------------------------------------
    // SHRINK VECTOR (keep only newL elements of the array)
    //-----------------------------------------------------
    void shrink(int newL)
    {
        std::shared_ptr<prec[]> newP(new prec[newL]);

        copy(P, newP, newL);
         
        P.reset(); P = 0;
        P = newP;
        length = newL;
    }
    //---
    static void shrink(std::shared_ptr<prec[]> &vec, int newL)
    {
        std::shared_ptr<prec[]> newP(new prec[newL]);

        copy(vec, newP, newL);
        
        vec.reset(); vec = 0;
        vec = newP;
    }
     //-----------------------------------------------------
    // ENLARGE VECTOR 
    //-----------------------------------------------------
    void enlarge(int newL)
    {
        if (newL < length) throw_line("ERROR: enlarging vector into one with smaller size\n");
        std::shared_ptr<prec[]> newP(new prec[newL]);
        copyTo(newP);
        
        P.reset(); P = 0;
        P = newP;
        length = newL;
    }

    //-------------------------------------
    // UNIFY VECTORS
    //-------------------------------------
    static std::shared_ptr<prec[]> merge (std::shared_ptr<prec[]> &vec1, int vec1L, std::shared_ptr<prec[]> &vec2, int vec2L)
    {
        int newL = vec1L + vec2L;

        std::shared_ptr<prec[]> res(new prec[newL]);

        copy(vec1, res, vec1L);
        copy(vec2, res , vec2L, vec1L);

        return res;
    }
    //---  


    //-------------------------------------------------------
    // PRINT VECTOR
    //-------------------------------------------------------
    void print()
    {
        std::cout << "\n-------\n";
        for (int i = 0; i < length; i++)
        {
            std::cout << " " << P[i] << "\n";
        }
        std::cout << "-------\n";
    }
    //---
    void print(std::string vecName)
    {
        std::cout << "\n-------\n";
        std::cout << vecName;
        std::cout << "\n-------\n";
        for (int i = 0; i < length; i++)
        {
            std::cout << " " << P[i] << "\n";
        }
        std::cout << "-------\n";
    }
    //---
    static void print(std::shared_ptr<prec[]> &inVec, int N)
    {
        printf("\n-------\n");
        for (int i = 0; i < N; i++)
        {
            printf("%10.8" format "\n", inVec[i]);
        }
        printf("-------\n");
    }
    //---
    static void print(prec* &inVec, int N)
    {
        printf("\n-------\n");
        for (int i = 0; i < N; i++)
        {
            printf("%10.8" format "\n", inVec[i]);
        }
        printf("-------\n");
    }
    //---
    static void print(int* &inVec, int N)
    {
        printf("\n-------\n");
        for (int i = 0; i < N; i++)
        {
            printf("%d\n", inVec[i]);
        }
        printf("-------\n");
    }
    //---
    static void print(std::shared_ptr<prec[]> &inVec, int N, std::string vecName)
    {
        std::cout << "\n-------\n";
        std::cout << vecName;
        std::cout << "\n-------\n";
        for (int i = 0; i < N; i++)
        {
            printf("%10.8" format_e "\n", inVec[i]);
        }
        printf("-------\n");
    }
    //---
    static void print(std::shared_ptr<int[]> &inVec, int N)
    {
        std::cout << "\n-------\n";
        for (int i = 0; i < N; i++)
        {
            std::cout << " " << inVec[i] << "\n";
        }
        std::cout << "-------\n";
    }
    //---
    static void print(std::shared_ptr<int[]> &inVec, int N, std::string vecName)
    {
        std::cout << "\n-------\n";
        std::cout << vecName;
        std::cout << "\n-------\n";
        for (int i = 0; i < N; i++)
        {
            std::cout << " " << inVec[i] << "\n";
        }
        std::cout << "-------\n";
    }
    //---
    static void printFile(std::shared_ptr<prec[]> &vec, int N, const char* outFileName)
    {
        FILE* outFile = fopen(outFileName,"w");
        for (int i = 0; i < N; i++)
        {
            fprintf(outFile, "%20.16" format "\n", vec[i]);
        }
        fclose(outFile);
    }
    //---
    static void printFile(std::shared_ptr<int[]> &vec, int N, const char* outFileName)
    {
        FILE* outFile = fopen(outFileName,"w");
        for (int i = 0; i < N; i++)
        {
            fprintf(outFile, "%5d \n", vec[i]);
        }
        fclose(outFile);
    }
    //---
    void print(const char* outFileName)
    {

        FILE* outFile = fopen(outFileName,"w");
        
        //Print matrix;
        for (int i = 0; i < length; i++)
        {
            fprintf(outFile, "%15.11" format "\n", P[i]);
        }
        fclose(outFile);
    }
    //---
    
    void printRow()
    {
        for (int i = 0; i < length; i++)
        {
            std::cout << P[i] << " ";
        }
        std::cout << "\n";
    }
    //---
    void printRow(std::string nameVec)
    {
        std::cout << nameVec << ": ";
        for (int i = 0; i < length; i++)
        {
            std::cout << P[i] << " ";
        }
        std::cout << "\n";
    }
    //---
    void printRowMatlab()
    {
        if (length == 0) std::cout << "[]\n";
        else
        {
            std::cout << "[";
            for (int i = 0; i < length-1; i++)
            {
                std::cout << P[i] << ", ";
            }
            std::cout << P[length-1] << "]\n";
        }
    }
    //---
    void printRowMatlab(std::string nameVec)
    {
        std::cout << nameVec << ": ";
        if (length == 0) std::cout << "[]\n";
        else
        {
            std::cout << "[";
            for (int i = 0; i < length-1; i++)
            {
                std::cout << P[i] << ", ";
            }
            std::cout << P[length-1] << "]\n";
        }
    }
    //---
    //---
    static void printRowMatlab(std::shared_ptr<prec[]> &vec, int N)
    {
        if (N == 0) std::cout << "[]\n";
        else
        {
            std::cout << "[";
            for (int i = 0; i < N-1; i++)
            {
                std::cout << vec[i] << ", ";
            }
            std::cout << vec[N-1] << "]\n";
        }
    }
    //---
    static void printRowMatlab(std::shared_ptr<prec[]> &vec, int N, std::string nameVec)
    {
        std::cout << nameVec << ": ";
        if (N == 0) std::cout << "[]\n";
        else
        {
            std::cout << "[";
            for (int i = 0; i < N-1; i++)
            {
                std::cout << vec[i] << ", ";
            }
            std::cout << vec[N-1] << "]\n";
        }
    }
    //---
    void printRowMatlabInLine()
    {
        if (length == 0) std::cout << "[]";
        else
        {
            std::cout << "[";
            for (int i = 0; i < length-1; i++)
            {
                std::cout << P[i] << ", ";
            }
            std::cout << P[length-1] << "]";
        }
    }
    //---
    void printRowMatlabInLine(std::string nameVec)
    {
        std::cout << nameVec << ": ";
        if (length == 0) std::cout << "[]";
        else
        {
            std::cout << "[";
            for (int i = 0; i < length-1; i++)
            {
                std::cout << P[i] << ", ";
            }
            std::cout << P[length-1] << "]";
        }
    }
    //---
    static void printRowMatlabInLine(prec* &vec, int N)
    {
        if (N == 0) std::cout << "[]";
        else
        {
            std::cout << "[";
            for (int i = 0; i < N-1; i++)
            {
                std::cout << vec[i] << ", ";
            }
            std::cout << vec[N-1] << "]";
        }
    }
    //---
    static void printRowMatlabInLine(prec* &vec, int N, std::string nameVec)
    {
        std::cout << nameVec << ": ";
        if (N == 0) std::cout << "[]";
        else
        {
            std::cout << "[";
            for (int i = 0; i < N-1; i++)
            {
                std::cout << vec[i] << ", ";
            }
            std::cout << vec[N-1] << "]";
        }
    }
    //---

    //-----------------------------------------------
    // PRINT 2 POINTERS SIDE BY SIDE (4 int, prec dispositions)
    //-----------------------------------------------
    static void print2(std::shared_ptr<int[]> &vec1, std::shared_ptr<int[]>  &vec2, int N, std::string name1, std::string name2)
    {
        std::cout << "\n-------\t\t" << "-------\n";
        std::cout << name1 << "\t\t" << name2 << "\n";
        std::cout << "-------\t\t" << "-------\n";
        for (int i = 0; i < N; i++)
        {
            //std::cout << " " << vec1[i] << "\t\t " << vec2[i] << "\n";
            printf(" %4.2d \t\t %4.2d\n", vec1[i], vec2[i]);
        }
        std::cout << "-------\t\t" << "-------\n";
    }
    //---
    static void print2(std::shared_ptr<int[]>  &vec1, std::shared_ptr<prec[]>  &vec2, int N, std::string name1, std::string name2)
    {
        std::cout << "\n-------\t\t" << "-------\n";
        std::cout << name1 << "\t\t" << name2 << "\n";
        std::cout << "-------\t\t" << "-------\n";
        for (int i = 0; i < N; i++)
        {
            //std::cout << " " << vec1[i] << "\t\t " << vec2[i] << "\n";
            printf(" %4.2d\t\t %4.2" format "\n", vec1[i], vec2[i]);
        }
        std::cout << "-------\t\t" << "-------\n";
    }
    //---
    static void print2(std::shared_ptr<prec[]>  &vec1, std::shared_ptr<int[]>  &vec2, int N, std::string name1, std::string name2)
    {
        std::cout << "\n-------\t\t" << "-------\n";
        std::cout << name1 << "\t\t" << name2 << "\n";
        std::cout << "-------\t\t" << "-------\n";
        for (int i = 0; i < N; i++)
        {
            //std::cout << " " << vec1[i] << "\t\t " << vec2[i] << "\n";
            printf(" %4.2" format "\t\t %4.2d \n", vec1[i], vec2[i]);
        }
        std::cout << "-------\t\t" << "-------\n";
    }
    //---
    static void print2(std::shared_ptr<prec[]> &vec1, std::shared_ptr<prec[]> &vec2, int N, std::string name1, std::string name2)
    {
        std::cout << "\n-------\t\t" << "-------\n";
        std::cout << name1 << "\t\t" << name2 << "\n";
        std::cout << "-------\t\t" << "-------\n";
        for (int i = 0; i < N; i++)
        {
            //std::cout << " " << vec1[i] << "\t\t " << vec2[i] << "\n";
            printf(" %9.7" format "\t\t %9.7" format "\n", vec1[i], vec2[i]);
        }
        std::cout << "-------\t\t" << "-------\n";
    }
    //---
    
    //----------------
    // PRINT IN ROW
    //----------------
    void printRowPure()
    {
        for (int i = 0; i < length-1; i++)
        {
            printf("%.2" format "   ", P[i]);
        }
        printf("%.2" format, P[length-1]);
    }

    void print_with_ids()
    {
        for (int i = 0; i < length-1; i++)
        {
            printf("\n%d: %.9" format , i, P[i]);
        }
        printf("\n%d: %.9" format "\n", length-1, P[length-1]);
    }
    
    //--------------------------------------------
    // CLEAR VECTOR
    //--------------------------------------------
    void clear()
    {
        P.reset(); P = 0;
        length = 0;
    }
    //---

    void completeReset()
    {
        P.reset(); P = 0;
        length = 0;
    }
    //---

    void complete_reset()
    {
        P.reset(); P = 0;
        length = 0;
    }
    //---



};