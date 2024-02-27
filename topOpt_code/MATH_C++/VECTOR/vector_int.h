#pragma once

#include "../../CODE_HEADERS/codeHeader.h"

class VECTOR_INT
{
public:
    std::shared_ptr<int[]> P = 0;
    int length = 0;
    //-------------------------------------------------------
    // CONSTRUCTOR
    //-------------------------------------------------------
    // EMPTY CONSTRUCTOR
    VECTOR_INT() 
    {
        P = 0;
    }	
    //---
    VECTOR_INT(int N)
    {   
        length = N;
        P = std::shared_ptr<int[]>(new int[N]);
    }
    //---
    VECTOR_INT(int* &vecPointer, int N)
    {  
        P = std::shared_ptr<int[]>(new(vecPointer) int[N], [](int* vecPointer){free(vecPointer);});
        length = N;
    }
    //---
    VECTOR_INT( std::shared_ptr<int[]> &vecPointer, int N)
    {  
        P = vecPointer;
        length = N;
    }
    //---
    // VECTOR_INT(VECTOR_INT &inVec)
    // {
    //     length = inVec.length;
    //     int* P_in = inVec.P;
    //     P = (int*) malloc(length * sizeof(int));
    //     for (int i = 0; i < length; i++)
    //     {
    //         P[i] = P_in[i];
    //     }
    // }
    ~VECTOR_INT()
    {}

    //-----------------------
    // INITIALIZE
    //-----------------------
    void initialize(int N)
    {
        length = N;
        P = std::shared_ptr<int[]>(new int[N]);
    }
    //---
    void initializeZero(int N)
    {
        length = N;
        P = std::shared_ptr<int[]>(new int[N]);
        for (int i = 0; i < N; i++) P[i] = 0;
    }
    //---
    static std::shared_ptr<int[]> makePointer(int N)
    {
        std::shared_ptr<int[]> Pout(new int[N]);

        return Pout;
    }
    //---
    static std::shared_ptr<int*[]> makeDoublePointer(int N)
    {
        std::shared_ptr<int*[]> Pout(new int*[N]);
        return Pout;
    }
    //---
    static VECTOR_INT zeros(int N)
    {
        VECTOR_INT vec(N);
        for (int i = 0; i < N; i++) vec[i] = 0;
        return vec;
    }
    //---
    static void zeros(VECTOR_INT & vec, int N)
    {
        vec.initialize(N);
        for (int i = 0; i < N; i++) vec[i] = 0;
    }
    //---
    void zeros()
    {
        for (int i = 0; i < length; i++) P[i] = 0;
    }
    //---
    void setZeros(int N)
    {
        P.reset(); P = 0;
        length = N;
        P = std::shared_ptr<int[]>(new int[N]);
        for (int i = 0; i < length; i++) P[i] = 0;
    }
    //---
    void reset(int coef)
    {
        for (int i = 0; i < length; i++) P[i] = coef;
    }
    //---
    void define(int N, std::shared_ptr<int[]> &inVec)
    {
        initialize(N);
        for (int i = 0; i < length; i++)
        {
            P[i] = inVec[i];
        }
    }
    //-------------------------------------------------------
    // GET SINGLE ENTRY
    //-------------------------------------------------------
    int& operator [] (int index)
    {
        if (index >= length || index < 0) throw_line("ERROR: index out of bounds\n");
        return P[index];
    }
    //---
    void set3Values(int val0, int val1, int val2)
    {
        P.reset(); P = 0;
        initialize(3);
        P[0] = val0; P[1] = val1; P[2] = val2;
    }

    void set4Values(int val0, int val1, int val2, int val3)
    {
        P.reset(); P = 0;
        initialize(4);
        P[0] = val0; P[1] = val1; P[2] = val2; P[3] = val3;
    }
    //-----------------
    static void set2Values(std::shared_ptr<int[]> &P, int val0, int val1)
    {
        P[0] = val0; P[1] = val1;
    }

    static void set2Values(int* &P, int val0, int val1)
    {
        P[0] = val0; P[1] = val1;
    }

    static void set3Values(std::shared_ptr<int[]> &P, int val0, int val1, int val2)
    {
        P[0] = val0; P[1] = val1; P[2] = val2;
    }
    static void set3Values(int* &P, int val0, int val1, int val2)
    {
        P[0] = val0; P[1] = val1; P[2] = val2;
    }

    static void set4Values(std::shared_ptr<int[]> &P, int val0, int val1, int val2, int val3)
    {
        P[0] = val0; P[1] = val1; P[2] = val2; P[3] = val3;
    }
    static void set4Values(int* &P, int val0, int val1, int val2, int val3)
    {
        P[0] = val0; P[1] = val1; P[2] = val2; P[3] = val3;
    }

    //-------------------------------------------------------
    // DOT PRODUCT
    //-------------------------------------------------------
    int dot(VECTOR_INT &inVec)
    {
        if (length != inVec.length) printf("ERROR: the vectors have different lengths \n;");
        int res = 0;
        std::shared_ptr<int[]> pIn = inVec.P;
        for (int i = 0; i < length; i++)
        {
            res += P[i] * pIn[i];
        }
        return res;
    } 
    //---
    int dot(std::shared_ptr<int[]> &inVec)
    {
        int res = 0;
        for (int i = 0; i < length; i++)
        {
            res += P[i] * inVec[i];
        }
        return res;
    } 

    //-----------------------------------------------------
    // DOT PRODUCT OF POINTERS UNTIL COMPONENT N
    //-----------------------------------------------------  
    static int dot(std::shared_ptr<int[]> &vec1, std::shared_ptr<int[]> &vec2, int N)
    {
        int res = 0;
        for (int i = 0; i < N; i++)
        {
            res += vec1[i] * vec2[i];
        }
        return res;
    }

    //-----------------------------------------------------
    // DOT PRODUCT WITH BEGIN & END
    //-----------------------------------------------------
    int dot(VECTOR_INT &inVec, int begin, int end)
    {
        if (length != inVec.length) printf("ERROR: the vectors have different lengths \n;");
        if (begin < 0 || end > length) printf("ERROR: The indices are not valid \n");
        int res = 0;
        std::shared_ptr<int[]>  P_in = inVec.P;
        for (int i = begin; i < end; i++)
        {
            res += P[i] * P_in[i];
        }
        return res;
    } 
    
    //-----------------------------------------------------
    // POINTWISE VECTOR PRODUCT
    //-----------------------------------------------------
    VECTOR_INT operator * (VECTOR_INT &inVec)
    {
        if (length != inVec.length) printf("ERROR: the vectors have different lengths \n;");
        VECTOR_INT res(length);    
        std::shared_ptr<int[]> resP = res.P;        
        std::shared_ptr<int[]>  P_in = inVec.P;
        for (int i = 0; i < length; i++)
        {
            resP[i] = P[i] * P_in[i];
        }
        return res;
    }
    //---
    void pointdot(int* &inVec, VECTOR_INT &res)
    {   
        if (length != res.length) printf("ERROR: the vectors have different lengths \n;");
        for (int i = 0; i < length; i++)
        {
            res.P[i] = P[i] * inVec[i];
        }
    }
    //---
    void pointdot(std::shared_ptr<int[]>  &inVec, std::shared_ptr<int[]>  &res)
    {   
        if(res == 0)
        {
            res = std::shared_ptr<int[]>(new int[length]);
        }

        for (int i = 0; i < length; i++)
        {
            res[i] = P[i] * inVec[i];
        }
    }
    //---
    static void pointdot(std::shared_ptr<int[]> &vec1, std::shared_ptr<int[]> &vec2, int N, std::shared_ptr<int[]> &res)
    {
        if(res == 0)
        {
            res = std::shared_ptr<int[]>(new int[N]);
        } 
        for (int i = 0; i < N; i++)
        {
            res[i] = vec1[i] * vec2[i];
        }
    }

    //----------------------------------------------------------
    // SELF OPERATIONS
    //----------------------------------------------------------
    void operator += (VECTOR_INT &inVec)
    {
        if (length != inVec.length) throw_line("ERROR: the vectors have different lengths \n;")
        std::shared_ptr<int[]> P_in = inVec.P;
        for (int i = 0; i < length; i++)
        {
            P[i] += P_in[i];
        }
    }
    //------------------------
    void operator += (int inScalar)
    {
        for (int i = 0; i < length; i++)
        {
            P[i] += inScalar;
        }
    }

    void operator -= (VECTOR_INT &inVec)
    {
        if (length != inVec.length) throw_line("ERROR: the vectors have different lengths \n;")
        std::shared_ptr<int[]> P_in = inVec.P;
        for (int i = 0; i < length; i++)
        {
            P[i] -= P_in[i];
        }
    }

    void operator /= (int coef)
    {
        for (int i = 0; i < length; i++)
        {
            P[i] /= coef;
        }
    }

    void operator *= (int coef)
    {
        for (int i = 0; i < length; i++)
        {
            P[i] *= coef;
        }
    }

    void operator *= (VECTOR_INT &vec)
    {
        if (vec.length != length) throw_line("ERROR: Pointwise product of vectors of different lenght\n");
        for (int i = 0; i < length; i++)
        {
            P[i] *= vec[i];
        }
    }

    void operator /= (VECTOR_INT &vec)
    {
        if (vec.length != length) throw_line("ERROR: Pointwise product of vectors of different lenght\n");
        for (int i = 0; i < length; i++)
        {
            if (vec[i] == 0) throw_line("ERROR: Dividing by zero in pointwise divsion of vectors\n");
            P[i] /= vec[i];
        }
    }
    
    //----------------------------------------------------------
    // VECTOR_INT SUM
    //----------------------------------------------------------
    VECTOR_INT operator + (VECTOR_INT &inVec)
    {
        if (length != inVec.length) throw_line("ERROR: the vectors have different lengths \n;")
        VECTOR_INT res(length);
        std::shared_ptr<int[]> resP = res.P;
        std::shared_ptr<int[]> P_in = inVec.P;
        for (int i = 0; i < length; i++)
        {
            resP[i] = P[i] + P_in[i];
        }
        return res;
    }
    //----------------------------------------------------------
    VECTOR_INT operator + (int coef)
    {
        VECTOR_INT res(length);
        std::shared_ptr<int[]> resP = res.P;
        for (int i = 0; i < length; i++)
        {
            resP[i] = P[i] + coef;
        }
        return res;
    }
    //---
    int sum()
    {
        int res = 0;
        for (int i = 0; i < length; i++) res += P[i];
        return res;
    }
    //---
    void sum(std::shared_ptr<int[]> &inVec, VECTOR_INT &res)
    {
        for (int i = 0; i < length; i++)
        {
            res.P[i] = P[i] + inVec[i];
        }
    }
    //---
    void sum(std::shared_ptr<int[]> &inVec, std::shared_ptr<int[]> &res)
    {
        if(res == 0) 
        {
            res = std::shared_ptr<int[]>(new int[length]);
        }
        for (int i = 0; i < length; i++)
        {
            res[i] = P[i] + inVec[i];
        }
    }
    //---
    static void sum(std::shared_ptr<int[]> &vec1, std::shared_ptr<int[]> &vec2, int N, std::shared_ptr<int[]> &res)
    {
        if(res == 0) 
        {
            res = std::shared_ptr<int[]>(new int[N]);
        }
        for (int i = 0; i < N; i++)
        {
            res[i] = vec1[i] + vec2[i];
        }
    }

    //----------------------------------------------------------
    // VECTOR DIFFERENCE
    //----------------------------------------------------------
    VECTOR_INT operator - (VECTOR_INT &inVec)
    {
        if (length != inVec.length) throw_line("ERROR in VECTOR DIFFERENCE: the vectors have different lengths \n;")
        VECTOR_INT res(length);
        std::shared_ptr<int[]> resP = res.P;
        std::shared_ptr<int[]> P_in = inVec.P;
        for (int i = 0; i < length; i++)
        {
            resP[i] = P[i] - P_in[i];
        }
        return res;
    }
    //---
    VECTOR_INT operator - (int coef)
    {
        VECTOR_INT res(length);
        std::shared_ptr<int[]> resP = res.P;
        for (int i = 0; i < length; i++)
        {
            resP[i] = P[i] - coef;
        }
        return res;
    }
    //---
    void diff(std::shared_ptr<int[]> &inVec, VECTOR_INT &res)
    {
        for (int i = 0; i < length; i++)
        {
            res.P[i] = P[i] - inVec[i];
        }
    }
    //---
    void diff(std::shared_ptr<int[]> &inVec, std::shared_ptr<int[]> &res)
    {
        if(res == 0) 
        {
            res = std::shared_ptr<int[]>(new int[length]);
        }
        for (int i = 0; i < length; i++)
        {
            res[i] = P[i] - inVec[i];
        }
    }
    //---
    static void diff(std::shared_ptr<int[]> &vec1, std::shared_ptr<int[]> &vec2, int N, std::shared_ptr<int[]> &res)
    {
        if(res == 0) 
        {
            res = std::shared_ptr<int[]>(new int[N]);
        }
        for (int i = 0; i < N; i++)
        {
            res[i] = vec1[i] - vec2[i];
        }
    }

    //---------------------------------------------------------------
    // DIVIDE BY SCALAR
    //---------------------------------------------------------------
    VECTOR_INT operator / (int scalar)
    {
        VECTOR_INT res(length);
        std::shared_ptr<int[]> resP = res.P;
        for (int i = 0; i < length; i++)
        {
            resP[i] = P[i] / scalar;
        }
        return res;
    }
    //---
    void divideByScalar(int scalar, VECTOR_INT &res)
    {
        for (int i = 0; i < length; i++)
        {
            res.P[i] = P[i]/scalar;
        }
    }
    //---
    void divideByScalar(int scalar, std::shared_ptr<int[]> &res)
    {
        if(res == 0) 
        {
            res = std::shared_ptr<int[]>(new int[length]);
        }      
        for (int i = 0; i < length; i++)
        {
            res[i] = P[i]/scalar;
        }
    }
    //---
    static void divideByScalar( std::shared_ptr<int[]> &vec, int scalar, int N,  std::shared_ptr<int[]> & res)
    {
        if(res == 0) 
        {
            res = std::shared_ptr<int[]>(new int[N]);
        }
        for (int i = 0; i < N; i++)
        {
            res[i] = vec[i]/scalar;
        }
    }

    // //---------------------------------------------------
    // // MULTIPLY TRANSPOSE VECTOR_INT BY MATRIX
    // //---------------------------------------------------
    // VECTOR_INT operator & (MATRIX mat)
    // {
    //     int nRows = mat.nRow;
    //     if (length != nRows) throw_line("ERROR: VECTOR_INT and matrix sizes are not compatible \n")        
    //     int nCols = mat.nCol;
    //     VECTOR_INT res(nCols);
    //     int* resP = res.P;
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
    VECTOR_INT operator * (int scalar)
    {
        VECTOR_INT res(length);
        std::shared_ptr<int[]> resP = res.P;
        for (int i = 0; i < length; i++)
        {
            resP[i] = P[i] * scalar;
        }
        return res;
    }
    //---
    void multiplyByScalar(int scalar, VECTOR_INT &res)
    {
        for (int i = 0; i < length; i++)
        {
            res.P[i] = P[i]*scalar;
        }
    }
    //---
    void multiplyByScalar(int scalar,  std::shared_ptr<int[]> &res)
    {
        if(res == 0) 
        {
            res = std::shared_ptr<int[]>(new int[length]);
        }
        for (int i = 0; i < length; i++)
        {
            res[i] = P[i]*scalar;
        }
    }
    //---
    static void multiplyByScalar(std::shared_ptr<int[]> &vec, int scalar, int N, std::shared_ptr<int[]> &res)
    {
        if(res == 0) 
        {
            res = std::shared_ptr<int[]>(new int[N]);
        }
        for (int i = 0; i < N; i++)
        {
            res[i] = vec[i]*scalar;
        }
    }
    //---
    static void multiplyByScalar(std::shared_ptr<int[]> &vec, int scalar, int N, VECTOR_INT &res)
    {
        std::shared_ptr<int[]> P_in = res.P;
        for (int i = 0; i < N; i++)
        {
            P_in[i] = vec[i]*scalar;
        }
    }


    //------------------------------------------------------
    // VECTOR_INT NORM
    //------------------------------------------------------
    int norm()
    {
        int res = 0;
        for (int i = 0; i < length; i++)
        {
                res += P[i] * P[i];
        }
        res = sqrt(res);
        return res;
    }
    //---
    static int norm(VECTOR_INT &vec)
    {
        int res = 0;
        std::shared_ptr<int[]> P_in = vec.P;
        for (int i = 0; i < vec.length; i++)
        {
                res += P_in[i] * P_in[i];
        }
        res = sqrt(res);
        return res;
    }
    //
    static int norm(std::shared_ptr<int[]> &vec, int N)
    {
        int res = 0;
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
    int norm(int begin, int end)
    {
        if (begin < 0 || end > length) throw_line("ERROR: The indices are not valid \n");
        int res = 0;
        for (int i = begin; i < end; i++)
        {
                res += P[i] * P[i];
        }
        res = sqrt(res);
        return res;
    }

    //------------------------------------------------------
    // COPY VECTOR_INT INTO OTHER VECTOR_INT
    //------------------------------------------------------
    void operator = (VECTOR_INT &inVec) //(works also to copy vectors of different sizes)
    {
        //if (length != inVec.length) throw_line("ERROR: the vectors have different lengths \n;")
        if (length == 0)
        {
            length = inVec.length;
            P = std::shared_ptr<int[]>(new int[length]);
        }
        // length = inVec.length;
        // if (P != nullptr) free(P); 
        // P = (int*) malloc(length * sizeof(int));

        std::shared_ptr<int[]> P_in = inVec.P;
        for (int i = 0; i < length; i++)
        {
            P[i] = P_in[i];
        }
    }
    //---
    void operator = (int* &inVec)
    {
        for (int i = 0; i < length; i++)
        {
            P[i] = inVec[i];
        }
    }
    //---
    void operator = (std::shared_ptr<int[]> &inVec)
    {
        for (int i = 0; i < length; i++)
        {
            P[i] = inVec[i];
        }
    }
    //---
    void copyTo(VECTOR_INT &res)
    {
        std::shared_ptr<int[]> P_in = res.P;
        for (int i = 0; i < length; i++) P_in[i] = P[i];
    }
    //---
    void copyTo(std::shared_ptr<int[]> &res)
    {
        if(res == 0)
        {
            res = std::shared_ptr<int[]>(new int[length]);
        } 
        for (int i = 0; i < length; i++) res[i] = P[i];
    }
    //---
    void copy(std::shared_ptr<int[]> &vec)
    {
        if(vec == nullptr) throw_line("ERROR: COPYING NULL PTR\n");
        for (int i = 0; i < length; i++) P[i] = vec[i];
    }
    //---
    static void copy(std::shared_ptr<int[]> &inVec, std::shared_ptr<int[]> &outVec, int N, int startPos = 0)
    {
        for (int i = 0; i < N; i++) outVec[i+startPos] = inVec[i];
    }
    //---
    static void exact_copy(VECTOR_INT &in_vec, VECTOR_INT &out_vec)
    {
        out_vec.clear();
        out_vec.length = in_vec.length;
        out_vec.P = in_vec.P;
    }

    //----------------------------------------
    // GET MAX AND MIN
    //----------------------------------------
    int max()
    {
        int max = int(-1e50);
        for (int i = 0; i < length; i++)
        {
            if (P[i] > max)
            {
                max = P[i];
            }
        }
        return max;
    }

    int min()
    {
        int min = int(1e50);
        for (int i = 0; i < length; i++)
        {
            if (P[i] < min)
            {
                min = P[i];
            }
        }
        return min;
    }

    //----------------------------------------
    // APPEND VECTOR_INT
    //----------------------------------------
    void append(VECTOR_INT &vecIn)
    {
        int newL = length + vecIn.length;
        std::shared_ptr<int[]> newP(new int[newL]);

        for (int i = 0; i < length; i++) newP[i] = P[i];

        copy(vecIn.P, newP, vecIn.length, length);

        P.reset(); P = 0;
        P = newP;
        length = newL;
    }
    //---
    void append (std::shared_ptr<int[]> &vecIn, int vecInL)
    {
        int newL = length + vecInL;

        std::shared_ptr<int[]> newP(new int[newL]);

        for (int i = 0; i < length; i++) newP[i] = P[i];
        copy(vecIn, newP, vecInL, length);

        P.reset(); P = 0;
        P = newP;
        length = newL;
    }
    //---
    static void append (VECTOR_INT &vec, VECTOR_INT &vecIn)
    {
        vec.append(vecIn);
    }
    //---
    static void append (std::shared_ptr<int[]> &vec, int vecL, std::shared_ptr<int[]> &vecIn, int vecInL)
    {
        int newL = vecL + vecInL;
        std::shared_ptr<int[]> newP(new int[newL]);
        for (int i = 0; i < vecL; i++) newP[i] = vec[i];
        copy(vecIn, newP, vecInL, vecL);

        vec.reset(); vec = 0;
        vec = newP;
    }

    //-------------------------------
    // APPEND ONE VALUE
    //-------------------------------
    void append (int val)
    {
        int newL = length + 1;
        std::shared_ptr<int[]> newP(new int[newL]);
        for (int i = 0; i < length; i++) newP[i] = P[i];
        newP[length]= val;
        P.reset(); P = 0;
        P = newP;
        length = newL;
    }
    //---

    //-----------------------------------------------------
    // SHRINK VECTOR (keep only newL elements of the array)
    //-----------------------------------------------------
    void shrink(int newL)
    {
        std::shared_ptr<int[]> newP(new int[newL]);

        copy(P, newP, newL);
        
        P.reset(); P = 0;
        P = newP;
        length = newL;
    }
    //---
    static void shrink(std::shared_ptr<int[]> &vec, int newL)
    {
        std::shared_ptr<int[]> newP(new int[newL]);

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
        std::shared_ptr<int[]> newP(new int[newL]);
        copyTo(newP);
        
        P.reset(); P = 0;
        P = newP;
        length = newL;
    }
    //---
    
    //-------------------------------------------------------
    // UNIFY
    //-------------------------------------------------------
    //---
    static std::shared_ptr<int[]> merge (std::shared_ptr<int[]> &vec1, int vec1L, std::shared_ptr<int[]> &vec2, int vec2L)
    {
        int newL = vec1L + vec2L;

        std::shared_ptr<int[]> res(new int[newL]);

        copy(vec1, res, vec1L);
        copy(vec2, res , vec2L, vec1L);

        return res;
    }
    //------------------------------------------------------
    // CHECK IF THERE IS AN ELEMENT IN A VECTOR
    //------------------------------------------------------
    bool hasInOrd(int value, int &nElPassed)
    {
        bool answer = false;
        int i;
        for (i = 0; i < length; i++)
        {
            int tempVal = value - P[i];
            if (tempVal > 0) continue;

            if (tempVal == 0)
            {
                answer = true;
                nElPassed = i;
                return answer;
            } 
            if (tempVal < 0)
            {
                nElPassed = i;
                return false;
            } 
        }
        nElPassed = length;
        return answer;
    }
    //---
    bool hasIn(int value, int &index)
    {
        bool answer = false;
        index = -1;
        for (int i = 0; i < length; i++)
        {
            if (P[i] == value)
            {
                answer = true;
                index = i;
                return answer;
            } 
        }
        return answer;
    }
    //---
    bool hasIn(int value)
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
    bool isIn(int value, int* idAddress)
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
    static bool isIn(std::shared_ptr<int[]> &V_in, int N, int value)
    {
        bool answer = false;
        int i;
        for (i = 0; i < N; i++)
        {
            if (V_in[i] == value)
            {
                answer = true;
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
    //-------------------------------------------------------
    // PRINT VECTOR_INT
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
    void print(const char* outFileName, bool length_flag = false)
    {

        FILE* outFile = fopen(outFileName,"w");
        
        if (length_flag) fprintf(outFile, "%d\n", length);

        //Print matrix;
        for (int i = 0; i < length; i++)
        {
            fprintf(outFile, "%d\n", P[i]);
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
    

    //-----------------------------------------------
    // PRINT 2 POINTERS SIDE BY SIDE (4 int, int dispositions)
    //-----------------------------------------------
    static void print2(std::shared_ptr<int[]> &vec1, std::shared_ptr<int[]> &vec2, int N, std::string name1, std::string name2)
    {
        std::cout << "\n-------\t\t" << "-------\n";
        std::cout << name1 << "\t\t" << name2 << "\n";
        std::cout << "-------\t\t" << "-------\n";
        for (int i = 0; i < N; i++)
        {
            std::cout << " " << vec1[i] << "\t\t " << vec2[i] << "\n";
        }
        std::cout << "-------\t\t" << "-------\n";
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
};