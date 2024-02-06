#pragma once
#include "codeHeader.h"

class STREAM
{
public:
    //-------------------
    // GET LINES
    //-------------------
    static void getLines(std::ifstream &stream, std::string line, int nLines){
        // if the next value to read is after #n lines, the use this with nLines = n-1
        for (int i = 0; i < nLines; i++){
            getline(stream, line);
        }
    }

    //-------------------
    // GET FUNCTION
    //-------------------
    static void getFunctionCoef(std::ifstream &stream, std::string line, std::istringstream &iss, prec* &coeff, int &nFuncCoef){
        // read coefficients for the following function form
        //    f =   A    + 
        //          B  x + C  x^2  + D cos(x) + E sin(x) +
        //          F  y + G  y^2  + H cos(y) + I  sin(y) + J  x*y +
        //          K  z + L  z^2  + M cos(z)  + N sin(z)  + O  x*z + P y*z
        getLines(stream, line, 2);
        getline(stream, line);
        iss.str(line);
        for (int i = 0; i < nFuncCoef; i++){
            iss >> coeff[i];
        }
        iss.clear();
    }
    //---
    // static void getFunction(std::ifstream &stream, std::string line, std::istringstream &iss, int dim, FUNCTION &func)
    // {
    //     func.initialize(dim);
    //     getFunctionCoef(stream, line, iss, func.coef, func.nFuncCoef);
    // }
    
    //-------------------
    // GET VALUE
    //-------------------
    static void getValue(std::ifstream &stream, std::string line, std::istringstream &iss, bool &val){
        // get a "prec" value using a pointer
        getline(stream, line);
        iss.str(line);
        iss >> val;
        iss.clear();
    }
    //---
    static void getValue(std::ifstream &stream, std::string line, std::istringstream &iss, prec &val){
        // get a "prec" value using a pointer
        getline(stream, line);
        iss.str(line);
        iss >> val;
        iss.clear();
    }
    //---
    static void getValue(std::ifstream &stream, std::string line, std::istringstream &iss, int &val){
        // get an "int" value using a pointer
        getline(stream, line);
        iss.str(line);
        iss >> val;
        iss.clear();
    }
    //---
    static void getValue(std::ifstream &stream, std::string line, std::istringstream &iss, std::string &val){
        // get a "prec" value using a pointer
        getline(stream, line);
        iss.str(line);
        iss >> val;
        iss.clear();
    }
    //---

    //-------------------
    // GET ROW VECTOR
    //-------------------
    static void getRowVector(std::ifstream &stream, std::string line, std::istringstream &iss, prec* &vect, int vLength){
        // get a "row" vector of "prec" of length "vLength" using a pointer 
        getline(stream, line);
        iss.str(line);
        for (int i = 0; i < vLength; i++){
        iss >> vect[i]; 
        }
        iss.clear();
    }
    //---

    static void getRowVector(std::ifstream &stream, std::string line, std::istringstream &iss, int* &vect, int vLength){
        // get a "row" vector of "int" of length "vLength" using a pointer 
        getline(stream, line);
        iss.str(line);
        for (int i = 0; i < vLength; i++){
        iss >> vect[i]; 
        }
        iss.clear();
    }
    //---

    static void getRowVector(std::ifstream &stream, std::string line, std::istringstream &iss, VECTOR &vec){
        // get a "row" vector of "int" of length "vLength" using a pointer 
        getline(stream, line);
        iss.str(line);
        int vLength = vec.length;
        for (int i = 0; i < vLength; i++){
        iss >> vec[i]; 
        }
        iss.clear();
    }
    //---

    static void getRowVector(std::ifstream &stream, std::string line, std::istringstream &iss, VECTOR_INT &vec)
    {
        // get a "row" vector of "int" of length "vLength" using a pointer 
        getline(stream, line);
        iss.str(line);
        int vLength = vec.length;
        for (int i = 0; i < vLength; i++){
        iss >> vec[i]; 
        }
        iss.clear();
    }
    //---

    static void getRowVector(std::ifstream &stream, std::string line, std::istringstream &iss, std::vector<std::string> &vec){
        // get a "row" vector of "int" of length "vLength" using a pointer 
        getline(stream, line);
        iss.str(line);
        for (unsigned int i = 0; i < vec.size(); i++){
        iss >> vec[i]; 
        }
        iss.clear();
    }

    //-------------------
    // GET COLUMN VECTOR
    //-------------------
    static void getColVector(std::ifstream &stream, std::string line, std::istringstream &iss, prec* &vect, int vLength){
        // get a "column" vector of "prec" of length "vLength" using a pointer     
        for (int i = 0; i < vLength; i++)
        {
            getValue(stream, line, iss, vect[i]);
        }
        iss.clear();
    }
    //---
    static void getColVector(std::ifstream &stream, std::string line, std::istringstream &iss, VECTOR &vect, int vLength){
        // get a "column" vector of "prec" of length "vLength" using a pointer     
        for (int i = 0; i < vLength; i++)
        {
            getValue(stream, line, iss, vect[i]);
        }
        iss.clear();
    }
    //---

    //---
    static void getColVector(std::ifstream &stream, std::string line, std::istringstream &iss, VECTOR_INT &vect, int vLength){
        // get a "column" vector of "int" of length "vLength" using a pointer     
        for (int i = 0; i < vLength; i++)
        {
            getValue(stream, line, iss, vect[i]);
        }
        iss.clear();
    }
    //---

    static void getColVector(std::ifstream &stream, std::string line, std::istringstream &iss, int* &vect, int vLength){    
        // get a "column" vector of "int" of length "vLength" using a pointer 
        for (int i = 0; i < vLength; i++)
        {
            getValue(stream, line, iss, vect[i]); 
        }
        iss.clear();
    }
    //---
    //----------------------
    // GET MATRIX
    //----------------------    
    static void getMatrix(std::ifstream &stream, std::string line, std::istringstream &iss, MATRIX_INT &mat)
    {
        try
        {
            int nRow = mat.nRow; int nCol = mat.nCol;
            std::shared_ptr<int[]> p = mat.P; int count = 0;
            for (int i = 0; i < nRow; i++)
            {
                getline(stream, line);
                iss.str(line);
                for (int j = 0; j < nCol; j++)
                {
                    iss >> p[count]; count++;
                }
                iss.clear();
            }
        }
        catch(const std::exception& e)
        {
            std::cerr << "ERROR IN STREAM: GET MATRIX"<< '\n';
        }
    }
    //---
    static void getMatrix(std::ifstream &stream, std::string line, std::istringstream &iss, MATRIX &mat)
    {
        try
        {
            int nRow = mat.nRow; int nCol = mat.nCol;
            std::shared_ptr<prec[]> p = mat.P; int count = 0;
            for (int i = 0; i < nRow; i++)
            {
                getline(stream, line);
                iss.str(line);
                for (int j = 0; j < nCol; j++)
                {
                    iss >> p[count]; count++;
                }
                iss.clear();
            }
        }
        catch(const std::exception& e)
        {
            std::cerr << "ERROR IN STREAM: GET MATRIX"<< '\n';
        }
    }

    //--------------------------------------------------
    // GET ALL INFO ABOUT A TIME DEPENDENT DIRICHLET BC
    //--------------------------------------------------
    static void getDirTime(std::ifstream &stream, std::string line, std::istringstream &iss, int dim, int &boundId, int &nCases, PIECEWISE_FUNCTION &func)
    {
        getline(stream, line);
        getValue(stream, line, iss, boundId);
        getline(stream, line);
        getValue(stream, line, iss, nCases);
        getline(stream, line);
        VECTOR startTimes(nCases);
        getColVector(stream, line, iss, startTimes, nCases);
        getline(stream, line);
        std::vector<std::vector<std::string>> funcs;
        funcs.resize(nCases);
        for (int icase = 0; icase < nCases; icase++)
        {
            funcs[icase].resize(dim);
            getRowVector(stream, line, iss, funcs[icase]);
        }
        func.define(dim, nCases, startTimes, funcs);
    }

    // static void getDirTimeCase(std::ifstream &stream, std::string line, std::istringstream &iss, int dim, prec &startTime, SPACE_TIME_FUNCTION &h)
    // {
    //     int nCoef1D = 5;
    //     int nFuncCoef = 1 + 4*(dim) + (dim*(dim-1))/2;
    //     getline(stream, line);
    //     getValue(stream, line, iss, startTime);
    //     getline(stream, line);
    //     prec* spaceCoef = (prec*) malloc(nFuncCoef * sizeof(prec));
    //     prec* timeCoef  = (prec*) malloc(nCoef1D * sizeof(prec));
    //     getRowVector(stream, line, iss, spaceCoef, nFuncCoef);
    //     getline(stream, line);
    //     getRowVector(stream, line, iss, timeCoef, nCoef1D);
    //     FUNCTION spaceFunc(dim, spaceCoef);
    //     FUNCTION timeFunc(1, timeCoef);
    //     h.define(spaceFunc, timeFunc);

    //     free(spaceCoef); free(timeCoef);
    // }
    //---

    //--------------------------------------------------
    // GET ALL INFO ABOUT A TIME DEPENDENT NEUMANN BC
    //--------------------------------------------------
    static void getNeuTime(std::ifstream &stream, std::string line, std::istringstream &iss, int &boundId, int &nCases, PIECEWISE_FUNCTION &matFunc)
    {
        getline(stream, line);
        getValue(stream, line, iss, boundId);
        getline(stream, line);
        getValue(stream, line, iss, nCases);
        getline(stream, line);
        matFunc.initialize(1, nCases);
        getColVector(stream, line, iss, matFunc.startTimeCases, nCases);
        getline(stream, line);
        for (int icase = 0; icase < nCases; icase++) 
        {
            getValue(stream, line, iss, matFunc.function[icase][0]); // [0] cause we're dealing with a 1d function of time  
        }
    }

    // static void getNeuTimeCase(std::ifstream &stream, std::string line, std::istringstream &iss, int dim, prec &startTime, prec* neuVec, FUNCTION &timeFunc)
    // {
    //     int nCoef1D = 5;
    //     getline(stream, line);
    //     getValue(stream, line, iss, startTime);
    //     getline(stream, line);
    //     getRowVector(stream, line, iss, neuVec, dim);
    //     getline(stream, line);
    //     getRowVector(stream, line, iss, timeFunc.coef, nCoef1D);
    // }
};