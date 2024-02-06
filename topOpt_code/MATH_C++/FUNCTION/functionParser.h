// 1    SPLIT SU []
// 2    SPLIT SU + E -
// 3    SPLIT SU * E /
// 4    SPLIT SU ()
// RICERCA FUNZIONI SUL PEZZO ATTACCATO ALLE ()
// VALUTIAMO E TORNIAMO INDIETRO
#pragma once
#include "../../CODE_HEADERS/codeHeader.h"

class FUNCTION_PARSER
{
    public:
    //---
    static void evalCoord(std::string str, MATRIX &coord, prec t, VECTOR &value)
    {
        int nNodes = coord.nRow; 
        // FIRST STEP: SPLIT ON + - SIGNS
        value.initialize(nNodes);
        std::string signs; 
        std::vector<std::string> resSplit;

        std::string delimiters = "+-";

        splitDelimiters(&str[0], str.length(), &delimiters[0], delimiters.length(), resSplit, signs);
        //-------------------
        delimiters = "*/";
        int size = resSplit.size();
        std::vector<VECTOR> prodValues(size);
        int i_start = 0;
        if (resSplit[0].length() == 0) 
        {
            prodValues[0].setZeros(nNodes); i_start++;
        }

        for (int i = i_start; i < size; i++)
        {
            VECTOR tempValue(nNodes);
            std::string tempString = resSplit[i];
            std::vector<std::string> tempResSplit;
            std::string tempSigns; 
            splitDelimiters(&tempString[0], tempString.length(), &delimiters[0], delimiters.length(), tempResSplit, tempSigns);

            int locSize = tempResSplit.size();
            if (locSize == 1)
            {
                tempString = tempResSplit[0];
                tempString = removeSquareBracket(tempString);
                eval(tempString, coord, t, tempValue);
                prodValues[i] = tempValue;
                continue;
            }
            tempString = tempResSplit[0];
            tempString = removeSquareBracket(tempString);
            evalCoord(tempString, coord, t, tempValue);
            for (int j = 1; j < locSize; j++)
            {
                tempString = tempResSplit[j];
                tempString = removeSquareBracket(tempString);

                char tempSign = tempSigns[j-1];
                VECTOR tVal;
                if (tempSign == '*') 
                {
                    evalCoord(tempString, coord, t, tVal);
                    tempValue *= tVal;
                }
                else if (tempSign == '/')
                {
                    evalCoord(tempString, coord, t, tVal);
                    tempValue /= tVal;
                } 
                else throw_line("ERROR: wrong operators in input function\n");
            }
            prodValues[i] = tempValue;
        }
        //--------
        value = prodValues[0];
        for (int i = 1; i < size; i++) 
        {
            char tempSign = signs[i-1];
            if (tempSign == '+') value  += prodValues[i];
            else if (tempSign == '-') value -= prodValues[i];
            else throw_line("ERROR: wrong operators in input function\n");
            prodValues[i].P.reset();;
        }
    }
    //---
    static prec evalScalar(std::string str, prec x, prec y, prec z, prec t)
    {
        // FIRST STEP: SPLIT ON + - SIGNS
        prec value = 0;
        std::string signs; 
        std::vector<std::string> resSplit;

        std::string delimiters = "+-";

        splitDelimiters(&str[0], str.length(), &delimiters[0], delimiters.length(), resSplit, signs);
        //-------------------
        delimiters = "*/";
        int size = resSplit.size();
        std::vector<prec> prodValues(size);
        int i_start = 0;
        if (resSplit[0].length() == 0) 
        {
            prodValues[0] = 0; i_start++;
        }
        for (int i = i_start; i < size; i++)
        {
            prec tempValue = 0;
            std::string tempString = resSplit[i];
            std::vector<std::string> tempResSplit;
            std::string tempSigns; 
            splitDelimiters(&tempString[0], tempString.length(), &delimiters[0], delimiters.length(), tempResSplit, tempSigns);

            int locSize = tempResSplit.size();
            if (locSize == 1)
            {
                tempString = tempResSplit[0];
                tempString = removeSquareBracket(tempString);
                tempValue = eval(tempString, x, y, z, t);
                prodValues[i] = (tempValue);
                continue;
            }
            tempString = tempResSplit[0];
            tempString = removeSquareBracket(tempString);
            tempValue = evalScalar(tempString, x, y, z, t);
            for (int j = 1; j < locSize; j++)
            {
                tempString = tempResSplit[j];
                tempString = removeSquareBracket(tempString);

                char tempSign = tempSigns[j-1];
                if (tempSign == '*') tempValue*= evalScalar(tempString, x, y, z, t);
                else if (tempSign == '/') tempValue/= evalScalar(tempString, x, y, z, t);
                else throw_line("ERROR: wrong operators in input function\n");
            }
            prodValues[i] = tempValue;
        }
        //--------
        value = prodValues[0];
        for (int i = 1; i < size; i++) 
        {
            prec tempValue = prodValues[i];
            char tempSign = signs[i-1];
            if (tempSign == '+') value  += tempValue;
            else if (tempSign == '-') value -= tempValue;
            else throw_line("ERROR: wrong operators in input function\n");
        }
        return value;
    }
    private:
    //----------------------------------------------------
    // EVAL
    //----------------------------------------------------
    static void eval(std::string str, MATRIX &coord, prec t, VECTOR & value)
    {
        int nNodes = coord.nRow; int dim = coord.nCol;
        value.initialize(nNodes);
        switch (dim)
        {
        case 2:
        {
            for (int i = 0; i < nNodes; i++)
            {
                prec* point = coord[i];
                prec x = point[0]; prec y = point[1];
                // parse
                if (str == "x")
                {
                    value[i] = x;
                }else if (str == "y")
                {
                    value[i] = y;
                }
                else if (str == "z") {throw_line("ERROR: evaluating 2D function in z-coord\n");}
                else if (str == "t")
                {
                    value[i] = t;
                }
                else if (str == "pi") value[i] = 3.141592653589793238462643383279502884;
                else if (std::isdigit(str[0])) value[i] = std::stod(str);
                else 
                {
                    // search for round brackets
                    std::string funcStr; std::string argStr;
                    bool found = manageRoundBrackets(str, funcStr, argStr);
                    if (!found) 
                    {
                        std::string msg = "ERROR: wrong data in input function, found symbol: " + str + "\n";
                        throw_line(&msg[0]);
                    }
                    prec argVal = evalScalar(argStr, x, y, 0, t);
                    if (funcStr == "cos") value[i] = cos(argVal);
                    else if (funcStr == "sin") value[i] = sin(argVal);
                    else if (funcStr == "log") value[i] = log(argVal);
                    else if (funcStr == "exp") value[i] = exp(argVal);
                    else if (funcStr == "sqrt") value[i] = sqrt(argVal);
                    else 
                    {
                        std::string msg = "ERROR: Evaluating a function not defined in the database " + funcStr + "\n";
                        throw_line(msg);
                    }
                }
            }
            break;
        } 
        case 3:
        {
            for (int i = 0; i < nNodes; i++)
            {
                prec* point = coord[i];
                prec x = point[0]; prec y = point[1]; prec z = point[2];
                // parse
                if (str == "x")
                {
                    value[i] = x;
                } 
                else if (str == "y")
                {
                    value[i] = y;
                }
                else if (str == "z")
                {
                    value[i] = z;
                }
                else if (str == "t")
                {
                    value[i] = t;
                }
                else if (str == "pi") value[i] = 3.141592653589793238462643383279502884;
                else if (std::isdigit(str[0])) value[i] = std::stod(str);
                else 
                {
                    // search for round brackets
                    std::string funcStr; std::string argStr;
                    bool found = manageRoundBrackets(str, funcStr, argStr);
                    if (!found) throw_line("ERROR: wrong data in input function\n");
                    prec argVal = evalScalar(argStr, x, y, z, t);
                    if (funcStr == "cos") value[i] = cos(argVal);
                    else if (funcStr == "sin") value[i] = sin(argVal);
                    else if (funcStr == "log") value[i] = log(argVal);
                    else if (funcStr == "exp") value[i] = exp(argVal);
                    else if (funcStr == "sqrt") value[i] = sqrt(argVal);
                    else 
                    {
                        throw_line("ERROR: Evaluating a function not defined in the database\n");
                    }
                }
            }
            break;
        }
        }
    }
    //-----------------
    static prec eval(std::string str, prec x, prec y, prec z, prec t)
    {
        prec value; 
        if (str == "x")
        {
            value = x;
        } 
        else if (str == "y")
        {
            value = y;
        }  
        else if (str == "z")
        {
            value = z;
        }  
        else if (str == "t")
        {
            value = t;
        }  
        else if (str == "pi") value = 3.141592653589793238462643383279502884;
        else if (std::isdigit(str[0])) value = std::stod(str);
        else 
        {
            // search for round brackets
            std::string funcStr; std::string argStr;
            bool found = manageRoundBrackets(str, funcStr, argStr);
            if (!found) throw_line("ERROR: wrong data in input function\n");
            prec argVal = evalScalar(argStr, x, y, z, t);
            if (funcStr == "cos") value = cos(argVal);
            else if (funcStr == "sin") value = sin(argVal);
            else if (funcStr == "log") value = log(argVal);
            else if (funcStr == "exp") value = exp(argVal);
            else 
            {
                throw_line("ERROR: Evaluating a function not defined in the database\n");
                return 0;
            }

        }
        return value;
    }
    //-----------------------------------------------------------------------------------
    // MANAGE BRACKETS
    //-----------------------------------------------------------------------------------
    static std::string removeSquareBracket(std::string str)
    {
        int length = str.length(); int i = 0;
        if (str[0] == '[')
        {
            std::string outStr;
            int count = 1;
            while (count != 0 && i < length)
            {
                i++;
                if (str[i] == ']') count--;
                else if (str[i] ==  '[') count++;

                if (count != 0) outStr += str[i];
            }
            return outStr;
        }
        else return str;
    }
    //---
    static bool manageRoundBrackets(std::string str, std::string &funcStr, std::string &argStr)
    {
        int length = str.length(); int i = 0; bool found = false;
        while (i < length)
        {
            char tempChar = str[i];
            if (tempChar == '(')
            {
                int count = 1;
                while (i < length)
                {
                    i++;
                    if (str[i] == '(') count++;
                    else if (str[i] == ')') count--;
                    if (count == 0) break;
                    argStr += str[i];
                }
                i++;
                found = true;
                // if (i != length) throw_line("error in function evaluation with round brackets\n");
            } else
            {
                funcStr += tempChar;
                i++;
            }
        }
        return found;
    }
    //------------------------------------------------------------------------
    // GENERAL SPLIT
    //------------------------------------------------------------------------
    //---
    static bool splitDelimiters(char* str, int length, char* delimiters, int nDelimiters, std::vector<std::string> &resSplit, std::string &signsSplit)
    {
        int countSBrk = 0; int countRBrk = 0; int i  = 0;
        bool found = false;
        std::string tempString;
        while (i < length)
        {
            char tempChar = str[i];
            if (str[i] == '[') 
            {
                countSBrk++;
                tempString += str[i];
                while (!(countSBrk == 0 && countRBrk == 0) && i < length)
                {
                    i++;
                    tempString += str[i];
                    if (str[i] == '[') countSBrk++;
                    else if (str[i] == ']') countSBrk--;
                    else if (str[i] == '(') countRBrk++;
                    else if (str[i] == ')') countRBrk--;
                }
            } else if (str[i] == '(') 
            {
                countRBrk++;
                tempString += str[i];
                while (!(countSBrk == 0 && countRBrk == 0) && i < length)
                {
                    i++;
                    tempString += str[i];
                    if (str[i] == '[') countSBrk++;
                    else if (str[i] == ']') countSBrk--;
                    else if (str[i] == '(') countRBrk++;
                    else if (str[i] == ')') countRBrk--;
                }
            } else if (isIn(delimiters, nDelimiters, tempChar))
            {
                signsSplit += str[i];
                resSplit.push_back(tempString);
                tempString = "";
                found = true;
            }
            else tempString += tempChar;
            i++;
        }
        resSplit.push_back(tempString);
        return found;
    }
    //---
    static bool isIn(char* str, int length, char val)
    {
        for (int i = 0; i < length; i++)
        {
            if (str[i] == val) return true;
        }
        return false;
    }

};

class PIECEWISE_FUNCTION
{
    public:
    int dim;
    int nTimeCases;
    VECTOR startTimeCases;
    std::vector<std::vector<std::string>> function;

    PIECEWISE_FUNCTION(){};

    PIECEWISE_FUNCTION(int dimension, int nCases)
    {
        dim = dimension;
        nTimeCases = nCases;
        startTimeCases.initialize(nCases);
        function.resize(nCases);
        for (int icase = 0; icase < nCases; icase++)
        {
            function[icase].resize(dim);
        }
    }

    //--------------------------------
    void initialize(int dimension, int nCases)
    {
        dim = dimension;
        nTimeCases = nCases;
        startTimeCases.initialize(nCases);
        function.resize(nCases);
        for (int icase = 0; icase < nCases; icase++)
        {
            function[icase].resize(dim);
        }
    }
    //---
    void define(int dimension, int nCases, VECTOR startTimes, std::vector<std::vector<std::string>> funcsIn)
    {
        dim = dimension;
        nTimeCases = nCases;
        startTimeCases.initialize(nCases);
        startTimeCases = startTimes;
        function = funcsIn;
    }
    //-------------------------------
    prec* operator () (prec x, prec y, prec z, prec t)
    {
        prec* res = (prec*) malloc(dim * sizeof(prec));
        int timeCase = -1;
        for (int icase = 0; icase < nTimeCases-1; icase++)
        {
            if (t < startTimeCases[icase+1])
            {
                timeCase = icase;
                break;
            }
        }  
        if (timeCase == -1) timeCase = nTimeCases-1;
        for (int icomp = 0; icomp < dim; icomp++)
        {
            std::string currFunc = function[timeCase][icomp];
            res[icomp] = FUNCTION_PARSER::evalScalar(currFunc, x, y, z, t);
        }  
        return res;  
    }
    //-------------------------------
    std::vector<VECTOR> operator () (MATRIX &coord, prec t)
    {
        std::vector<VECTOR> res(dim);
        int timeCase = -1;
        for (int icase = 0; icase < nTimeCases-1; icase++)
        {
            if (t < startTimeCases[icase+1])
            {
                timeCase = icase;
                break;
            }
        }  
        if (timeCase == -1) timeCase = nTimeCases-1; 
        for (int icomp = 0; icomp < dim; icomp++)
        {
            std::string currFunc = function[timeCase][icomp];
            VECTOR tVal;
            FUNCTION_PARSER::evalCoord(currFunc, coord, t, tVal);
            res[icomp] = tVal;
        }   
        return res;  
    }
    //-----------------------------------
    std::vector<VECTOR> eval(MATRIX &coord, prec t)
    {
        return (*this)(coord, t);
    }
    //-----------------------------------
    // PRINT
    //-----------------------------------
    void print()
    {
        std::cout << "\n........" << "Piecewise Function" << "........\n";
        std::cout << "   --\t"; std::cout << "\n";
        std::cout << "  /\t"; std::cout << "\n";
        for (int icase = 0; icase < nTimeCases-1; icase++)
        {
            std::cout << " {\t"; 
            CUSTOM::printRowStdInLine(function[icase]);
            std::cout << "\t for t in [" << startTimeCases[icase] << " , " << startTimeCases[icase+1] << "[ \n";
        }
        std::cout << " {\t"; 
        CUSTOM::printRowStdInLine(function[nTimeCases-1]);
        std::cout << "\t for t in [" << startTimeCases[nTimeCases-1] << " , " << "+oo[ \n";
        std::cout << "  \\ \t"; std::cout << "\n";
        std::cout << "    --\t"; std::cout << "\n";
        printf("\n.........................\n");
    }
    //---
    void print(std::string funcName)
    {
        std::cout << "\n==============| " << funcName << " |==============\n";
        std::cout << "  /--\t"; std::cout << "\n";
        for (int icase = 0; icase < nTimeCases-1; icase++)
        {
            std::cout << " {\t"; 
            CUSTOM::printRowStdInLine(function[icase]);
            std::cout << "\t for t in [" << startTimeCases[icase] << " , " << startTimeCases[icase+1] << "[ \n";
        }
        std::cout << " {\t"; 
        CUSTOM::printRowStdInLine(function[nTimeCases-1]);
        std::cout << "     for t in [" << startTimeCases[nTimeCases-1] << "," << "+oo[ \n";
        std::cout << "  \\-- \t"; std::cout << "\n";
        printf("========================================================\n\n"); 
    }
};

class PIECEWISE_MATRIX
{
    public:
    int dim;
    int nTimeCases;
    VECTOR startTimeCases;
    MATRIX constant;
    std::vector<std::string> function;

    PIECEWISE_MATRIX(){};

    PIECEWISE_MATRIX(int dimension, int nCases)
    {
        dim = dimension;
        nTimeCases = nCases;
    }
    //--------------------------------
    void initialize( int nCases, int dimension)
    {
        dim = dimension;
        nTimeCases = nCases;
        startTimeCases.initialize(nCases); 
        constant.initialize(nCases, dim);
        function.resize(nCases);
    }
    //---
    void define(int dimension, int nCases, VECTOR startTimes, MATRIX constants, std::vector<std::string> funcsIn)
    {
        dim = dimension;
        nTimeCases = nCases;
        startTimeCases.initialize(nCases); 
        startTimeCases = startTimes;
        constant = constants;
        function = funcsIn;
    }
    //-------------------------------
    prec* operator () (prec t)
    {
        prec* res = (prec*) malloc(dim * sizeof(prec));
        int timeCase = -1;
        for (int icase = 0; icase < nTimeCases-1; icase++)
        {
            if (t < startTimeCases[icase+1])
            {
                timeCase = icase;
                break;
            }
        }  
        if (timeCase == -1) timeCase = nTimeCases-1;
        prec* constantCase = constant[timeCase];
        for (int icomp = 0; icomp < dim; icomp++) res[icomp] = constantCase[icomp] * FUNCTION_PARSER::evalScalar(function[timeCase], 0, 0, 0, t); 
        return res;  
    }

    //-----------------------------------
    // PRINT
    //-----------------------------------
    void print()
    {
        std::cout << "\no-o-o-o-o-o-o-| Piecewise Matrix |-o-o-o-o-o-o-o\n";
        //std::cout << "   --\t"; std::cout << "\n";
        std::cout << "  /--\t"; std::cout << "\n";
        for (int icase = 0; icase < nTimeCases-1; icase++)
        {
            std::cout << " {\t"; 
            prec* tempP = constant[icase];
            VECTOR::printRowMatlabInLine(tempP, dim);
            std::cout << " * (" << function[icase] << ")";
            std::cout << "\t for t in [" << startTimeCases[icase] << " , " << startTimeCases[icase+1] << "[ \n";
        }
        std::cout << " {\t"; 
        prec* tempP = constant[nTimeCases-1];
        VECTOR::printRowMatlabInLine(tempP, dim);
        std::cout << " * (" << function[nTimeCases-1] << ")";
        std::cout << "\t for t in [" << startTimeCases[nTimeCases-1] << " , " << "+oo[ \n";
        std::cout << "  \\-- \t"; std::cout << "\n";
        //std::cout << "    --\t"; std::cout << "\n";
        printf("o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o\n");
    }
    //---
    void print(std::string funcName)
    {                   
        std::cout << "\no-o-o-o-o-o-o-| " << funcName << " |-o-o-o-o-o-o-o\n";
        //std::cout << "   --\t"; std::cout << "\n";
        std::cout << "  /--\t"; std::cout << "\n";
        for (int icase = 0; icase < nTimeCases-1; icase++)
        {
            std::cout << " {\t"; 
            prec* tempP = constant[icase];
            VECTOR::printRowMatlabInLine(tempP, dim);
            std::cout << " * (" << function[icase] << ")";
            std::cout << "\t for t in [" << startTimeCases[icase] << " , " << startTimeCases[icase+1] << "[ \n";
        }
        std::cout << " {\t"; 
        prec* tempP = constant[nTimeCases-1];
        VECTOR::printRowMatlabInLine(tempP, dim);
        std::cout << " * (" << function[nTimeCases-1] << ")";
        std::cout << "\t for t in [" << startTimeCases[nTimeCases-1] << " , " << "+oo[ \n";
        std::cout << "  \\-- \t"; std::cout << "\n";
        //std::cout << "    --\t"; std::cout << "\n";
        printf("o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o\n");
    }
    //---
};