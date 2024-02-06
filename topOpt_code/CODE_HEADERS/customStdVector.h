#pragma once

#include "codeHeader.h"

class CUSTOM
{
public:
    // INT vector
    static void print(std::vector<int> vec)
    {
        std::cout << "\n-------\n";
        for (unsigned int i = 0; i < vec.size(); i++)
        {
            std::cout << " " << vec[i] << "\n";
        }
        std::cout << "-------\n";
    }
    //---
    static void print(std::vector<int> vec, std::string vecName)
    {
        std::cout << "\n-------\n";
        std::cout << vecName;
        std::cout << "\n-------\n";
        for (unsigned int i = 0; i < vec.size(); i++)
        {
            std::cout << " " << vec[i] << "\n";
        }
        std::cout << "-------\n";
    }
    //---
    static void printRow(std::vector<int> vec)
    {
        for (unsigned int i = 0; i < vec.size(); i++)
        {
            std::cout << vec[i] << " ";
        }
        std::cout << "\n";
    }
    //---
    static void printRow(std::vector<int> vec, std::string nameVec)
    {
        std::cout << nameVec << ": ";
        for (unsigned int i = 0; i < vec.size(); i++)
        {
            std::cout << vec[i] << " ";
        }
        std::cout << "\n";
    }
    //---
    static void printRowMatlab(std::vector<int> vec)
    {
        if (vec.size() == 0) std::cout << "[]\n";
        else
        {
            std::cout << "[";
            for (unsigned int i = 0; i < vec.size()-1; i++)
            {
                std::cout << vec[i] << ", ";
            }
            std::cout << vec[vec.size()-1] << "]\n";
        }
    }
    //---
    static void printRowMatlab(std::vector<int> vec, std::string nameVec)
    {
        std::cout << nameVec << ": ";
        if (vec.size() == 0) std::cout << "[]\n";
        else
        {
            std::cout << "[";
            for (unsigned int i = 0; i < vec.size()-1; i++)
            {
                std::cout << vec[i] << ", ";
            }
            std::cout << vec[vec.size()-1] << "]\n";
        }
    }
    //---
    static void printRowStd(std::vector<int> vec)
    {
        if (vec.size() == 0) std::cout << "{}\n";
        else
        {
            std::cout << "{";
            for (unsigned int i = 0; i < vec.size()-1; i++)
            {
                std::cout << vec[i] << ", ";
            }
            std::cout << vec[vec.size()-1] << "}\n";
        }
    }
    //---
    static void printRowStd(std::vector<int> vec, std::string nameVec)
    {
        std::cout << nameVec << ": ";
        if (vec.size() == 0) std::cout << "{}\n";
        else
        {
            std::cout << "{";
            for (unsigned int i = 0; i < vec.size()-1; i++)
            {
                std::cout << vec[i] << ", ";
            }
            std::cout << vec[vec.size()-1] << "}\n";
        }
    }
    //---
    static void printRowStdInLine(std::vector<int> vec)
    {
        if (vec.size() == 0) std::cout << "{}";
        else
        {
            std::cout << "{";
            for (unsigned int i = 0; i < vec.size()-1; i++)
            {
                std::cout << vec[i] << ", ";
            }
            std::cout << vec[vec.size()-1] << "}";
        }
    }
    //---
    static void printRowStdInLine(std::vector<int> vec, std::string nameVec)
    {
        std::cout << nameVec << ": ";
        if (vec.size() == 0) std::cout << "{}";
        else
        {
            std::cout << "{";
            for (unsigned int i = 0; i < vec.size()-1; i++)
            {
                std::cout << vec[i] << ", ";
            }
            std::cout << vec[vec.size()-1] << "}";
        }
    }
    //---

    // PREC vector
    static void print(std::vector<prec> vec)
    {
        std::cout << "\n-------\n";
        for (unsigned int i = 0; i < vec.size(); i++)
        {
            std::cout << " " << vec[i] << "\n";
        }
        std::cout << "-------\n";
    }
    //---
    static void print(std::vector<prec> vec, std::string vecName)
    {
        std::cout << "\n-------\n";
        std::cout << vecName;
        std::cout << "\n-------\n";
        for (unsigned int i = 0; i < vec.size(); i++)
        {
            std::cout << " " << vec[i] << "\n";
        }
        std::cout << "-------\n";
    }
    //---
    static void printRow(std::vector<prec> vec)
    {
        for (unsigned int i = 0; i < vec.size(); i++)
        {
            std::cout << vec[i] << " ";
        }
        std::cout << "\n";
    }
    //---
    static void printRow(std::vector<prec> vec, std::string nameVec)
    {
        std::cout << nameVec << ": ";
        for (unsigned int i = 0; i < vec.size(); i++)
        {
            std::cout << vec[i] << " ";
        }
        std::cout << "\n";
    }
    //---
    static void printRowMatlab(std::vector<prec> vec)
    {
        if (vec.size() == 0) std::cout << "[]\n";
        else
        {
            std::cout << "[";
            for (unsigned int i = 0; i < vec.size()-1; i++)
            {
                std::cout << vec[i] << ", ";
            }
            std::cout << vec[vec.size()-1] << "]\n";
        }
    }
    //---
    static void printRowMatlab(std::vector<prec> vec, std::string nameVec)
    {
        std::cout << nameVec << ": ";
        if (vec.size() == 0) std::cout << "[]\n";
        else
        {
            std::cout << "[";
            for (unsigned int i = 0; i < vec.size()-1; i++)
            {
                std::cout << vec[i] << ", ";
            }
            std::cout << vec[vec.size()-1] << "]\n";
        }
    }
    //---
    static void printRowStd(std::vector<prec> vec)
    {
        if (vec.size() == 0) std::cout << "{}\n";
        else
        {
            std::cout << "{";
            for (unsigned int i = 0; i < vec.size()-1; i++)
            {
                std::cout << vec[i] << ", ";
            }
            std::cout << vec[vec.size()-1] << "}\n";
        }
    }
    //---
    static void printRowStd(std::vector<prec> vec, std::string nameVec)
    {
        std::cout << nameVec << ": ";
        if (vec.size() == 0) std::cout << "{}\n";
        else
        {
            std::cout << "{";
            for (unsigned int i = 0; i < vec.size()-1; i++)
            {
                std::cout << vec[i] << ", ";
            }
            std::cout << vec[vec.size()-1] << "}\n";
        }
    }
    //---
    static void printRowStdInLine(std::vector<prec> vec)
    {
        if (vec.size() == 0) std::cout << "{}";
        else
        {
            std::cout << "{";
            for (unsigned int i = 0; i < vec.size()-1; i++)
            {
                std::cout << vec[i] << ", ";
            }
            std::cout << vec[vec.size()-1] << "}";
        }
    }
    //---
    static void printRowStdInLine(std::vector<prec> vec, std::string nameVec)
    {
        std::cout << nameVec << ": ";
        if (vec.size() == 0) std::cout << "{}";
        else
        {
            std::cout << "{";
            for (unsigned int i = 0; i < vec.size()-1; i++)
            {
                std::cout << vec[i] << ", ";
            }
            std::cout << vec[vec.size()-1] << "}";
        }
    }
    //---

    // STRING vector
    static void print(std::vector<std::string> vec)
    {
        std::cout << "\n-------\n";
        for (unsigned int i = 0; i < vec.size(); i++)
        {
            std::cout << " " << vec[i] << "\n";
        }
        std::cout << "-------\n";
    }
    //---
    static void print(std::vector<std::string> vec, std::string vecName)
    {
        std::cout << "\n-------\n";
        std::cout << vecName;
        std::cout << "\n-------\n";
        for (unsigned int i = 0; i < vec.size(); i++)
        {
            std::cout << " " << vec[i] << "\n";
        }
        std::cout << "-------\n";
    }
    //---
    static void printRow(std::vector<std::string> vec)
    {
        for (unsigned int i = 0; i < vec.size(); i++)
        {
            std::cout << vec[i] << " ";
        }
        std::cout << "\n";
    }
    //---
    static void printRow(std::vector<std::string> vec, std::string nameVec)
    {
        std::cout << nameVec << ": ";
        for (unsigned int i = 0; i < vec.size(); i++)
        {
            std::cout << vec[i] << " ";
        }
        std::cout << "\n";
    }
    //---
    static void printRowMatlab(std::vector<std::string> vec)
    {
        if (vec.size() == 0) std::cout << "[]\n";
        else
        {
            std::cout << "[";
            for (unsigned int i = 0; i < vec.size()-1; i++)
            {
                std::cout << vec[i] << ", ";
            }
            std::cout << vec[vec.size()-1] << "]\n";
        }
    }
    //---
    static void printRowMatlab(std::vector<std::string> vec, std::string nameVec)
    {
        std::cout << nameVec << ": ";
        if (vec.size() == 0) std::cout << "[]\n";
        else
        {
            std::cout << "[";
            for (unsigned int i = 0; i < vec.size()-1; i++)
            {
                std::cout << vec[i] << ", ";
            }
            std::cout << vec[vec.size()-1] << "]\n";
        }
    }
    //---
    static void printRowStd(std::vector<std::string> vec)
    {
        if (vec.size() == 0) std::cout << "{}\n";
        else
        {
            std::cout << "{";
            for (unsigned int i = 0; i < vec.size()-1; i++)
            {
                std::cout << vec[i] << ", ";
            }
            std::cout << vec[vec.size()-1] << "}\n";
        }
    }
    //---
    static void printRowStd(std::vector<std::string> vec, std::string nameVec)
    {
        std::cout << nameVec << ": ";
        if (vec.size() == 0) std::cout << "{}\n";
        else
        {
            std::cout << "{";
            for (unsigned int i = 0; i < vec.size()-1; i++)
            {
                std::cout << vec[i] << ", ";
            }
            std::cout << vec[vec.size()-1] << "}\n";
        }
    }
    //---
    static void printRowStdInLine(std::vector<std::string> vec)
    {
        if (vec.size() == 0) std::cout << "{}";
        else
        {
            std::cout << "{";
            for (unsigned int i = 0; i < vec.size()-1; i++)
            {
                std::cout << vec[i] << ", ";
            }
            std::cout << vec[vec.size()-1] << "}";
        }
    }
    //---
    static void printRowStdInLine(std::vector<std::string> vec, std::string nameVec)
    {
        std::cout << nameVec << ": ";
        if (vec.size() == 0) std::cout << "{}";
        else
        {
            std::cout << "{";
            for (unsigned int i = 0; i < vec.size()-1; i++)
            {
                std::cout << vec[i] << ", ";
            }
            std::cout << vec[vec.size()-1] << "}";
        }
    }
    //---
};