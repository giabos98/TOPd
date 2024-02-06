#pragma once
#include "codeHeader.h"

class exception : public std::runtime_error {
    std::string msg;
public:
    exception(const std::string &arg, const char *file, int line) :
    std::runtime_error(arg) {
        std::ostringstream o;
        o << file << ":" << line << ": " << arg;
        msg = o.str();
    }
    exception(const char* message) :
    std::runtime_error(message) {
        std::ostringstream o;
        o << message;
        msg = o.str();
    }

    ~exception() throw() {}
    const char *what() const throw() {
        return msg.c_str();
    }
    //---
};
#define throw_line(arg) throw exception(arg, __FILE__, __LINE__);
#define checkPtr(arg) CHECK_ALLOC::checkPtrAllocation(arg, __FILE__, __LINE__);
#define pause() do {std::cout << '\n' << "Press the Enter key to continue.\n";} while (std::cin.get() != '\n');

#define prec double
#define format "f"
#define format_e "e"

class CHECK_ALLOC
{
    public:
    static void checkPtrAllocation(long double* p, const char *file, int line)
    {
        if (p == nullptr) throw exception("ERROR IN POINTER ALLOCATION", file, line);
    }
    static void checkPtrAllocation(prec* p, const char *file, int line)
    {
        if (p == nullptr) throw exception("ERROR IN POINTER ALLOCATION", file, line);
    }
    static void checkPtrAllocation(int* p, const char *file, int line)
    {
        if (p == nullptr) throw exception("ERROR IN POINTER ALLOCATION", file, line);
    }
    //------
    static void checkPtrAllocation(long double** p, const char *file, int line)
    {
        if (p == nullptr) throw exception("ERROR IN POINTER ALLOCATION", file, line);
    }
    static void checkPtrAllocation(double** p, const char *file, int line)
    {
        if (p == nullptr) throw exception("ERROR IN POINTER ALLOCATION", file, line);
    }
    static void checkPtrAllocation(int** p, const char *file, int line)
    {
        if (p == nullptr) throw exception("ERROR IN POINTER ALLOCATION", file, line);
    }
    //------
    static void checkPtrAllocation(long double*** p, const char *file, int line)
    {
        if (p == nullptr) throw exception("ERROR IN POINTER ALLOCATION", file, line);
    }
    static void checkPtrAllocation(double*** p, const char *file, int line)
    {
        if (p == nullptr) throw exception("ERROR IN POINTER ALLOCATION", file, line);
    }
    static void checkPtrAllocation(int*** p, const char *file, int line)
    {
        if (p == nullptr) throw exception("ERROR IN POINTER ALLOCATION", file, line);
    }
};