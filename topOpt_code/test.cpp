#include <iostream>
#include <vector>
#include "CODE_HEADERS/codeHeader.h"

class A
{
    public: 
    int a;
    A(){};
    A(int num)
    {
        a = num;
    }
    int get()
    {
        return a;
    }
    void print()
    {
        std::cout << "\nA print";
        std::cout << "\nA: " << a << "\n";
    }
};

class B : public A
{
    public: 
    int b;
    B(){};
    B(int num1,int num2) : A(num1)
    {
        b = num2;
    }
    int get()
    {
        return b;
    }
    void print()
    {
        A::print();
        std::cout << "B print";
        std::cout << "\nB: " << b << "\n";
    }
};

int main()
{
    int num1 = 3;
    int num2 = 5;
    A a(num1);
    B b(num1, num2);
    a.print();
    b.print();

    int n = 1;
    std::vector<std::shared_ptr<A>> vect;
    for (int i = 0; i < n; i++)
    {
        vect.emplace_back(std::make_shared<B>(b));
        (*vect[i]).print();
    }
    return 0;
}