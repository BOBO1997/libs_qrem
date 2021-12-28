#include <iostream>

using namespace std;

class A {
public:
    int a;
    A(){a = 1;};
    virtual void f(int i, ...) = 0;
};

class B: public A {
public:
    int b;
    B(): A() {b = 2;};
    void f(int i, ...) {}
    // void f(int i, int j, int k) { this->b = i + j + k; }
};

int main() {
    B* c = new B();
    c->f(1,2,3,4,5);
    cout << c->a << endl;
    cout << c->b << endl;
}