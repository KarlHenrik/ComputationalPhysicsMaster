//Variables and math functions
#include <iostream>
#include <cmath>

using namespace std;

int main(int argc, char const *argv[]) {
    string name = "John";
    int age = 36.99;
    float myNum;
    myNum = age + 0.5;
    double better = 1;

    bool isMale = true;
    int myGod = name.find("uhn", 0);
    cout << name.length() << endl;
    cout << myGod << endl;
    cout << name[3] << endl;

    cout << pow(3, 2) << endl;
    cout << fmax(10, -22) << endl;
    std::cout << abs(-23) << '\n';
    // round(), ceil(), floor(), fmax(), fmin() ,abs()
    return 0;
}
