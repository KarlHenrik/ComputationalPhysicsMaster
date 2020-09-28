#include <iostream>
#include <cmath>

using namespace std;


int main(int argc, char const *argv[]) {

    int age = 19;
    double gpa = 2.7;
    string name = "Mike";
    int *pAge = &age;
    int **ppAge = &pAge;

    cout << *&ppAge << endl;
    cout << ppAge;

    return 0;
}
