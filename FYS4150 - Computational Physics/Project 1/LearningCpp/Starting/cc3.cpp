#include <iostream>
#include <cmath>

using namespace std;

int main(int argc, char const *argv[]) {
    /*
    int num1, num2;
    cout << "Enter first number: ";
    cin >> num1;

    cout << "Enter second number: ";
    cin >> num2;

    cout << num1 + num2;
    */

    string coolLine;
    cout << "Give me: ";
    getline(cin, coolLine);
    std::cout << coolLine << '\n';

    return 0;
}
