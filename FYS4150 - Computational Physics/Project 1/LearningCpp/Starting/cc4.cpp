#include <iostream>
#include <cmath>

using namespace std;

void sayHi(string name, int age=10) {
    cout << "Hello " << name << "!" << endl;
    // &&, ||, ! work
    if (age > 10 && name=="Karl") {
        cout << "Oldie" << endl;
    }
    cout << age << endl;
}

int main(int argc, char const *argv[]) {

    int luckyNums[5] = {40, 15, 22};

    cout << luckyNums[0] << endl;

    sayHi("Karl", 30);

    return 0;
}
