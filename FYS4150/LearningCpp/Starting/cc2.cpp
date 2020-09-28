//Input
#include <iostream>
#include <cmath>

using namespace std;

int main(int argc, char const *argv[]) {

    int age;

    cout << "Enter your age: ";
    cin >> age;
    cout << "You are " << age << " years old!" << endl;

    string name;

    cout << "Enter your name: ";
    cin.ignore(); //If you're using getline() after cin >> something,
                  //you need to flush the newline
                  //character out of the buffer in between.
                  //You can do it by using cin.ignore().
                  //https://stackoverflow.com/questions/18786575/using-getline-in-c
    getline(cin, name);
    cout << "You are " << name << endl;

    return 0;
}
