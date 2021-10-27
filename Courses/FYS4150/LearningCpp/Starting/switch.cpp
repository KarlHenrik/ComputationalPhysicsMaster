#include <iostream>
#include <cmath>

using namespace std;

string getDay(int dayNum) {
    string dayName;

    switch (dayNum) {
        case 0:
            dayName = "Mandag";
            break;
        case 1:
            dayName = "Tirsdag";
            break;
        case 2:
            dayName = "Onsdag";
            break;
        case 3:
            dayName = "Torsdag";
            break;
        case 4:
            dayName = "Fredag";
            break;
        default:
            dayName = "Helg!";
    }

    return dayName;
}

int main(int argc, char const *argv[]) {

    cout << getDay(3);

    return 0;
}
