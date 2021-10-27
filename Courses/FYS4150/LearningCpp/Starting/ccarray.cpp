#include <iostream>
#include <cmath>

using namespace std;


int main(int argc, char const *argv[]) {

    //only have to include bounds for the first dimention
    int grid[3][2] = {
                    {1, 2},
                    {3, 4},
                    {5, 6}
                };

    cout << grid[2][0] << '\n';

    return 0;
}
