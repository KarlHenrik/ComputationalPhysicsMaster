// A comment line begins like this in C++ programs
// Standard ANSI-C++ include files
#include <cstdlib> // atof function
#include <iostream>   // input and output
#include <cmath>      // math functions

using namespace std;

int main (int argc, char* argv[])
{
  // convert the text argv[1] to double using atof:
  double r = atof(argv[1]);  // convert the text argv[1] to double
  double s = sin(r);
  cout << "Hello, World nice! sin(" << r << ") = " << s << endl;
  return 0;           // success execution of the program
}
