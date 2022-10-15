#include <iostream>

#include "io.hpp"
#include "RadModel.hpp"

using namespace std;

int main(int argc, char* argv[])
{
    string fileName = "./config/params";
    myLib::radParams params = myLib::readParams(fileName);
    myLib::RadModel radModel(params);

    // Gen spectrum here

    cout << "Complete." << endl;

    return 0;
}
