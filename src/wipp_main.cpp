#include <stdio.h>
#include <math.h>
#include "time.h"
#include <unistd.h>
#include <stdlib.h>
#include <iostream>

#include <consts.h>
#include <wipp.h>

// #include <boost/program_options.hpp>
// #include <cassert>

using namespace std;

// namespace po = boost::program_options;


 
int main(int argc, char *argv[]) 
{
    // char *fileName;
    string fileName;

    if (argc != 2) {
        fileName = "python/four_adjacent.ray";
    } else {
        fileName = argv[1];
    }




    // fileName= "rayfile.ray";
    cout << "---- 3D WIPP ----\n";
    read_rayfile(fileName);
    return 0; // Return statement.
} // Closing Main.