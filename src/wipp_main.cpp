#include <stdio.h>
#include <math.h>
#include "time.h"
#include <unistd.h>
#include <stdlib.h>
#include <iostream>

#include <consts.h>
#include <wipp.h>

using namespace std;
 
int main(int argc, char *argv[]) 
{
    char *fileName;

    fileName= "rayfile.ray";
    cout << "---- 3D WIPP ----\n";
    read_rayfile(fileName);
    return 0; // Return statement.
} // Closing Main.