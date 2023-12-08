#include <iostream>
#include <vector>
using namespace std;

#include "ProgressBar.h"
#include "CommandLine.h"

int main(int argc, char *argv[])
{
   CommandLine CL(argc, argv);

   string InputFileName = CL.Get("Input");
   string OutputFileName = CL.Get("Output");

   return 0;
}





