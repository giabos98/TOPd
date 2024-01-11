#include "CODE_HEADERS/codeHeader.h"
#include "nsProblem.h"
#include "TopOpt.h"
int PARALLEL::nThread = 4;
int main()
{    

    std::string inputFileNS = "INPUT_FILES/readProblemNS.txt";
    TOP_OPT topOpt(inputFileNS);

    //*-*-*--*-*-*-*-*-*-*
    topOpt.solve();
    //*-*-*-*-*-*-*-*-*-*-*-

    return 0;
}