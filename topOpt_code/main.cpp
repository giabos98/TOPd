#include "CODE_HEADERS/codeHeader.h"
#include "nsProblem.h"
#include "TopOpt.h"
int PARALLEL::nThread = 4;
int main()
{    
    prec startTime = omp_get_wtime();

    // std::string inputFileNS = "INPUT_FILES/readProblemNS.txt";
    // TOP_OPT topOpt(inputFileNS);
    std::string inputFileDARCY = "INPUT_FILES/readProblemDarcy.txt";
    TOP_OPT topOpt(inputFileDARCY);
    

    //*-*-*--*-*-*-*-*-*-*
    topOpt.solve();
    //*-*-*-*-*-*-*-*-*-*-*-

    prec endTime = omp_get_wtime();

    prec totalTime = endTime - startTime;
    std::cout << "\n ---| Solution Time: " << totalTime << "\n";
    // topOpt.print_stats(totalTime);

    return 0;

}